########### Function that:
#1) Trains the HMM parameters on the current ML model-processed target regions (train_hmm)
#2) Saves the resulting transition and emission probabilities of the HMM (TODO: if already present, skip training and reload unless otherwise specified)
#3) Gets the most probable sequence of HMM states given the current ML model observation and the estimated HMM (get_hmm_states_by_chromosome)
#4) Adds a concordant ML model calls- HMM state column (add_concordant_calls_column)
#5) Returns the sequence dataframe

#' Apply an HMM to a single chromosome for prediction and analysis.
#'
#' This function applies the Hidden Markov Model (HMM) to a single chromosome for prediction and analysis. It trains
#' the HMM on the specified chromosome using the Baum-Welch algorithm and ML caller's predictions as observations
#' of the hidden states. It saves the trained HMM transition and emission matrices, obtains the most probable
#' sequence of states and states posterior probabilities, and performs additional analysis.
#'
#' @param curr_data A data frame containing the sample predictions.
#' @param base_hmm An initial HMM object to use as a starting point for training.
#' @param hmm_bw_max_iter The maximum number of Baum-Welch iterations to perform during training.
#' @param bw_delta The convergence threshold for the Baum-Welch algorithm.
#' @param chr The chromosome on which to train the HMM.
#' @param predictions_column_name The name of the column in the data frame `curr_data` that contains the ML caller's predictions.
#' @param post_prob_col_names A character vector containing the column names for the posterior probabilities.
#' @param unique_observations A character vector containing the unique observations for the HMM.
#' @param verbose A logical value indicating whether to print verbose output.
#'
#' @return A data frame representing the results of applying the HMM to the chromosome.
#'
#'

#find the resource folder
resource_folder = file.path(
  dirname(sub("--file=", "", commandArgs(trailingOnly = FALSE)[4])),
  "HMM_resources")
source(normalizePath(file.path(resource_folder, "get_untrained_hmm.R")))
source(normalizePath(file.path(resource_folder, "get_untrained_hmm_MChrX.R")))

apply_hmm_single_chrom <- function(curr_data, hmm_bw_max_iter, bw_delta,
                                   predictions_column_name,
                                   unique_observations, is_male, prefilter_low_mapp, mapp_thr,
                                   verbose = FALSE){
  
  curr_data <- droplevels(curr_data)
  curr_data_copy <- curr_data
  if (prefilter_low_mapp == TRUE) {
    #print("Prefiltering SVM observations: removing low mappability positions for HMM training and inference")
    curr_data <- subset(curr_data, Mappability > mapp_thr) # Keep only rows with good quality
  }
  chr <- as.character(unique(curr_data$Chr))
  if (verbose){
    print(paste("Training the HMM on", chr) )
  }
  #instantiate base hmm object
  if ( (chr %in% c("chrX", "chrx", "ChrX", "Chrx", "X", "x")) && (is_male == TRUE)){
    base_hmm <- get_untrained_hmm_MChrX(unique_observations)
  } else{
    base_hmm <- get_untrained_hmm(unique_observations)
  }
  trained_hmm <- train_hmm(curr_data, base_hmm, hmm_bw_max_iter, bw_delta, chr, predictions_column_name)
  #get sequence of states and states posterior probabilities
  if (verbose){
    print(paste("Getting most probable sequence of states for chromosome", chr))
  }
  result_df <- get_hmm_states_by_chromosome(curr_data, trained_hmm, predictions_column_name, is_male)
  #add concordant (HMM state - ML caller) column
  result_df <- add_concordant_calls_column(result_df, predictions_column_name, unique_observations)
  #casting factor columns to character (no idea why factors pop out for some columns, TODO check)
  i <- sapply(result_df, is.factor)
  result_df[i] <- lapply(result_df[i], as.character)
  i <- sapply(curr_data_copy, is.factor)
  curr_data_copy[i] <- lapply(curr_data_copy[i], as.character)
  # Ensure the "Chr" column has consistent types in both data frames
  curr_data_copy$Chr <- as.character(curr_data_copy$Chr)
  result_df$Chr <- as.character(result_df$Chr)
  # Find rows in curr_data_copy that are not in result_df based on Chr, Start, and End
  missing_rows <- curr_data_copy[
    !paste(curr_data_copy$Chr, curr_data_copy$Start, curr_data_copy$End) %in%
      paste(result_df$Chr, result_df$Start, result_df$End),
  ]
  # Initialize an NA-filled data frame for missing rows with the structure of result_df
  additional_columns <- setdiff(names(result_df), names(curr_data_copy))  # Find extra columns in result_df
  if (nrow(missing_rows) > 0 && length(additional_columns) > 0) {
    missing_rows_na <- missing_rows  # Start with the full structure of curr_data_copy
    for (col in additional_columns) {
      missing_rows_na[[col]] <- rep(NA, nrow(missing_rows_na))
    }
    # Combine result_df with the NA-filled missing rows
    result_df <- rbind(result_df, missing_rows_na)
  }
  #missing_rows_na[additional_columns] <- NA  # Add the extra columns with NA values
  # Sort by Chr, Start, and End
  result_df <- result_df[order(result_df$Chr, result_df$Start, result_df$End), ]
  return(result_df)
}