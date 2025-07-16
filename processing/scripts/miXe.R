set.seed(42)

suppressPackageStartupMessages({
  library(HMM)
  library(optparse)
  library(parallel)
  library(jsonlite)
})

std_wd <- commandArgs()[1]
def_hmm_bw_max_iter <- 20
def_bw_delta <- 1E-9

option_list = list(
  make_option(c("-j", "--json"), type="character", help="Path to the miXer json file"),
  make_option(c("-w", "--work_directory"), type="character", default=file.path(dirname(sub("--file=", "", commandArgs(trailingOnly = FALSE)[4])), "HMM_resources"), help="This script's working directory."),
  make_option(c("-D", "--dataset_directory"), type="character", default=NULL, help="Path to sample folders."),
  make_option(c("-o", "--output_directory"), type="character", default=file.path(getwd(), "HMM_processed_output"), help="Specify output directory."),
  # make_option(c("-z", "--script_resource_dir"), type="character", default=file.path(getwd(), 'HMM_resources'), help="Path pointing to script resources."),
  make_option(c("-n", "--use_noise"), type="logical", default=TRUE, help="Load dataframes processed with miXer SVM trained on noisy XLR dataset."),
  make_option(c("-b", "--max_baum_welch_iterations"), type="numeric", default=def_hmm_bw_max_iter, help="Maximum Baum-Welch algorithm iterations."),
  make_option(c("-d", "--baum_welch_delta"), type="numeric", default=def_bw_delta, help="Minimum parameter delta for Baum-Welch algorithm."),
  make_option(c("-y", "--resume_execution"), type="logical", default=TRUE, help="Whether to resume partial execution."),
  make_option(c("-l", "--show_wt_windows"), type="logical", default=FALSE, help="Whether to show WT windows in window calls file.")
)

opt <- parse_args(OptionParser(option_list=option_list))
source(file.path(opt$work_directory, "get_model_id.R"))
source(file.path(opt$work_directory, "train_hmm.R"))
source(file.path(opt$work_directory, "get_hmm_states_by_chromosome.R"))
source(file.path(opt$work_directory, "add_concordant_calls_column.R"))
source(file.path(opt$work_directory, "apply_hmm_single_chrom.R"))
source(file.path(opt$work_directory, "get_untrained_hmm.R"))
source(file.path(opt$work_directory, "get_untrained_hmm_MChrX.R"))
source(file.path(opt$work_directory, "window_maker_X3.R"))

get_py_bool_string <- function(bool_val) {
  if (isTRUE(bool_val)) "True" else "False"
}

guess_23 <- function(sample_pred, pred_col, sample_name, output_dir, x_aliases, par_regions) {
  chrx <- subset(sample_pred, Chr %in% x_aliases)
  non_par <- subset(chrx, (End < par_regions[1, ]$Start) | (Start > par_regions[1, ]$End & End < par_regions[2, ]$End))

  no_X <- nrow(non_par) == 0
  one_X_copy <- FALSE
  is_problematic <- FALSE
  most_frequent_prediction <- NA
  most_frequent_num_occurrencies <- NA
  total_occurrencies <- NA
  fraction_of_total <- NA

  if (!no_X) {
    ml_pred <- as.numeric(non_par[[pred_col]])
    pred_counts <- as.data.frame(table(ml_pred))
    pred_counts_ordered <- pred_counts[order(-pred_counts$Freq), ]

    most_frequent_prediction <- as.character(pred_counts_ordered[1,]$ml_pred)
    most_frequent_num_occurrencies <- pred_counts_ordered[1,]$Freq
    total_occurrencies <- sum(pred_counts_ordered$Freq)
    fraction_of_total <- most_frequent_num_occurrencies / total_occurrencies
    
    if (most_frequent_prediction == -2) {
      one_X_copy <- TRUE
      is_problematic <- FALSE
    } else {
      one_X_copy <- FALSE
      is_problematic <- TRUE
    }
  }

  summary_line <- paste(sample_name, no_X, one_X_copy, most_frequent_prediction, most_frequent_num_occurrencies, total_occurrencies, fraction_of_total, sep = "\t")
  txt_file_path <- file.path(output_dir, "SVM_guess_ploidy_results.txt")

  if (!file.exists(txt_file_path)) {
    header <- "Sample\tNO_X_CHR\tOne_X_copy\tMost_frequent_CopyNumb_prediction\tOccurrencies\tTotal_predictions_made\tMPP_fraction_of_total_occurrencies\n"
    cat(header, file = txt_file_path)
  }
  cat(summary_line, file = txt_file_path, sep = "\n", append = TRUE)

  return(list(one_X_copy = one_X_copy, is_problematic = is_problematic))
}

`%||%` <- function(x, y) if (is.null(x)) y else x
print_window_number <- function(windows, sample_name) {
  counts <- table(windows$Call)
  cat(sprintf("Sample %s CNV Windows | DDEL(-2): %d | DEL(-1): %d | DUP(+1): %d | MDUP(2+): %d\n",
              sample_name,
              counts["-2"] %||% 0, counts["-1"] %||% 0,
              counts["+1"] %||% 0, counts["2+"] %||% 0))
}

process_sample <- function(sample_name, sample_info, config) {
  
  cat("--------------------------------------------------\n")
  cat(sprintf("Processing sample: %s\n", sample_name))
  
  result_df_save_folder <- file.path(config$output_directory, sample_name)
  dir.create(result_df_save_folder, showWarnings = FALSE, recursive = TRUE)
  
  dataframe_path <- list.files(sample_info$path, pattern = "*_pred.txt*", full.names = TRUE)[1]
  
  sample_pred <- read.table(dataframe_path, sep="\t", header=TRUE, fill=TRUE, quote="\"")
  pred_col_name <- grep("_pred$", colnames(sample_pred), value = TRUE)
  
  autosomes <- sample_pred[!sample_pred$Chr %in% config$x_aliases, ]
  median_nrc <- median(autosomes$NRC_poolNorm, na.rm = TRUE)
  sample_pred$NRC_poolNorm <- sample_pred$NRC_poolNorm - median_nrc
  
  ploidy_guess <- guess_23(sample_pred, pred_col_name, sample_name, config$output_directory, config$x_aliases, config$par_regions)
  is_male <- ploidy_guess$one_X_copy
  u_obs <- unique(sample_pred[[pred_col_name]])
  
  single_chroms <- split(sample_pred, sample_pred$Chr)
  
  cores <- config$json_data$threads
  if (.Platform$OS.type == "unix") {
    results <- mclapply(single_chroms, function(chrom_df) {
        tryCatch({
            apply_hmm_single_chrom(
              chrom_df, 
              hmm_bw_max_iter = opt$max_baum_welch_iterations, 
              bw_delta = opt$baum_welch_delta,
              predictions_column_name = pred_col_name,
              unique_observations = u_obs, 
              is_male = is_male,
              prefilter_low_mapp=TRUE, 
              mapp_thr=0.6
            )
          }, error = function(e) {
            message("Error in chromosome: ", e$message)
            return(NULL)
          })
        }, mc.cores = cores)
  } else {
    cl <- makeCluster(cores)
    clusterEvalQ(cl, {
      library(HMM)
    })
    exported_vars <- c(
      "apply_hmm_single_chrom", "get_model_id", "train_hmm", "get_hmm_states_by_chromosome",
      "add_concordant_calls_column", "get_untrained_hmm", "get_untrained_hmm_MChrX",
      "get_hmm_states", "get_hmm_states_MChrX", "get_hmm_emission_symbols",
      "get_hmm_emission_symbols_MChrX", "get_hmm_starting_probabilities", "get_hmm_starting_probabilities_MChrX",
      "get_hmm_emission_probabilities", "get_hmm_emission_probabilities_MChrX",
      "get_hmm_transition_probabilities", "get_hmm_transition_probabilities_MChrX",
      "config", "pred_col_name", "u_obs", "is_male"
    )
    clusterExport(cl, varlist = exported_vars, envir = environment())
    
    results <- parLapply(cl, single_chroms, function(chrom_df) {
      apply_hmm_single_chrom(
        chrom_df, 
        hmm_bw_max_iter = config$max_baum_welch_iterations, 
        bw_delta = config$baum_welch_delta,
        predictions_column_name = pred_col_name,
        u_obs = u_obs, 
        is_male = is_male,
        prefilter_low_mapp=TRUE, 
        mapp_thr=0.6
      )
    })
    stopCluster(cl)
  }
  
  result_df <- do.call(rbind, results)
  
  cat("Making windows from HMM results...\n")
  windows <- window_maker(result_df, is_male, pred_col_name, config$par_regions, show_wt_windows = config$show_wt_windows)
  hc_windows <- subset(windows, ProbCall > 0.9)
  
  cat(sprintf("Sample %s | Identified %d total CNV Windows\n", sample_name, nrow(windows)))
  print_window_number(windows, sample_name)
  cat(sprintf("Sample %s | Identified %d PASS (Prob > 0.9) CNV Windows\n", sample_name, nrow(hc_windows)))
  print_window_number(hc_windows, sample_name)
  
  filepath_allwindows <- file.path(result_df_save_folder, sprintf("%s_hmm_bw%d_windows.bed", sample_name, config$max_baum_welch_iterations))
  filepath_hc_windows <- file.path(result_df_save_folder, sprintf("%s_hmm_bw%d_PASS_ONLY_windows.bed", sample_name, config$max_baum_welch_iterations))
  
  write.table(windows, file=filepath_allwindows, quote=FALSE, sep="\t", row.names=FALSE)
  write.table(hc_windows, file=filepath_hc_windows, quote=FALSE, sep="\t", row.names=FALSE)
  
  return(TRUE)
}

json_data <- fromJSON(opt$json)
par_regions <- read.table(json_data$par, header=FALSE, col.names=c("Chr", "Start", "End"))
x_aliases <- c("chrX", "ChrX", "chrx", "X", "x")

dir.create(opt$output_directory, showWarnings = FALSE, recursive = TRUE)

model_id <- "SVC"
all_samples <- list.dirs(opt$dataset_directory, full.names = TRUE, recursive = FALSE)

if (isTRUE(opt$resume_execution)) {
  already_processed <- list.files(opt$output_directory, full.names = FALSE)
  all_samples <- setdiff(all_samples, already_processed)
}
samples_to_process <- list()
for (sample_path in all_samples) {
    sample <- basename(sample_path)
    noise_str <- get_py_bool_string(opt$use_noise)
    print(noise_str)
    pattern <- paste0("^", model_id, ".*Noise_", noise_str)
    subfolder_path <- list.files(sample_path, 
                                 pattern = pattern, 
                                 full.names = TRUE, 
                                 include.dirs = TRUE)
    if (length(subfolder_path) == 1) {
        samples_to_process[[sample]] <- list(path = subfolder_path[1])
    }
}
cat(sprintf("Found %d samples to process.\n", length(samples_to_process)))

config <- c(opt, list(json_data = json_data, par_regions = par_regions, x_aliases = x_aliases))

if (length(samples_to_process) > 0) {
    for (i in seq_along(samples_to_process)) {
        sample_name <- names(samples_to_process)[i]
        sample_info <- samples_to_process[[i]]
        process_sample(sample_name, sample_info, config)
    }
}

cat("--------------------------------------------------\n")
cat("HMM processing complete.\n")