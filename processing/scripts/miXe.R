######Script that applies the hmm to the samples processed with the SVM

set.seed(42)

#importing libraries
library(HMM)
library(optparse)
library(parallel)

options(warn=1)
######Defining default values
std_wd <-commandArgs()[1]
default_hmm_states <- 3
def_hmm_bw_max_iter <- 20
def_bw_delta <- 1E-9
def_train_chrom <- "chr1"

###### Parsing arguments

#get command line arguments, uses optparse, it's similar to Python's ArgParse
option_list = list(
  make_option(c("-w", "--work_directory"), type="character", default=getwd(), 
              help="This script's working directory. Please use the script's path. [Default = %default.]", metavar="character"),
  
  make_option(c("-D", "--dataset_directory"), type="character", default=NULL, 
              help="Path to sample folders (that contain each sample's TrainX-processed dataframes).", metavar="character"),
  
  make_option(c("-t", "--num_threads"), type = "numeric", default = 2,
              help="Number of threads to spawn to train/decode sequences with HMMs. [Default = %default.]", metavar="character"),
  
  make_option(c("-m", "--model_identifier"), type="character", default="svc", 
              help="Identifier for model (e.g. SVC or NN) [Default = %default. TODO: fixed to SVM, remove]", metavar="character"),
  
  make_option(c("-n", "--use_noise"), type="character", default=TRUE, 
              help="Load dataframes processed with miXer SVM trained on noisy XLR dataset [Default = %default. TODO: eliminare in qualche modo].", metavar="character"),
  
  make_option(c("-b", "--max_baum_welch_iterations"), type="numeric", default=def_hmm_bw_max_iter, 
              help="Maximum Baum-Welch algorithm iterations. [Default = %default].", metavar="character"),
  
  make_option(c("-d", "--baum_welch_delta"), type="numeric", default=def_bw_delta, 
              help="Minimum parameter delta for Baum-Welch algorithm. Stop if change is smaller than [Default = %default].", metavar="character"),

  make_option(c("-o", "--output_directory"), type="character", default=paste0(std_wd, "/HMM_processed_output/"), 
              help="Specify output directory. [Default = %default].", metavar="character"),

  make_option(c("-v", "--training_verbose"), type = "logical", default = FALSE, 
              help = "Whether to print training info (current chromosome processing, when it's making the windows). [Default = %default]", metavar = "character"),
  
  make_option(c("-r", "--pseudo_autosomal_regions"), type = "character", default = NULL,
              help = "Path to current reference's X chromosome's pseudo autosomal regions.", metavar = "character"),
  
  make_option(c("-y", "--resume_execution"), type = "logical", default = FALSE, 
              help = "Whether to resume partial execution.", metavar = "character"),
  make_option(c("-l", "--show_wt_windows"), type = "logical", default = FALSE,
                help = "Whether to show WT windows in window calls file. [Default = %default]", metavar = "character"),
  make_option(c("-z", "--script_resource_dir"), type = "character", default = dirname(std_wd),
              help = "Path pointing to script resources. [Default = %default]", metavar="character")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

################ Function definition (TODO: mettere in script separato)

get_true_false_string <- function(bool){
  #needed since trainX training script renames folders with booleans like this
  #Noise_True
  #while when using as.character(bool) in R returns "TRUE" or "FALSE"
  #I need this to recover the right model folder (parameters are written in the folder name)
  #TODO FIND ANOTHER SOLUTION
  if (bool == T){
    return("True")
  }
  else{
    return("False")
  }
}

guess_23 <- function(sample_pred, predictions_column_name, sample, output_dir, x_aliases, par_regions){

  #take only chrX target regions
  chrx <- subset(sample_pred, sample_pred$Chr %in% x_aliases)
  #take only non PAR regions
  non_par <- subset(chrx, (chrx$End < par_regions[1, ]$Start) | (chrx$Start > par_regions[1,]$End & chrx$End < par_regions[2,]$End))
  # Check if chrx or non_par is empty
  if (nrow(chrx) == 0 || nrow(non_par) == 0) {
    one_X_copy = TRUE
    is_problematic = TRUE
    no_X = TRUE
  } else {
    no_X = FALSE
    #get ML-model predictions
    ml_pred <- as.numeric(non_par[[predictions_column_name]])
    pred_occurrencies <- as.data.frame(table(ml_pred))
    #vote by majority to estimate number of X chromosomes
    pred_occurrencies_ordered <- pred_occurrencies[order(-pred_occurrencies$Freq), ]
    most_frequent_prediction <- as.character(pred_occurrencies_ordered[1,]$ml_pred)
    most_frequent_num_occurrencies <- pred_occurrencies_ordered[1,]$Freq
    total_occurrencies <- sum(pred_occurrencies_ordered$Freq)
    fraction_of_total <- most_frequent_num_occurrencies/total_occurrencies
    
    if (most_frequent_prediction == -1){
      one_X_copy <- TRUE
      is_problematic <- TRUE
    }
    else if (most_frequent_prediction == 0){
      one_X_copy <- FALSE
      is_problematic <- TRUE
    }
    else if (most_frequent_prediction == 1 | most_frequent_prediction == 2 ){
      one_X_copy <- FALSE
      is_problematic <- TRUE
    }
    else if (most_frequent_prediction == -2){
      one_X_copy <- FALSE
      is_problematic <- FALSE
    }
  }
  
  out <- cbind(as.data.frame(one_X_copy), as.data.frame(is_problematic), as.data.frame(no_X))
  colnames(out) <- c("one_X_copy", "is_problematic", "no_X_chr")
  
  #saving SVM prediction occurrencies
  l <- paste(sample, no_X, one_X_copy, most_frequent_prediction, most_frequent_num_occurrencies, total_occurrencies, fraction_of_total, sep = "\t" )
  txt_file_path <- file.path(output_dir, "SVM_guess_ploidy_results.txt")

  if (file.exists(txt_file_path)){
    cat(l, file= txt_file_path, sep = "\n", append=TRUE)
    
  }
  else{
    cat(paste0("Sample\t", "NO_X_CHR\t", "One_X_copy\t", "Most_frequent_CopyNumb_prediction\t", "Occurrencies\t", "Total_predictions_made\t", "MPP_fraction_of_total_occurrencies", "\n"), file = txt_file_path, sep = "\t" )
    cat(l, file= txt_file_path, sep="\n", append = TRUE)
  }
  
  return(out)
}
print_window_number <- function(windows){
  windows_ddels <- subset(windows, windows$Call == "-2")
  windows_dels <- subset(windows, windows$Call == "-1")
  windows_dups <- subset(windows, windows$Call == "+1")
  windows_ddups <- subset(windows, windows$Call == "2+")
  print(paste("Windows by CNV type: DDEL(-2):", nrow(windows_ddels), " DEL(-1): ", nrow(windows_dels), " DUP(+1): ", nrow(windows_dups), " MDUP(2+): ", nrow(windows_ddups)))
  return()
}
#########


################ Processing part of the script
#### Setting everything up
#check if no arguments have been passed
# if (is.null(opt$file)){
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }

#check if working directory exists, if it doesn't, create it
if (!dir.exists(opt$work_directory)){
  dir.create(opt$work_directory)
}

#changing working directory to that of the script
setwd(opt$work_directory)

#importing self defined functions
#find the resource folder
if (opt$script_resource_dir != std_wd){
  resource_folder = opt$script_resource_dir
}else{
  resource_folder = file.path(getwd(), "HMM_resources")
}

#get the function files
source(file.path(resource_folder, "get_model_id.R"))
source(file.path(resource_folder, "train_hmm.R"))
source(file.path(resource_folder, "get_hmm_states_by_chromosome.R"))
source(file.path(resource_folder, "add_concordant_calls_column.R"))
source(file.path(resource_folder, "apply_hmm_single_chrom.R"))
source(file.path(resource_folder, "get_untrained_hmm.R"))
source(file.path(resource_folder, "get_untrained_hmm_MChrX.R"))
source(file.path(resource_folder, "window_maker_X3.R"))


#################
#TODO SVM mappability filter parameters: keep? set as arguments? who knows
prefilter_low_mapp=TRUE
mapp_thr=0.6
#################

#define main output dir
main_output_dir <- opt$output_directory

if (!dir.exists(main_output_dir)){
  dir.create(main_output_dir)
}

#get correct model id (just a precaution)
model_id <- get_model_id(opt$model_identifier) #TODO remove, it is always an SVM now

#Scan SVM-processed dataset folder to find the list of samples (one sample for each subfolder)
trainx_processed_samples <- list.dirs(opt$dataset_directory, full.names = F, recursive = F)

#Actually load the datasetss 
##########  TODO set up a function for this
samples_to_process = list()
for (sample in trainx_processed_samples){
  sample_path <- file.path(opt$dataset_directory, sample)

  subfolders <- list.dirs(sample_path, full.names = F, recursive = F)
  
  pattern <- paste0("^", paste0(model_id, "*"))
  subfolders <- subfolders[grepl(pattern, subfolders)]
  pattern <- paste0("Noise_", get_true_false_string(as.logical(opt$use_noise)))
  subfolders <- subfolders[grepl(pattern, subfolders)]
  subfolder = file.path(sample_path, as.character(subfolders[grepl(pattern, subfolders)]) )
  l <- list()
  l[[as.character(sample)]] = subfolder
  samples_to_process <- append(samples_to_process, l)
}

#### Preparing for hmm training
#set max number of iterations for the baum-welch algorithm
hmm_bw_max_iter <- opt$max_baum_welch_iterations

#set baum-welch algorithm's delta
#if parameters change less than the defined delta, halt the hmm training algorithm
bw_delta <- opt$baum_welch_delta

#set hmm training chromosome(s)
train_chrom <- opt$training_chromosome

#### Processing the samples with the HMM(s)
iteration_count <- 1
savecols <- c("Chr", "Start", "End")

x_aliases <- c("chrX", "ChrX", "chrx", "X", "x")

#load pseudo-autosomal region coordinates file
par_regions <- read.table(opt$pseudo_autosomal_regions)
colnames(par_regions) <- c("Chr", "Start", "End")
#keep only chrx regions
par_regions <- subset(par_regions, par_regions$Chr %in% x_aliases)

svm_ploidy_dir <- main_output_dir

if (opt$resume_execution == T){
  #Get the list of the names of the files in the specified folder
  already_processed <- list.files(main_output_dir, full.names = FALSE)
  #Remove the names of the subfolders from the list of the filenames
  #TODO check if last file is completed
  trainx_processed_samples <- setdiff(trainx_processed_samples, already_processed)
}
tot_samples <- length(trainx_processed_samples)

for(sample in trainx_processed_samples){
  print(paste("Processing sample", paste(sample, paste("| Sample", iteration_count, "out of", tot_samples))))
  
  #get dataframes path
  data_path = samples_to_process[[sample]]
  
  pat <- "*_pred.txt*"
  filename <- sample
  
  #must save the processed dataframes with the script's result
  result_df_save_folder <- file.path(main_output_dir, paste0(sample))

  if (!dir.exists(result_df_save_folder)){
    dir.create(result_df_save_folder)
  }
  
  ##hmm matrices output folder
  #out_hmm <- file.path(result_df_save_folder, "trained_hmm_matrices")

  #if not present, make directory
  #if (!dir.exists(out_hmm)){
  #  dir.create(out_hmm)
  #}
  
  #there will be only one dataframe name that will match the pattern
  dataframe_name <- list.files(data_path, pattern=pat, all.files=FALSE,full.names=FALSE)
  dataframe_path <- file.path(data_path, dataframe_name)
  
  #load samples predictions
  sample_pred <- read.table(file = dataframe_path, sep="\t",quote="\"",fill=T,header=T)
  all_chrs <- as.character(unique(sample_pred$Chr))
  dataframe_colnames <- colnames(sample_pred)
  
  #normalize again NRC_poolNorm with median just to be safe
  #filter the dataframe to include only rows where "Chr" is not the X chromosome
  filtered_df <- sample_pred[!sample_pred$Chr %in%  c("chrx", "chrX", "ChrX", "Chrx", "X", "x"), ]
  
  #calculate the median of the "NRC_poolNorm" column for the filtered rows
  median_value <- median(filtered_df$NRC_poolNorm, na.rm = TRUE)

  #apply median normalization to the "NRC_poolNorm" column for the entire dataframe
  sample_pred$NRC_poolNorm <- sample_pred$NRC_poolNorm - median_value

  predictions_column_name <- dataframe_colnames[grepl("*_pred$", dataframe_colnames)]
  
  #print(paste("Training the HMM with max.", paste(hmm_bw_max_iter, paste("Baum-Welch algorithm iterations and", paste(bw_delta, "minimum parameter update value") ) )))
  
  #result_df <- data.frame()
  
  #print(paste("Training one HMM for each chromosome with max.", paste(hmm_bw_max_iter, paste("Baum-Welch algorithm iterations and", paste(bw_delta, "minimum parameter update value") ) )))
  
  #get the unique values of observations that the model made
  u_obs <- unique(sample_pred[[predictions_column_name]])
  
  #estimating if M or F using the SVM
  #TODO if sex is specified, override this
  estim_male <- guess_23(sample_pred, predictions_column_name, sample, svm_ploidy_dir, x_aliases, par_regions)
  is_male = estim_male[1,]$one_X_copy
  #set the number of cores you want to use for parallel execution
  num_cores <- opt$num_cores
  
  #define posterior probability column names
  post_prob_col_names <- c("DDEL_post_prob", "DEL_post_prob", "WT_post_prob", "DUP_post_prob")
  single_chroms <- split(sample_pred, sample_pred$Chr)
  
  ##### REMOVE COMMENT TO DEBUG WITH TWO CHROMOSOMES
  #dbg <- droplevels(sample_pred[sample_pred$Chr %in% c("chr16"), ]) # c("chr16","chrX")
  #single_chroms <- split(dbg,dbg$Chr)
  
  # data_list = list of objects to which apply the function in parallel
  # function_doing_something = function to apply to each element of data_list
  cores <- opt$num_threads # number of cores
  if (.Platform$OS.type == "unix") {
    results <- parallel::mclapply(1:length(single_chroms),# *
                                  function(idx) {
                                    res <- apply_hmm_single_chrom(single_chroms[[idx]], hmm_bw_max_iter, bw_delta,
                                                                  predictions_column_name,
                                                                  u_obs, is_male, prefilter_low_mapp, mapp_thr)
                                    return(res)
                                  },
                                  mc.cores = cores)
  } else{
    # when working in Windows you need to import all the data and libraries you need for your function
    cl <- parallel::makeCluster(cores)
    clusterEvalQ(cl, {
      # import here the libraries
      library(HMM)
    })
    clusterExport(cl,
                  c(
                    # import here (as strings) the variable name of the objects you are working with
                    "apply_hmm_single_chrom", "get_model_id", "train_hmm", "get_hmm_states_by_chromosome",
                    "add_concordant_calls_column", "get_untrained_hmm", "get_untrained_hmm_MChrX",
                    "get_hmm_states", "get_hmm_states_MChrX", "get_hmm_emission_symbols",
                    "get_hmm_emission_symbols_MChrX", "get_hmm_starting_probabilities", "get_hmm_starting_probabilities_MChrX",
                    "get_hmm_emission_probabilities", "get_hmm_emission_probabilities_MChrX",
                    "get_hmm_transition_probabilities", "get_hmm_transition_probabilities_MChrX",

                    "single_chroms", "hmm_bw_max_iter", "bw_delta", "predictions_column_name",
                    "hmm_matrices_output_folder", "u_obs", "is_male", "prefilter_low_mapp", "mapp_thr"
                  ),
                  envir = environment())
    results <- parLapply(cl,
                         1:length(single_chroms), # *
                         function(idx) {
                           res <- apply_hmm_single_chrom(single_chroms[[idx]], hmm_bw_max_iter, bw_delta,
                                                         predictions_column_name, hmm_matrices_output_folder,
                                                         u_obs, is_male, prefilter_low_mapp, mapp_thr)
                           return(res)
                         })
    stopCluster(cl)
  }
  result_df <- do.call(rbind, results)
  
  list2env(result_df, envir = environment())
  
  #write.csv(result_df, file = "result_df.csv", row.names = FALSE)
  
  #save all predictions and probability results
  savecols <- c(savecols, c(predictions_column_name, "hmm_state", "Mappability"))
  savecols <- c(savecols, post_prob_col_names)
  
  #make windows where the hmm state is constant
  print("Making windows")
  windows <- window_maker(result_df, estim_male[1,]$one_X_copy, predictions_column_name, par_regions, show_wt_windows = opt$show_wt_windows)
  
  #print("Saving CSVs containing both singleTR + HMM predictions.")
  all_results_filepath <- file.path(result_df_save_folder, paste0(sample, paste0(paste0(paste0("_", model_id), "_hmm_bw"), paste0(hmm_bw_max_iter,"_processed.tsv"))) )
  
  result_df <- result_df[savecols]
  
  #write.table(result_df, file = all_results_filepath, sep = "\t", row.names = FALSE)
  
  print(paste("Sample", paste(sample, paste("| ", paste(nrow(windows), "identified CNV Windows")))))
  print_window_number(windows)
  
  #filtering for window confidence (mean_post_prob*mean_mapp)
  hc09Windows <- subset(windows, windows$ProbCall > 0.9)
  print(paste("Sample", paste(sample, paste("| ", paste(nrow(hc09Windows), " identified PASS CNV\n ")))))
  print_window_number(hc09Windows)
  
  #save HMM windows  
  filepath_allwindows <- file.path(result_df_save_folder, paste0(filename, paste0("_hmm_bw", paste0(hmm_bw_max_iter, "_windows.bed"))))
  filepath_hc09windows <-file.path(result_df_save_folder, paste0(filename, paste0("_hmm_bw", paste0(hmm_bw_max_iter, "_PASS_ONLY_windows.bed"))))
  
  #print("Saving windows")
  write.table(windows, file=filepath_allwindows, quote = FALSE, sep = "\t", row.names = FALSE)
  
  #print("Saving high confidence windows")
  write.table(hc09Windows, file= filepath_hc09windows, quote = FALSE, sep = "\t", row.names = FALSE)
  
  iteration_count <- iteration_count + 1
  
}





















