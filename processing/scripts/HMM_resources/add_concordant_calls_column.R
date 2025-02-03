########### Function that adds a column to the processed dataframe containing concordant SVM-HMM single target region calls

add_concordant_calls_column <- function(result_df, predictions_column_name,
                                        unique_observations) {
  
  result_df$concordant <- NA
  del_states <- c("DEL")
  wt_states <- c("WT")
  dup_states <- c("DUP")
  
  if ("-2" %in% unique_observations){
    doubledel_states <- c("DDEL")
    result_df$concordant[(result_df$hmm_state  %in% doubledel_states & result_df[[predictions_column_name]] == -2)] <- "DDEL"
  }
  
  #if concordant, keep the prediction
  result_df$concordant[(result_df$hmm_state  %in% del_states & result_df[[predictions_column_name]] == -1)] <- "DEL"
  result_df$concordant[(result_df$hmm_state %in% dup_states & result_df[[predictions_column_name]]  == 1)] <- "DUP"
  result_df$concordant[(result_df$hmm_state %in% wt_states & result_df[[predictions_column_name]]  == 0)] <- "WT"
  
  #if discordant, call WT
  result_df$concordant <- as.character(result_df$concordant)
  result_df$concordant[is.na(result_df$concordant)] <- "WT"
  
  return(result_df)
}