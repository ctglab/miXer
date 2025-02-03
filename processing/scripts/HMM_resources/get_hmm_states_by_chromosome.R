########### Function that uses the trained HMM on the exome sequence and obtains the most probable sequence of HMM states

#' Apply the Viterbi algorithm to a sequence of observations and estimate the most probable
#' sequence of states using a trained Hidden Markov Model (HMM).
#'
#' This function takes a trained HMM and applies the Viterbi algorithm to obtain the most probable
#' sequence of states that explains the current observation sequence. It also estimates the posterior
#' probability of the state sequence.
#'
#' @param y A data frame containing the observations and chromosome information.
#' @param trained_hmm An object representing the trained HMM.
#' @param predictions_column_name The name of the column in the data frame `y` that contains the predictions.
#'
#' @return A data frame containing the original observations (`y`), the most probable sequence of states,
#'   and the corresponding posterior probabilities.
#'
get_hmm_states_by_chromosome <- function(y, trained_hmm, predictions_column_name, is_male) {
  
  #function that takes the trained hmm and applies the viterbi algorithm to obtain the most probable
  #sequence of states that explains the current observation sequence
  #and also estimates the posterior probability of the state sequence
  all_states <- data.frame()
  all_post_probs <- data.frame()
  for (chr in unique(y$Chr)) {
    #get current chromosome positions
    curr_chrom_pos <- subset(y, y$Chr == chr)
    curr_chrom_pos <- curr_chrom_pos[order(curr_chrom_pos$Start), ]
    #calculate distance between exons
    curr_chrom_pos$Position <- (curr_chrom_pos$End - curr_chrom_pos$Start)*0.5
    curr_chrom_pos$Position <- curr_chrom_pos$Start + curr_chrom_pos$Position
    curr_chrom_pos$Dist_from_next <- c(0, diff(curr_chrom_pos$Position))
    #get exons which are far from their next one
    id <- which(curr_chrom_pos$Dist_from_next > 0.25*10^6, arr.ind = TRUE)
    id <- c(id, nrow(curr_chrom_pos))
    #check if there are two adjacent ids (e.g. 9430 and 9431)
    #if this happens, merge the first with the second
    id2 = c()
    for (i in 1:(length(id) -1)){
      curr_id = id[i]
      next_id = id[i+1]
      id2 = c(id2, curr_id)
      if (curr_id == next_id | curr_id == (next_id -1) ){
        id2 <- id2[-length(id2)]
      }
    }
    id2 <- c(id2, nrow(curr_chrom_pos))
    
    #split in chunks if distance between exons is greater than a threshold
    for (i in 1:length(id2)){
      if(i==1){
        subseq <- curr_chrom_pos[ 1:(id2[i]-1), ]
      }else if (i < length(id2)){
        subseq <- curr_chrom_pos[ (id2[i-1]):(id2[i]-1), ]
      }
      else{
        subseq <- curr_chrom_pos[(id2[i-1]):(id2[i]), ]
      }
      #get current chromosome's miXer predictions
      curr_observations <- as.character(subseq[[predictions_column_name]])
      curr_observations <- append(curr_observations, "0", 1)
      
      #estimate the most probable sequence of states using the Viterbi algorithm
      mpss <- as.data.frame(viterbi(trained_hmm$hmm, curr_observations))
      #estimate each state's posterior probability
      post_probs <- as.data.frame(t(posterior(trained_hmm$hmm, curr_observations)))
      mpss = as.data.frame(mpss[-1,])
      post_probs = as.data.frame(post_probs[-1,])
      colnames(mpss) <- c("hmm_state")
      if (is_male == FALSE || !(chr %in% c("chrX", "chrx", "ChrX", "Chrx", "X", "x")) ){
        #3-state HMM has no DDEL_post_prob_column, must set to zero
        post_probs <- cbind(a = 0, post_probs)
      }
      colnames(post_probs) <- c("DDEL_post_prob", "DEL_post_prob", "WT_post_prob", "DUP_post_prob")
      all_states <- rbind(all_states, mpss)
      all_post_probs <- rbind(all_post_probs, post_probs)
    }
  }
  result_df <- cbind(y, all_states)
  result_df <- cbind(result_df, all_post_probs)
  return(result_df)
}