########### Function that trains an HMM on a chromosome and returns the trained HMM

#' Train a first-order Hidden Markov Model (HMM) on a selected chromosome using the Baum-Welch algorithm.
#'
#' This function trains a first-order HMM on the specified chromosome using the Baum-Welch algorithm. The
#' HMM is trained based on the Maximum Likelihood (ML) caller's predictions as observations of the hidden states.
#'
#' @param y A data frame containing the observations and chromosome information.
#' @param base_hmm An initial HMM object to use as a starting point for training.
#' @param hmm_bw_max_iter The maximum number of Baum-Welch iterations to perform during training.
#' @param bw_delta The convergence threshold for the Baum-Welch algorithm.
#' @param train_chrom The chromosome on which to train the HMM.
#' @param predictions_column_name The name of the column in the data frame `y` that contains the ML caller's predictions.
#'
#' @return An object representing the trained HMM after the Baum-Welch training.
#'
train_hmm <- function(y, base_hmm, hmm_bw_max_iter, bw_delta, train_chrom, predictions_column_name) {
  
  #function that trains a first order HMM on the selected chromosome
  #using the baum-welch algorithm 
  #and ML call's predictions as observations of the hidden states
  #taking only the train chromosome
  train <- subset(y, y$Chr == train_chrom)
  #order by chromosome and position
  train <- train[order(train$Chr, train$Start),]
  #the observer of the state is the ML caller, get its observations (ML call's prediction)
  train_observations <- as.character(train[[predictions_column_name]])
  #use baum-welch to estimate transition and emission probabilities given the observations of the hidden state
  #trace(baumWelch, edit="nano")
  bw = baumWelch(base_hmm,train_observations, maxIterations =  hmm_bw_max_iter, delta = bw_delta, pseudoCount = 1e-6)
  # print(bw$hmm)
  return(bw)
}