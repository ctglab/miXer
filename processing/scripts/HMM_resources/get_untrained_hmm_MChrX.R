########### Function that returns the correct HMM architecture
#Specific version for One ChrX samples

#unique observation: unique(curr_chrom[[predictions_column_name]])

#' Initialize and return an untrained Hidden Markov Model (HMM) object for a specific set of unique observations.
#'
#' @param unique_observations A vector of unique observations.
#'
#' @return An untrained HMM object with the specified states, symbols, starting probabilities,
#'   transition probabilities, and emission probabilities.
#'
get_untrained_hmm_MChrX <- function(unique_observations){
  #get hmm states
  states <- get_hmm_states_MChrX(unique_observations)
  #get the symbols for the states emissions
  symbols <- get_hmm_emission_symbols_MChrX()
  #get the starting probabilities of the hmm states
  start_probs <- get_hmm_starting_probabilities_MChrX()
  #get the base transition probability matrix
  trans_probs <- get_hmm_transition_probabilities_MChrX(symbols)
  #get the base emission probability matrix
  emiss_probs <- get_hmm_emission_probabilities_MChrX(symbols)
  #initialize the hmm object
  base_hmm <- initHMM(states, symbols, startProbs = start_probs, transProbs = trans_probs, emissionProbs = emiss_probs)
  return(base_hmm)
}

#' Get the Hidden Markov Model (HMM) states for the MChrX model.
#'
#' This function returns the states used in the Hidden Markov Model (HMM) for the MChrX model.
#' The states represent different types of observations: DDEL, DEL, WT, and DUP.
#'
#' @param unique_observations A vector of unique observations.
#'
#' @return A character vector containing the HMM states.
#'
get_hmm_states_MChrX <- function(unique_observations){
  
  states <- c("DDEL", "DEL", "WT", "DUP")
  
  return(states)
}

#' Get the emission symbols for the Hidden Markov Model (HMM) based on unique observations.
#'
#' This function determines the emission symbols for the HMM based on the provided unique observations.
#' The emission symbols are determined as follows:
#' - If "-2" is present in the unique observations, the symbols include "-2", "-1", "0", and "1".
#' - If both "-2" and "2" are present in the unique observations, the symbols include "-2", "-1", "0", "1", and "2".
#' - If "-2" is not present in the unique observations, the symbols include "-1", "0", and "1".
#'
#' @param unique_observations A vector of unique observations.
#'
#' @return A vector containing the emission symbols for the HMM.
#'
get_hmm_emission_symbols_MChrX <- function(){
      symbols = c("-2", "-1","0", "1", "2")
  return(symbols)
}

#' Retrieve the starting probabilities for the Hidden Markov Model (HMM) states.
#'
#' @return A numeric vector containing the starting probabilities for the HMM states.
#'
get_hmm_starting_probabilities_MChrX <- function(){

  start_probs <- c(0.0125,0.95, 0.0125, 0.025)
  
  return(start_probs)
}

#' Compute the emission probabilities for a specific set of symbols in the MChrX HMM.
#'
#' This function calculates the emission probabilities for a given set of symbols in the MChrX HMM.
#' The emission probabilities describe the likelihood of emitting a particular symbol from each state
#' of the HMM.
#'
#' @param symbols A vector of symbols for which the emission probabilities are computed.
#' @param r A parameter controlling the emission probabilities for the states.
#'           Default is 0.95.
#'
#' @return A matrix of emission probabilities, where rows represent states and columns represent symbols.
#'
get_hmm_emission_probabilities_MChrX <- function(symbols, r = 0.95){
  
    #working on a one-copy instance of the X chromosome (designated M or SVM-estimated M)
    #emiss_prob <- c(em_prob(-2),em_prob(-1), em_prob(0), em_prob(1))
    l1 <- (1-r)/3
    sm2_empr <- c(r, l1, l1, l1)
    sm1_empr <- c(l1, r, l1, l1)
    s0_empr <- c(l1, l1, r, l1)
    s1_empr <- c(l1, l1, l1, r)
    if ("2" %in% symbols){ #should be consequential to having -2 observations
      #emiss_prob <- c(em_prob(-2),em_prob(-1), em_prob(0), em_prob(1), em_prob(2))
      l1 <- (1-r)/4
      t = l1/8
      t2 = (l1 - t) + l1
      sm2_empr <- c(r, l1, l1, l1, l1)
      sm1_empr <- c(t, r+t2, l1, l1, l1)
      s0_empr <- c(l1, l1*2, r*3/6 +l1, r*2/6, r*1/6)
      s1_empr <- c(l1, l1, l1, l1, r)
      
    }
    
    emission_probs <-  matrix(c(sm2_empr, sm1_empr,s0_empr,s1_empr), byrow = TRUE, nrow =4)
  return(emission_probs)
  }

#' Calculate and return the transition probability matrix for a Hidden Markov Model (HMM) with specific symbols and transition parameter.
#'
#' @param symbols The symbols for the emissions of the HMM states.
#' @param s The transition parameter, representing the probability of transitioning to the same state in the next step (default: 0.9).
#'
#' @return The transition probability matrix for the HMM.
#'
get_hmm_transition_probabilities_MChrX <- function(symbols, s = 0.9){
  
  
    l2 <- (1-s)/3
    sm2_trpr0 <- c(s*2/5, s*3/5, l2*2/3, l2*1/3)
    sm1_trpr0 <- c(l2, s, l2, l2)
    s0_trpr0 <- c(l2, l2*2/3, s, l2*4/3)
    s1_trpr0 <- c(l2*1/3, l2*2/3, l2, s)
    
    
    transition_probs <- matrix(c(sm2_trpr0, sm1_trpr0, s0_trpr0,s1_trpr0), byrow = TRUE, nrow =4)
  
  
  return(transition_probs)
}
