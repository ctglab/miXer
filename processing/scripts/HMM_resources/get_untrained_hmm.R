########### Function that returns the correct HMM architecture

#' Initialize and return an untrained Hidden Markov Model (HMM) object for a specific set of unique observations.
#'
#' @param unique_observations A vector of unique observations.
#' @param is_male A boolean indicating whether the observations are from male individuals.
#'
#' @return An untrained HMM object with the specified states, symbols, starting probabilities,
#'   transition probabilities, and emission probabilities.
#'
get_untrained_hmm <- function(unique_observations){
  
  #get hmm states
  states <- get_hmm_states(unique_observations)
  #get the symbols for the states emissions
  symbols <- get_hmm_emission_symbols()
  #get the starting probabilities of the hmm states
  start_probs <- get_hmm_starting_probabilities()
  #get the base transition probability matrix
  trans_probs <- get_hmm_transition_probabilities(symbols)
  #get the base emission probability matrix
  emiss_probs <- get_hmm_emission_probabilities(symbols)
  #initialize the hmm object
  base_hmm <- initHMM(states, symbols, startProbs = start_probs, transProbs = trans_probs, emissionProbs = emiss_probs)
  return(base_hmm)
}

#' Get the Hidden Markov Model (HMM) states based on unique observations.
#'
#' @param unique_observations A vector of unique observations.
#'
#' @return A character vector representing the HMM states.
#'
get_hmm_states <- function(unique_observations){
  
  states <- c("DEL", "WT", "DUP")
  return(states)
}

#' Determine and return the emission symbols for a set of unique observations in a Hidden Markov Model (HMM).
#'
#' The emission symbols represent the possible values emitted by the states in the HMM.
#'
#' @param unique_observations A vector of unique observations.
#'
#' @return A vector of emission symbols representing the possible values emitted by the states.
#'
get_hmm_emission_symbols <- function(){
  symbols = c("-2", "-1","0", "1", "2")
  return(symbols)
}

#' Get the starting probabilities for the Hidden Markov Model (HMM) states.
#'
#' @param is_male A boolean indicating whether the observations are from male individuals.
#'
#' @return A numeric vector representing the starting probabilities for the HMM states.
#'
get_hmm_starting_probabilities <- function(){
  start_probs <- c(0.025, 0.95, 0.025)
  return(start_probs)
}

get_hmm_emission_probabilities <- function(symbols, r = 0.9){
  
  #SVM not capable of detecting double deletions and multiple duplications
  if (length(symbols) == 3){
    l1 <- (1-r)/2
    #emiss_prob <- c(em_prob(-1), em_prob(0), em_prob(1))
    sm1_empr <- c(r, l1, l1)
    s0_empr <- c(l1, r,l1)
    s1_empr <- c(l1,l1, r)
    
  }
  else{
      if ("-2" %in% symbols){
        l1 <- (1-r)/3
        sm1_empr <- c(r/2, r/2, l1*2, l1)
        s0_empr <- c(l1, l1, r, l1)
        s1_empr <- c(l1, l1, l1, r)
        
        if ("2" %in% symbols){

          l1 <- (1-r)/4
          sm1_empr <- c(r/2, r/2, l1*2, l1, l1)
          s0_empr <- c(l1, l1, r, l1, l1)
          s1_empr <- c(l1, l1, l1*2, r/2, r/2)
        }
      
    }
  }
  emission_probs <-  matrix(c(sm1_empr,s0_empr,s1_empr), byrow = TRUE, nrow =3)
  
  return(emission_probs)
}

#' Calculate and return the emission probability matrix for a Hidden Markov Model (HMM).
#'
#' This function calculates the emission probabilities for each symbol in the given set of symbols,
#' considering the characteristics of the observations (e.g., gender, chromosome).
#'
#' @param symbols A vector of symbols representing the observed emissions.
#' @param is_male A boolean indicating whether the observations are from male individuals.
#' @param r The emission probability for the reference symbol (0).
#'
#' @return A matrix representing the emission probabilities.
#'
get_hmm_transition_probabilities <- function(symbols, s = 0.9){
  
  l2 <- (1-s)/2
  # sm1_trpr0 <- c(s, l2, l2)
  # s0_trpr0 <- c(l2,s,l2)
  # s1_trpr0 <- c(l2,l2,s)
  sm1_trpr0 <- c(s, l2, l2)
  s0_trpr0 <- c(l2,s,l2)
  s1_trpr0 <- c(l2,l2,s)

  transition_probs <- matrix(c(sm1_trpr0, s0_trpr0,s1_trpr0), byrow = TRUE, nrow =3)
  return(transition_probs)
}
