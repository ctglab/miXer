########### Function that extracts windows where the HMM has constant (altered) state
#' Create Windows for Analysis
#'
#' This function takes a result dataframe, along with other parameters, and creates windows for analysis.
#'
#' @param result_df The dataframe containing the results.
#' @param one_X_copy A logical value indicating whether the specimen has one X chromosome.
#' @param predictions_column_name The name of the column containing the predictions.
#' @param par_windows A dataframe with PAR regions.
#' @param break_ties The method to break ties for window classification. Default is "majority".
#' @param sidestep_left A logical value indicating whether to sidestep regions on the left side of windows. Default is FALSE.
#' @param sidestep_right A logical value indicating whether to sidestep regions on the right side of windows. Default is FALSE.
#' @param median_thr_delToDDel The median threshold for transitioning from DEL to DDEL. Default is -4.
#' @param median_thr_wtToDup The median threshold for transitioning from WT to DUP. Default is 0.5.
#' @param median_thr_DupToMdup The median threshold for transitioning from DUP to MDUP. Default is 0.81.
#'
#' @return A dataframe containing the windows for analysis.
#'
#' @examples
#' window_maker(result_df, one_X_copy, predictions_column_name, par_windows,
#'              break_ties = "majority", sidestep_left = FALSE, sidestep_right = FALSE,
#'              median_thr_delToDDel = -4, median_thr_wtToDup = 0.5, median_thr_DupToMdup = 0.81)
#'
#'

window_maker <- function(result_df, one_X_copy, predictions_column_name, par_windows,
                         show_wt_windows = F,
                         break_ties = "majority", sidestep_left = F, sidestep_right = F,
                         median_thr_delToDDel = -4, median_thr_wtToDup= 0.5, median_thr_DupToMdup = 0.81,
                         sliding_window_length = 3, score_threshold = 0, break_on_na =FALSE, min_contiguous_nas = 5){

  cols <-c("Chr", "Start", "End",
           "State", "Call", "CN",
           "ProbCall_noPenalties", "ProbCall", "p_error", "Median_NRC",  "Number_of_TR ",
           "window_length", "Mean_TR_length",
           "window_mean_mappability", "alt_post_prob_mean", "Flagged_Positions"
           )
  
  windows_df <- data.frame(matrix(ncol = length(cols), nrow = 0))
  
  colnames(windows_df) <- cols
  i <- sapply(windows_df, is.logical)
  windows_df[i] <- lapply(windows_df[i], as.character)
  chr_list <- unique(result_df$Chr)
  
  
  for (item in chr_list){
    #keep only 
    curr_chrom <- subset(result_df, result_df$Chr == item)
    row.names(curr_chrom) <- NULL
    #ordering dataframe by chromosome and start position
    curr_chrom <- curr_chrom[order(curr_chrom$Chr, curr_chrom$Start),]
    
    # Creating a copy of the chromosome to keep NA TRs for counting
    curr_chrom_copy <- curr_chrom
    
    if (!break_on_na){
      # Filter out rows where 'hmm_state' is NA
      curr_chrom <- curr_chrom[!is.na(curr_chrom$hmm_state), ]
    }
    
    if ( (one_X_copy == TRUE)&& (item %in% c("chrx", "chrX", "ChrX", "Chrx", "X", "x"))){
      alt_states <- c("DDEL","WT", "DUP")
      prob_names <- c("DDEL_post_prob", "DEL_post_prob", "WT_post_prob", "DUP_post_prob")
    } else{
      alt_states <- c("DEL", "DUP")
      prob_names <- c("DEL_post_prob", "WT_post_prob", "DUP_post_prob")
    }
    
    hmm_state_probabilities <- apply(curr_chrom[, prob_names], 1, max)
    curr_chrom$mps_prob <- hmm_state_probabilities
    
    # prefiltering problematic exons (i.e. those in regions where it's low mapp/low prob)
    scores <- numeric(nrow(curr_chrom))
    for (j in 1:(nrow(curr_chrom) - sliding_window_length + 1)) {
      window_probs <- curr_chrom[j:(j + sliding_window_length - 1), "mps_prob"]
      window_mapps <- curr_chrom[j:(j + sliding_window_length - 1), "Mappability"]
      window_score <- mean(window_probs * window_mapps)
      scores[j:(j + sliding_window_length - 1)] <- pmax(scores[j:(j + sliding_window_length - 1)], window_score, na.rm = TRUE)
    }
    
    # Assign the last score to remaining rows
    scores[(nrow(curr_chrom) - sliding_window_length + 2):nrow(curr_chrom)] <- scores[nrow(curr_chrom) - sliding_window_length + 1]
    curr_chrom$flagged <- scores > score_threshold #keeping only good quality ones
    
    #curr_chrom <- curr_chrom[curr_chrom$flagged == TRUE, ]
    # force set the states for subpar positions to WT
    curr_chrom$hmm_state[curr_chrom$flagged == FALSE] <- "WT"
    
    
    #creating a new column for altered state
    #this will be used to create ALT windows which will be "genotyped" later
    #as in: I can have windows with the following hmm state sequence ("DEL" "DEL" "DUP" "DEL")
    #which could be divided in three different windows or could be considered only one window
    #if they are to be considered one window, the cnv type will be decided by majority voting
    curr_chrom$altered_state <- curr_chrom$hmm_state
    curr_chrom$altered_state[curr_chrom$altered_state %in% alt_states] <- "ALT"

    #calculate the length of the sequences of consecutive values in hmm_state column using rle
    rle_lengths <- rle(curr_chrom$altered_state)$lengths
    
    #create a new vector that indicates the group to which each row of the df belongs using the rep function
    group <- rep(1:length(rle_lengths), rle_lengths)
    #add the vector to the dataframe as a column
    curr_chrom$group <- group
    #split the data frame in subgroups using group column <--- these are the windows! (WT windows included)
    subgroups <- split(curr_chrom, curr_chrom$group)

    #log current chromosome
    win_chrom = item
    
    
    if (break_on_na){
      # Breaking windows if too many TRs are flagged as NA (e.g. under mappability threshold)
      # Logic to split the subgroups again if more than 5 contiguous NAs are found
      i <- 1
      while (i <= length(subgroups)) {
        sg <- subgroups[[i]]
        
        # Identify indices of hmm_state that are NA
        na_indices <- which(is.na(sg$hmm_state))
        
        if (length(na_indices) > 0) {
          # Find contiguous sequences of NAs
          contiguous_na <- rle(c(diff(na_indices) == 1, FALSE))$lengths
          start_indices <- na_indices[cumsum(contiguous_na) - contiguous_na + 1]
          end_indices <- na_indices[cumsum(contiguous_na)]
          
          # Check if there is any sequence of NAs > 10
          long_na_seq <- which(contiguous_na > min_contiguous_nas)
          
          if (length(long_na_seq) > 0) {
            # For simplicity, split at the first such sequence found
            first_na_seq <- long_na_seq[1]
            split_start <- start_indices[first_na_seq]
            split_end <- end_indices[first_na_seq]
            
            # Split sg into two parts
            before_na <- sg[1:(split_start - 1), ]
            after_na <- sg[(split_end + 1):nrow(sg), ]
            
            # Update subgroups
            subgroups <- append(subgroups[1:(i-1)], list(before_na, after_na), after = i-1)
            subgroups <- subgroups[-(i+2)] # Remove the original sg
          } else {
            # Move to the next element if no long contiguous NA sequence is found
            i <- i + 1
          }
        } else {
          # Move to the next element if there are no NAs
          i <- i + 1
        }
      }
      
      # Logic to remove remaining NA values in the subgroups
      for (i in seq_along(subgroups)) {
        sg <- subgroups[[i]]
        
        # Remove rows with NAs in the hmm_state column
        sg <- sg[!is.na(sg$hmm_state), ]
        
        # Update the subgroup in the list
        subgroups[[i]] <- sg
      }
      # Remove empty subgroups (with zero rows)
      subgroups <- subgroups[sapply(subgroups, nrow) > 0]
    }
    
    #for each subgroup
    for (i in 1:length(subgroups)) {
      #print(paste("Window", i, "out of", length(subgroups)))
      sg <- subgroups[[i]]
      
      # Count the number of CONSIDERED TRs (i.e. those with sufficient mappability, will be used to calculate mean of probabilities, mean of mappabilities, etc etc)
      win_in_regions <- nrow(sg)
      # I can't consider the tot_win_in_regions since these can contain  a certain number of NA regions with insufficient mappabilities
      # Which do not contribute to the total metrics. This results in a wrong metrics calculation. 
      # As in. if i have a window with 10 TRs and 3 are low-mappability, if i consider tot_win_in_regions i get the mean on 10 TRs of which only 7 contribute to the accumulation of the metrics
      # Resulting in an undervaluing (and thus involuntary filtering)
      # However, for final tool metrics calculation i musst consider tot_win_in_regions
      # As this is the number of TRs considered in my window and i must consider ALL the TRs belonging to the sequencing target.
      # I think this is clear now, i'm going to rest after testing
      
      
      #get window start and end position
      win_start = as.numeric(sg[1, "Start"])
      win_end = as.numeric(sg[nrow(sg), "End"])
      #calculate window length
      win_length <- as.numeric(win_end - win_start)
      #calculate real number of TRs contained in window
      #also considering, if present, NA TRs (i.e. flagged for low mappaility) (see above for rational)
      # Subset the original copy of the dataframe to include rows with Start >= win_start and End <= win_end
      included_rows <- curr_chrom_copy[curr_chrom_copy$Start >= win_start & curr_chrom_copy$End <= win_end, ]
      tot_win_in_regions <- nrow(included_rows)
      #tot_win_in_regions <- nrow(sg)

      #calculate total exon length
      win_total_exon_length <- sum(sg$Length)
      #### ML model observations list
      ml_window_obs <- as.numeric(unlist(sg[[predictions_column_name]]))
      #### HMM state list 
      hmm_window_states <- c(as.character(unlist(sg[["hmm_state"]]) ) )
      #### HMM state probability sequence, handling NAs
      hmm_state_probabilities <- apply(sg[, prob_names], 1, max)
      #sg$mps_prob <- hmm_state_probabilities
      #### Call error probability sequence
      p_error = 1 - hmm_state_probabilities #TODO add column
      flagged_positions <- sum(sg$flagged)
      
      #accumulate window state probability
      win_total_posterior_probability <- sum(hmm_state_probabilities)
      #accumulate window error probability
      win_total_p_error <- sum(p_error)
      #calculate window state probability (mean of posterior probabilities)
      alt_post_prob_mean <- as.numeric(win_total_posterior_probability/win_in_regions)
      #calculate window error probability
      alt_post_p_error_mean <- as.numeric(win_total_p_error/win_in_regions)
      #accumulate the mappability of the window
      win_total_mappability <- sum(sg$Mappability)
      #calculate mean mappability value
      #win_mean_mapp <- as.numeric(win_total_mappability/win_in_regions)
      
      #if (curr_chrom$Chr == "chr16"){
      #  print(win_start)
      #  print(win_end)
      #  print(as.numeric(unlist(sg["Mappability"])))}
      
      win_mean_mapp <- mean(as.numeric(unlist(sg$Mappability)))
      
      #calculate mean exon length
      length_mean <- as.numeric(win_total_exon_length/win_in_regions)
      #### All NRC poolNorm values
      all_win_nrc <- na.omit(as.numeric(unlist(sg[["NRC_poolNorm"]])))
      #### Calculate Total distance between regions in window
      tot_dist <- sum(sg[-1, "Start"] - sg[-nrow(sg), "End"])
      #calculate mean distance between regions
      density <- win_total_exon_length/win_length
      #calculate the basic window confidence
      basic_win_confidence <- alt_post_prob_mean * win_mean_mapp
      # Floor the result at the third decimal point
      basic_win_confidence<- floor(basic_win_confidence * 1000) / 1000
      
      # #penalizing windows which have exons that are far away from each other
      # 
      # #TODO removed, option to keep and better tune it?
      # if (tot_dist > 0.999*win_length){
      #   win_confidence <- basic_win_confidence - exp(-300*(density))
      # }else if (tot_dist > 0.5*win_length){
      #   win_confidence <- basic_win_confidence - 0.8*exp(-350*(density))
      # }else{
      #   win_confidence <- basic_win_confidence
      # }

      win_confidence <- basic_win_confidence
      #calculate median nrc
      median_nrc <- median(all_win_nrc)
      #get number of unique observations for the SVM
      num_unique_ml_obs <- length(unique(na.omit(ml_window_obs)))
      
      #get number of unique states of the hmm
      unique_hmm_states <- length(unique(hmm_window_states))
      
      #definining a dataframe from the three vectors, for later usage
      tmp <- data.frame(ml_window_obs, hmm_window_states, hmm_state_probabilities)
      tmp$hmm_state_probabilities <- as.numeric(tmp$hmm_state_probabilities)
      
      #get mean posterior probability by hmm state
      ddels_mean_prob <- mean(subset(tmp, tmp$hmm_window_states == "DDEL")$hmm_state_probabilities)
      dels_mean_prob <- mean(subset(tmp, tmp$hmm_window_states == "DEL")$hmm_state_probabilities)
      dups_mean_prob <- mean(subset(tmp, tmp$hmm_window_states == "DUP")$hmm_state_probabilities)
      
      #fix missing probabilities by setting them to zero
      if (is.nan(ddels_mean_prob)){
        ddels_mean_prob = 0
      }
      if (is.nan(dels_mean_prob)){
        dels_mean_prob = 0
      }
      if (is.nan(dups_mean_prob)){
        dups_mean_prob = 0
      }
      
      #handling mixed state windows by means of majority voting
      #handling different HMM configurations for autosomes and X chromosome
      if ( (one_X_copy == TRUE)&& (item %in% c("chrx", "chrX", "ChrX", "Chrx", "X", "x"))){
        probs <- data.frame(alt_states, c(ddels_mean_prob, dels_mean_prob, dups_mean_prob))
      } else {
        probs <- data.frame(alt_states, c(dels_mean_prob, dups_mean_prob))
      }

      colnames(probs) <- c("State", "Mean_probs")
      if (unique_hmm_states > 1){
        #this means that the HMM shifted from DEL to DUP or viceversa
        #the window cnv type will be voted by majority
        #TODO other ways to treat it?
        #get number of occurrences by state
        
        state_frequencies <- as.data.frame(table(hmm_window_states))
        #order by number of appearances (decreasing)
        state_frequencies <- state_frequencies[order(-state_frequencies$Freq), ]
        #check if the first two rows have the same number of elements
        
        if (state_frequencies[1,]$Freq == state_frequencies[2,]$Freq){
          #print("Woah, found a window with equal number of DEL/DUP states, much rare")
          #check which state has higher mean posterior probability
          probs <- probs[order(-probs$Mean_probs),]
          if (probs[probs$State == as.character(state_frequencies[1,]$hmm_window_states), ]$Mean_probs >= probs[probs$State == as.character(state_frequencies[2,]$hmm_window_states), ]$Mean_probs){
            win_type <- as.character(state_frequencies[1,]$hmm_window_states)
          }
          else{
            win_type <- as.character(state_frequencies[2,]$hmm_window_states) #TODO should be state_frequencies[2,]$hmm_window_states
          }
        }
        else{
          #otherwise, assign to win_type the value of the most occurrent state
          win_type = as.character(state_frequencies[1, ]$hmm_window_states)
          #TODO: is it possible/useful to lower the confidence if this happens?
        }
      } else{ win_type = as.character(unique(sg$hmm_state) )} #No shifts between ALT states, business as usual
      
      
      #check if sex is specified as male
      #or, in the absence of a specification,
      #if one copy of the X chromosome has been detected by the guess_23 function with the SVM
      
      #ignoring PAR regions
      #I'm not in a pseudoautosomal region IF:
      #my window ends before the start of the first region
      #OR
      #my window starts after the first region and stops before the start of the second region
      #OR
      #my window starts after the end of the second region
      nonPar_cond <- (win_end <= par_windows[1,]$Start) || (win_start > par_windows[1,]$End && win_start < par_windows[2,]$Start) || (win_start > par_windows[2,]$End)
      #non-Pseudo autosomal region condition, as in nonPar_con == TRUE --> not in a pseudo-Autosomal region
      
      #prepare for "genotipization" 
      ml_frequencies <- as.data.frame(table(ml_window_obs))
      #order by number of appearances (decreasing)
      ml_frequencies <- ml_frequencies[order(-ml_frequencies$Freq), ]
      
      #If I'm working with a one X chromosome specimen (or perhaps explicitly labelled as M)
      #And I'm in chrX
      single_chrx_cond <- win_chrom %in% c("chrX", "chrx", "ChrX", "Chrx", "X", "x") && one_X_copy == TRUE
      
      if (single_chrx_cond == TRUE){
        #Check if I'm in a Pseudoautosomal region
        if (nonPar_cond == TRUE){ #Not in a pseudo-Autosomal region
          #If I'm in DEL state it actually means WT
          if (win_type == "DDEL"){#If 'm in DDEL it means one copy deletion
            win_type = "DEL"
          } else if (win_type == "DEL"){ #If I'm in DEL it means WT
            win_type = "WT"
          }else if (win_type == "WT"){ #If I'm in WT state it actually means one copy duplication 
            win_type = "DUP"
          } else if(win_type == "DUP"){ #If I'm in DUP state it actually means multiple copy duplication
            win_type = "MDUP"
          }
        }
        else{
          
          #I'm in a Pseudoautosomal region and considering a male! The shock
          #Must characterize the DUP windows
          #If the most frequent prediction is -1 ignore it for characterization
          if(ml_frequencies[1,]$ml_window_obs == -1){
            ml_frequencies = ml_frequencies[-1,]
          }
          
          if (nrow(ml_frequencies) > 0) { # As it turns out, there can be events where only -1 observations are present
            #In this case, no changes can be made to the window
            
          if (win_type == "DUP"){#if, at the end of the execution, the call is DUP, there's a possibility that it still needs to be changed to MDUP
            #check if the first two rows have the same number of elements
            if (nrow(ml_frequencies) > 1){
              if (ml_frequencies[1,]$Freq == ml_frequencies[2,]$Freq){
                #characterize the window using the NRC_poolNorm median threshold
                if (median_nrc > median_thr_DupToMdup){
                  win_type = "MDUP"
                }
              }
            } else if (ml_frequencies[1,]$ml_window_obs == 2){ #I have to change DUP to MDUP
              win_type = "MDUP"
            }
          } 
            
          }
          
        }
      } else {
        #I'm either in the autosomes or dealing with a female sample
        #normal "genotyping" ensues
        #If the most frequent prediction is 0 ignore it for characterization
        if(ml_frequencies[1,]$ml_window_obs == 0){
          ml_frequencies = ml_frequencies[-1,]
        }
        
        if (nrow(ml_frequencies) > 0) { # As it turns out, there can be events where only 0 observations are present
          #In this case, no changes can be made to the window
          
          if (win_type == "DUP"){#if, at the end of the execution, the call is DUP, there's a possibility that it still needs to be changed to MDUP
            #check if the first two rows have the same number of elements
            
            if (nrow(ml_frequencies) > 1){
              if (ml_frequencies[1,]$Freq == ml_frequencies[2,]$Freq){
                #characterize the window using the NRC_poolNorm median threshold
                if (median_nrc > median_thr_DupToMdup){
                  win_type = "MDUP"
                }
              }
              else if (ml_frequencies[1,]$ml_window_obs == 2){ #I have to change DUP to MDUP
                win_type = "MDUP"
              }
            }
            else if (ml_frequencies[1,]$ml_window_obs == 2){
              win_type = "MDUP"
            }
          }
          
          if (win_type == "DEL"){#if, at the end of the execution, the call is DEL, there's a possibility that it still needs to be changed to DDEL
            #check if the first two rows have the same number of elements
            if (nrow(ml_frequencies) > 1){
              if (ml_frequencies[1,]$Freq == ml_frequencies[2,]$Freq){
                #characterize the window using the NRC_poolNorm median threshold
                if (median_nrc < median_thr_delToDDel){
                  win_type = "DDEL"
                }
              }
              else if (ml_frequencies[1,]$ml_window_obs == -2){ # I have to change DEL to DDEL
                win_type = "DDEL"
              }
            }
            else if (ml_frequencies[1,]$ml_window_obs == -2){ #I have to change DEL to DDEL
              win_type = "DDEL"
            }
          }
        }
      }
      #check if the user wants to enlarge the window to the left
      #this should be done only if the SVM prediction is altered
      if (sidestep_left == TRUE && single_chrx_cond == FALSE){ #skipping chrX for this (no idea how to check for PAR regions)
        
        #get the SVM prediction of the exon before the one considered for the start
        #check that this is not the first subsequence (no exon before the starting one)
        if (i > 1){
          #if it is not, check that the window before this is WT (otherwise I may end up with overlapping altered windows)
          #this can be done using the previous info data_frame (which contains the state of the previous subsequence)
          prev_state = info["hmm_state"]
          if (prev_state == "WT"){
            #get the previous subsequence
            prev_subseq = subgroups[[i-1]]
            #get the last exon of the subsequence and check the SVM output
            svm_call = as.numeric(prev_subseq[nrow(prev_subseq), ][[predictions_column_name]])
            #get the starting coordinates of the last exon of the subsequence
            new_start = as.numeric(prev_subseq[nrow(prev_subseq), "Start"])
            #get the probability of the SVM call
            svm_call_prob = as.numeric(prev_subseq[nrow(prev_subseq), ][[paste0(predictions_column_name, "_proba")]])
            distance = new_start - win_start
            
            if (svm_call != 0 && svm_call_prob > 0.7 && distance < 10^5){ #TODO as of now altered goes, check for concordance with the HMM state?
              #change the window start to that of the last exon of the previous subsequence
              if (svm_call >= 1 && win_type %in% c("MDUP", "DUP")){
                win_start = new_start
                win_in_regions = win_in_regions +1
                tot_win_in_regions = tot_win_in_regions +1
              }
              if ( svm_call <= -1 && win_type %in% c("DDEL", "DEL")){
                win_start = new_start
                win_in_regions = win_in_regions +1
                tot_win_in_regions = tot_win_in_regions +1
              }
                  
            }
          }
          
        }
      }
      if (sidestep_right == TRUE && single_chrx_cond == FALSE){ #skipping chrX for this (no idea how to check for PAR regions)
        
        #get the SVM prediction of the exon after the one considered for the end
        #check that this is not the last subsequence (no exon after the ending one)
        if (i < length(subgroups)){
          #if it is not, check that the window after this is WT (otherwise I may end up with overlapping altered windows)
          #get the next subsequence
          next_subseq = subgroups[[i+1]]
          #get the state of the HMM in the first exon of the next window and pray it is WT
          next_state = next_subseq[1,"hmm_state"]
          if (next_state == "WT"){
            #get the first exon of the subsequence and check the SVM output
            svm_call = as.numeric(next_subseq[1, ][[predictions_column_name]])
            #get the ending coordinates of the first exon of the subsequence
            new_end = as.numeric(next_subseq[nrow(next_subseq), "End"])
            #get the probability of the SVM call
            svm_call_prob = as.numeric(next_subseq[nrow(next_subseq), ][[paste0(predictions_column_name, "_proba")]])
            distance = new_end - win_end
            
            if (svm_call != 0 && svm_call_prob > 0.7 && distance < 10^5){ #TODO as of now altered goes, check for concordance with the HMM state?
              #change the window start to that of the last exon of the previous subsequence
              if (svm_call >= 1 && win_type %in% c("MDUP", "DUP")){
                win_end = new_end
                win_in_regions = win_in_regions +1
                tot_win_in_regions = tot_win_in_regions +1
              }
              if ( svm_call <= -1 && win_type %in% c("DDEL", "DEL")){
                win_end = new_end
                win_in_regions = win_in_regions +1
                tot_win_in_regions = tot_win_in_regions +1
              }
              
            }
          }
          
        }
      }

      #creating State and Call df columns
      if (win_type %in% c("DDEL", "DEL")){
        State <- "DEL"
        if (win_type == "DDEL"){
          hmm_copynumber_estim <- "-2"
          cn_est <- -2
        }else{
            hmm_copynumber_estim <- "-1"
            cn_est <- -1
            }
        
      } else if (win_type %in% c("DUP", "MDUP")){
        State <- "DUP"
        if (win_type == "DUP"){
          hmm_copynumber_estim <- "+1"
          cn_est <- 1
        }else{
            hmm_copynumber_estim <- "2+"
            cn_est <- 2}
      } else{
        State <- "WT"
        hmm_copynumber_estim <- "0"
        cn_est <- 0
      }
      #creating Total copy number df column
      if (single_chrx_cond == FALSE || nonPar_cond == FALSE){ #Either 2 copies of the X-Chr or in a pseudo-Autosomal region
        total_copies <- 2 + cn_est
      } else{ #I'm in a male chrX non-pseudo autosomal region
          total_copies <- 1 + cn_est
      }
      #turning total copy number from int to str to change "4" to "4+"
      #Since we can't distinguish if there are 4 copies or more of a TR
      if (total_copies == 4) {
        total_copies <- "4+"
      } else {
        total_copies <- as.character(total_copies)
      }
      
      #end of window calculations
      #making a vector with all the window info

      info <- c(win_chrom, win_start, win_end,
                State, hmm_copynumber_estim, total_copies,
                basic_win_confidence, win_confidence, alt_post_p_error_mean, median_nrc,tot_win_in_regions,
                win_length, length_mean, 
                win_mean_mapp, alt_post_prob_mean, cn_est, flagged_positions)
      
      info <- data.frame(matrix(info, 1))
      
      colnames(info) <- colnames(windows_df)
      #adding the aforementioned vector to the windows dataframe
      windows_df <- rbind(windows_df, info)
      #moving on to next window
    }
  }
  
  #converting factors 
  windows_df$Chr <- as.character(windows_df$Chr)
  windows_df$Start <- as.numeric(as.character(windows_df$Start))
  windows_df$End <- as.numeric(as.character(windows_df$End))
  windows_df$window_length <- as.numeric(as.character(windows_df$window_length))
  windows_df$Mean_TR_length <- as.numeric(as.character(windows_df$Mean_TR_length))
  windows_df$Call <- as.character(windows_df$Call)
  windows_df$CN <- as.character(windows_df$CN)
  windows_df$Number_of_TR  <- as.numeric(as.character(windows_df$Number_of_TR ))
  windows_df$alt_post_prob_mean <- as.numeric(as.character(windows_df$alt_post_prob_mean))
  windows_df$ProbCall_noPenalties <- as.numeric(as.character(windows_df$ProbCall_noPenalties))
  windows_df$ProbCall  <- as.numeric(as.character(windows_df$ProbCall))
  windows_df$p_error <- as.numeric(as.character((windows_df$p_error)))
  #windows_df$Mean_NRC <- as.numeric(as.character(windows_df$Mean_NRC))
  windows_df$Median_NRC <- as.numeric(as.character(windows_df$Median_NRC))
  #windows_df$num_unique_ml_obs <- as.numeric(windows_df$num_unique_ml_obs)
  #windows_df$unique_hmm_states <- as.numeric((windows_df$unique_hmm_states))
  windows_df$Flagged_Positions <- as.numeric(as.character(windows_df$Flagged_Positions))
  #Removing WT windows if not explicitly requested
  if (show_wt_windows == F) {
    windows_df = subset(windows_df, windows_df$State != "WT")
  }
  return(windows_df)
}