########### Function that outputs the correctly formatted ML model identifier

get_model_id <- function(arg){
  arg <- tolower(arg)
  if (arg %in% c("mlp", "nn", "neuralnetwork", "neural_network")){
    model_id <- "NN"
  }
  else if (arg %in% c("rf", "randomforest", "random_forest")){
    model_id <- "RF"
  }
  else if (arg %in% c("svc", "svm", "supportvectorclassifier", "support_vector_classifier")){
    model_id <- "SVC"
  }
  else{
    stop("The supplied model identifier is wrong. Please check.")
  }
  return(model_id)
}