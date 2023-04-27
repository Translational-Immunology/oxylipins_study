library(tidyverse)
library(plotly)


#' Produce a basic statistics table
#' @param data dataset (with one column stimulus id) column 'sex', 'age', 'donor_id' are ignored
#' @param list_function list of function to apply c('variance', 'standard_deviation', 'median', 'mean')
get_basic_statistic_table <- function(data, list_function = c('variance', 'standard_deviation', 'median', 'mean', 'cv', 'q1', 'q3'), data_transform = NULL){

  cv <- function(x){if(mean(x) == 0){
    return(NaN)} else {
    return(100*(sd(x)/mean(x)))}
  }
  q1 <- function(x){
    return(as.numeric(as.list(summary(x)[['1st Qu.']])))}
  q3 <- function(x){return(as.numeric(as.list(summary(x)[['3rd Qu.']])))}

  allowed_func <- list('variance' = var,
                       'standard_deviation' = sd,
                       'median' =  median,
                       'mean' = mean,
                       'cv' = cv,
                       'q1' = q1,
                       'q3' = q3)

  data_table_ret <- NULL
  for(func in rev(list_function)){

    print('start : ')
    print(func)

    for(stimulus in unique(data$stimulus_id)){
      data_table_ret <- data %>%
        dplyr::mutate_at(dplyr::vars(- dplyr::any_of(c('donor_id','stimulus_id','sex','age'))), .funs = as.numeric) %>%
        {if(! is.null(data_transform)) dplyr::mutate_at(., dplyr::vars(- dplyr::any_of(c('donor_id','stimulus_id','sex','age'))), .funs = data_transform) else .} %>%
        dplyr::filter(stimulus_id == stimulus) %>%
        dplyr::select(- dplyr::any_of(c('donor_id','stimulus_id','sex','age'))) %>%
        dplyr::summarise_all(allowed_func[[func]]) %>%
        dplyr::mutate(stimulus = stimulus) %>%
        dplyr::mutate(metric = func) %>%
        dplyr::bind_rows(data_table_ret)
    }

    print('done -------')
  }

  return(data_table_ret)
}





