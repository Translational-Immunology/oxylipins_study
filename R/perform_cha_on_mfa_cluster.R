
#' Wrapper around the factominer::MFA() function
#' @param list_dataset list of the diffrent datablock we want to consider for the study, each data block should have a donor_id column to merge them together
#' @param list_dataset_name names used to identify each dataset (block)
#' @param list_dataset_type the type of each dataset : different type are : l (logged), c (centered), s (scaled), n (qualitative)
perform_mfa <- function(list_dataset,
                        list_dataset_name = NULL,
                        list_dataset_type = rep('l', length(list_dataset)),
                        stimulus_filter = NULL){


  ########### Check + sanitize arguments

  # each dataset should have a name
  if(is.null(list_dataset_name)){
    list_dataset_name <- paste0('dataset_', 1:length(list_dataset))
  }
  assertthat::assert_that(length(list_dataset) == length(list_dataset_name))

  #each dataset should have a donor_id column
  for(i in 1:length(list_dataset)){
    assertthat::assert_that(class(list_dataset[[i]]) %in% c("tbl_df","tbl","data.frame"), msg = paste0('Dataset ',list_dataset_name[[i]], ' is not a regular dataset object'))
    assertthat::assert_that('donor_id' %in% colnames(list_dataset[[i]]), msg = paste0('Dataset ',list_dataset_name[[i]], 'does not have a donor_id column'))
  }

  #dataset should have unique column id
  all_names <- lapply(list_dataset, FUN = colnames) %>% unlist() %>% setdiff('donor_id') %>% length()
  all_unique_names <- lapply(list_dataset, FUN = colnames) %>% unlist() %>% setdiff('donor_id') %>% unique() %>% length()
  assertthat::assert_that(all_names == all_unique_names, msg = "Non unique column identifier in the datasets")


  # if dataset type is 'l' it has to be logp1 transformed
  for(i in 1:length(list_dataset_type)){
    if(list_dataset_type[[i]] == 'l'){
      list_dataset_type[[i]] <- 'c'
      dataset_tmp <- list_dataset[[i]] %>%
        dplyr::mutate_at(.vars = dplyr::vars(- dplyr::any_of(c('stimulus_id','donor_id'))), .funs =  function(x){return(log(x+1))})
      list_dataset[[i]] <- dataset_tmp
    }
  }

  remove_non_variating_column <- function(.data, other_variable_to_keep = c('donor_id','stimulus_id'), verbose = FALSE) {
    n_unique_value <- function(x) { return(length(unique(x[!is.na(x)]))) }
    to_keep <- dplyr::summarise_each(tbl = .data, funs = n_unique_value) %>%
      tidyr::gather() %>%
      dplyr::filter(value > 1) %>%
      dplyr::pull("key")
    if (!is.null(other_variable_to_keep)) {
      to_keep <- c(to_keep, other_variable_to_keep)
    }
    if (verbose) {
      print("The following column(s) will be removed :")
      print(setdiff(colnames(.data), c(to_keep, other_variable_to_keep)))
    }
    return(dplyr::select(.data, dplyr::any_of(to_keep)))
  }
  list_dataset <- list_dataset %>% lapply(remove_non_variating_column)

  ### merge dataset
  dataset_mfa <- NULL
  for(i in 1:length(list_dataset)){
    if(i == 1){
      dataset_mfa <- list_dataset[[i]]
    } else {
      key = 'donor_id'
      if('stimulus_id' %in% unlist(lapply(list_dataset, FUN = colnames))){
        key = c(key, 'stimulus_id')
      }

      if(!is.null(stimulus_filter)){
        list_dataset[[i]] <- list_dataset[[i]] %>%
          dplyr::filter(stimulus_id %in% stimulus_filter) %>%
          dplyr::mutate_at(dplyr::vars(- dplyr::any_of(c('stimulus_id', 'donor_id'))), as.numeric) %>%
          tidyr::drop_na() %>%
          remove_non_variating_column()
      }
      dataset_mfa <- dataset_mfa %>% dplyr::inner_join(list_dataset[[i]], by = key)
    }

  }


  dataset_mfa <- dataset_mfa %>%
  {if('stimulus_id' %in% colnames(dataset_mfa)) dplyr::mutate(., donor_id = paste0(donor_id, '_', stimulus_id)) else .} %>%
    tibble::column_to_rownames(var = 'donor_id') %>%
    dplyr::select(- dplyr::any_of('stimulus_id')) #%>% remove_non_variating_column()

  group <- list_dataset %>% lapply(remove_non_variating_column) %>% lapply(function(x){dplyr::select(x, -dplyr::any_of(c('stimulus_id','donor_id')))}) %>% lapply(function(x){return(length((colnames(x))))}) %>% unlist() %>% as.vector()


  ret_mfa <- FactoMineR::MFA(dataset_mfa, # data set avec toutes les variables (des différents blocs) et sans les variables stables (ie les molécules qui ne réagissent pas au stimulus)
                             group = group,
                             # Permet de désigner les blocs dans l'ordre d'apparition des variables (colonnes)
                             # dans le data set rentré en 1er argumt
                             # 1er nombre = nb de cyto analysées => pour designer le 1er bloc = bloc cytokines,
                             # 2e nombre = nb d'oxy analysées => pour designer le 2e bloc = bloc oxylipines
                             type = list_dataset_type,
                             # "c": variable quanti, centrées et non réduites ("c" pr centered);
                             # autres possibilités: "s": quanti, centrées et réduites ("s" pr scaled); ou
                             # ou quali (n);
                             #  ici je veux que toutes les variables quanti soient centrées mais non réduites => type = "c"
                             name.group = list_dataset_name,
                             graph = FALSE,
                             ncp=10)

  return(ret_mfa)
}


#' Wrapper around the factominer::MFA() function + plotly scatter plot
#' @param list_dataset list of the diffrent datablock we want to consider for the study, each data block should have a donor_id column to merge them together
#' @param list_dataset_name names used to identify each dataset (block)
#' @param list_dataset_type the type of each dataset : different type are : l (logged), c (centered), s (scaled), n (qualitative)
#' @param dataset_color_code a dataframe composed of two columns dataset_name and color_hex to specifiy the color of the markers
#' @param stimulus_filter filter on stimulus, dataset should contains a stimulus_id column
#' @param x component to display on axis x
#' @param y component to display on axis y
#' @param z component to display on axis z, if not NULL then the scatter plot shift to a 3D scatter plot
#' @param marker_opacity the plotly parameters marker_opacity
#' @param marker_size the plotly parameters marker_size, need to be increased in case of 3D scatter plot
#' @param text_color_code a dataframe composed of two columns dataset_name and color_hex to specifiy the color of the text
get_mfa_scatter_plot <- function(list_dataset,
                                 list_dataset_name = NULL,
                                 list_dataset_type = rep('l', length(list_dataset)),
                                 dataset_color_code = NULL,
                                 stimulus_filter = NULL,
                                 x = 1,
                                 y = 2,
                                 z = NULL,
                                 marker_opacity = 1,
                                 marker_size = 15,
                                 text_size = 12,
                                 text_color = NULL){

  # each dataset should have a name
  if(is.null(list_dataset_name)){
    list_dataset_name <- paste0('dataset_', 1:length(list_dataset))
  }
  assertthat::assert_that(length(list_dataset) == length(list_dataset_name))

  # and a type
  assertthat::assert_that(length(list_dataset_type) == length(list_dataset))
  assertthat::assert_that(length(setdiff(unique(list_dataset_type),c('c','s','l','n'))) == 0)


  # each dataset should have a color
  if(is.null(dataset_color_code)){
    dataset_color_code <- data_frame(dataset_name = list_dataset_name) %>%
      dplyr::mutate(color_hex = lapply(list_dataset_name,function(x){randomcoloR::randomColor()}))
  } else {
    #assertthat::assert_that(class(dataset_color_code) %in% c("tbl_df","tbl","data.frame"))
    assertthat::assert_that("dataset_name" %in% colnames(dataset_color_code), msg = "dataset_color_code should contains a column 'dataset_name'")
    assertthat::assert_that("color_hex" %in% colnames(dataset_color_code), msg = "dataset_color_code should contains a column 'color_code'")
    assertthat::assert_that(length(setdiff(list_dataset_name, dataset_color_code$dataset_name)) == 0, msg = "all dataset should have a color code")
  }

  mfa_object <- perform_mfa(list_dataset = list_dataset,
                            list_dataset_name = list_dataset_name,
                            list_dataset_type = list_dataset_type,
                            stimulus_filter = stimulus_filter)

  all_analyte_name <- factoextra::get_mfa_var(mfa_object, "quanti.var")$coord %>%
    as.data.frame() %>%
    rownames()
  group_name <- NULL
  for(i in 1:length(list_dataset)){
    group_name <- c(group_name, rep(list_dataset_name[[i]], length(intersect(colnames(list_dataset[[i]]), all_analyte_name))))
  }
  coordinate <- factoextra::get_mfa_var(mfa_object, "quanti.var")$coord %>%
    as.data.frame() %>%
    dplyr::mutate(group = group_name) %>%
    tibble::rownames_to_column('analyte')

  fig <- plotly::plot_ly(type = {if(! is.null(z)) "scatter3d" else "scatter"}, mode = 'markers+text')

  for(grp in unique(coordinate$group)){
    coordinate_tmp <- coordinate %>% dplyr::filter(group == grp)
    fig <- fig %>% plotly::add_trace(x = dplyr::pull(coordinate_tmp, paste0('Dim.',x)),
                                     y = dplyr::pull(coordinate_tmp, paste0('Dim.',y)),
                                     z = {if(! is.null(z)) dplyr::pull(coordinate_tmp, paste0('Dim.',z))},
                                     name = grp, #dplyr::pull(coordinate_tmp,'analyte'),
                                     text = dplyr::pull(coordinate_tmp,'analyte'),
                                     textposition="top center",
                                     textfont=list(
                                       family="sans serif",
                                       size=text_size,
                                       color = {if(is.null(text_color))dataset_color_code[dataset_color_code$dataset_name == grp,]$color_hex[[1]] else text_color}
                                     ),
                                     #labels = dplyr::pull(coordinate_tmp,'analyte'),
                                     type = {if(! is.null(z)) "scatter3d" else "scatter"},
                                     mode = 'markers+text',
                                     color = grp,
                                     size = rep(marker_size, time(nrow(coordinate_tmp))),
                                     sizes = c(marker_size,marker_size),
                                     opacity = marker_opacity,
                                     legendgroup = grp,
                                     marker = list(color = dataset_color_code[dataset_color_code$dataset_name == grp,]$color_hex[[1]],
                                                   opacity = marker_opacity,
                                                   line = list(color = dataset_color_code[dataset_color_code$dataset_name == grp,]$color_hex[[1]])))
  }

  return(fig)

}



#' Visualize CHA on MFA as a dendogram
#' @param list_dataset list of the diffrent datablock we want to consider for the study, each data block should have a donor_id column to merge them together
#' @param list_dataset_name names used to identify each dataset (block)
#' @param list_dataset_type the type of each dataset : different type are : l (logged), c (centered), s (scaled), n (qualitative)
#' @param dataset_color_code a dataframe composed of two columns dataset_name and color_hex to specifiy the color of the markers
#' @param stimulus_filter filter on stimulus, dataset should contains a stimulus_id column
#' @param title plot title
get_mfa_cha_dendogram <- function(list_dataset,
                                  list_dataset_name = NULL,
                                  list_dataset_type = rep('l', length(list_dataset)),
                                  dataset_color_code = NULL,
                                  stimulus_filter = NULL,
                                  title = "CHA post MFA on all analytes"){

  # each dataset should have a name
  if(is.null(list_dataset_name)){
    list_dataset_name <- paste0('dataset_', 1:length(list_dataset))
  }
  assertthat::assert_that(length(list_dataset) == length(list_dataset_name))
  assertthat::assert_that(length(list_dataset_type) == length(list_dataset))
  assertthat::assert_that(length(setdiff(unique(list_dataset_type),c('c','s','l','n'))) == 0)


  mfa_object <- perform_mfa(list_dataset = list_dataset,
                            list_dataset_name = list_dataset_name,
                            list_dataset_type = list_dataset_type,
                            stimulus_filter = stimulus_filter)

  #each dataset should have a color
  if(is.null(dataset_color_code)){
    dataset_color_code <- data_frame(dataset_name = list_dataset_name) %>%
          dplyr::mutate(color_hex = lapply(list_dataset_name,function(x){randomcoloR::randomColor()}))
  } else {
    assertthat::assert_that("dataset_name" %in% colnames(dataset_color_code), msg = "dataset_color_code should contains a column 'dataset_name'")
    assertthat::assert_that("color_hex" %in% colnames(dataset_color_code), msg = "dataset_color_code should contains a column 'color_code'")
    assertthat::assert_that(length(setdiff(list_dataset_name, dataset_color_code$dataset_name)) == 0, msg = "all dataset should have a color code")
  }


  all_analyte_name <- factoextra::get_mfa_var(mfa_object, "quanti.var")$coord %>%
          as.data.frame() %>%
          rownames()

  group_name <- NULL
  for(i in 1:length(list_dataset)){
    group_name <- c(group_name, rep(list_dataset_name[[i]], length(intersect(colnames(list_dataset[[i]]), all_analyte_name))))
  }

  coordinate <- factoextra::get_mfa_var(mfa_object, "quanti.var")$coord %>%
          as.data.frame()

  hc <- fastcluster::hclust(dist(coordinate) ,method = "ward.D2" )


  set_analyte_color <- function(plot){
    for(i in 1:length(dataset_name)){
          label_target <- intersect(colnames(list_dataset[[i]]), all_analyte_name)
          color_target <- dataset_color_code[dataset_color_code$dataset_name == list_dataset_name[[i]],]$color_hex[[1]]
          plot <- plot %>%
            dendextend::set("labels_col", color_target, labels = labels(plot)[which(labels(plot) %in% label_target)])
    }
    return(plot)
  }

  dend <- as.dendrogram(hc)
  return_plot <- dend %>%
          dendextend::set("labels_cex", 0.65 ) %>%
          set_analyte_color() %>%
          plot(main=title)

  return(return_plot)

}




#' Wrapper around the factominer::MFA() + factoextra::fviz()
#' @param list_dataset list of the diffrent datablock we want to consider for the study, each data block should have a donor_id column to merge them together
#' @param list_dataset_name names used to identify each dataset (block)
#' @param list_dataset_type the type of each dataset : different type are : l (logged), c (centered), s (scaled), n (qualitative)
#' @param stimulus_filter filter on stimulus, dataset should contains a stimulus_id column
get_mfa_scree_plot <- function(list_dataset,
                               list_dataset_name = NULL,
                               list_dataset_type = rep('l', length(list_dataset)),
                               stimulus_filter = NULL){

  # each dataset should have a name
  if(is.null(list_dataset_name)){
    list_dataset_name <- paste0('dataset_', 1:length(list_dataset))
  }
  assertthat::assert_that(length(list_dataset) == length(list_dataset_name))

  # and a type
  assertthat::assert_that(length(list_dataset_type) == length(list_dataset))
  assertthat::assert_that(length(setdiff(unique(list_dataset_type),c('c','s','l','n'))) == 0)

  mfa_object <- perform_mfa(list_dataset = list_dataset,
                            list_dataset_name = list_dataset_name,
                            list_dataset_type = list_dataset_type,
                            stimulus_filter = stimulus_filter)

  return(factoextra::fviz_eig(mfa_object, addlabels = TRUE))

}