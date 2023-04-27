
visualize_as_pca <- function(data,
                             meta_data = NULL,
                             compo_x = 1,
                             compo_y = 2,
                             compo_z = NULL,
                             flip_x = FALSE,
                             flip_y = FALSE,
                             title_x = paste0('PC ', compo_x),
                             title_y = paste0('PC ', compo_y),
                             title_explained_variance = TRUE,
                             stimulus_filter = 'ALL',
                             color_palette = NULL,
                             color = NULL,
                             symbol = NULL,
                             log_fun = function(x){return(log(x+1))},
                             center = TRUE,
                             scale = FALSE,
                             size = 10,
                             marker_opacity = 0.75,
                             n_vect_loading_plot = 5,
                             fact_scale_loading_vect = 3,
                             loading_plot = FALSE,
                             margin_plot = FALSE){

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


    data <- data %>% dplyr::mutate_at(.vars = dplyr::vars(dplyr::all_of(c('stimulus_id','donor_id'))), .funs = as.factor) %>%
    dplyr::mutate_at(.vars = dplyr::vars(- dplyr::all_of(c('stimulus_id','donor_id'))), .funs = as.numeric) %>%
    remove_non_variating_column(other_variable_to_keep = c('stimulus_id','donor_id')) %>%
    {if(stimulus_filter[[1]] != 'ALL') dplyr::filter(., stimulus_id %in% stimulus_filter) else . }

    if(! is.null(meta_data)){
      meta_data <- meta_data %>% dplyr::mutate_at(.vars = dplyr::vars(dplyr::any_of(c('stimulus_id','donor_id'))), .funs = as.factor) %>%
        {if(stimulus_filter[[1]] != 'ALL' && ! is.null(meta_data[['stimulus_id']])) dplyr::filter(., stimulus_id %in% stimulus_filter) else . }
    } else {
      meta_data <- dplyr::select(data, dplyr::all_of(c('donor_id','stimulus_id')))
    }

  # If an argument to color the PCA is given in parameter it should be
  # either a dimension of the data used for the PCA or diretly provided in the meta_data argument
  if(!is.null(color) && !(color %in% colnames(meta_data)  || (color %in% colnames(data)))){
    warning(paste0(color," is neither present in data or metadata, color will be ignored"), classes = "aesthetic_warnng::color")
  }

  # Same logic for symbol : If an argument to color the PCA is given in parameter it should be
  # either a dimension of the data used for the PCA or diretly provided in the meta_data argument
  if(! is.null(symbol) && ! (symbol %in% colnames(meta_data)) ){
    warning(paste0(symbol," is neither present in data or metadata, color will be ignored"), classes = "aesthetic_warnng::symbol")
  }

  # If there is duplicated column name in data and metadata they are substracted from metadata
  common_column <- intersect(colnames(dplyr::select(data, - all_of(c('stimulus_id','donor_id')))),
                             colnames(dplyr::select(meta_data, - any_of(c('stimulus_id','donor_id')))))
  if(length(common_column) > 0){
    warning(common_column, "were both found in data and meta_data, only the one in data will be kept",  classes = "data_warnng::duplicated_columns")
    meta_data <- dplyr::select(meta_data, - all_of(common_column))
  }

  # generate components
  prin_comp <- data %>%
    dplyr::select(-dplyr::any_of(c('donor_id', 'stimulus_id', 'sex', 'cmv'))) %>%
    dplyr::mutate_all(as.numeric) %>%
    {if(is.function(log_fun)) as.data.frame(lapply(X = ., FUN = log_fun))
        else if (log_fun) log(.) else . } %>%
    prcomp(center = center, scale = scale)

  # explained variance
  explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
  explained_variance_ratio <- 100 * explained_variance_ratio

  if(title_explained_variance){
    title_x <- paste0(title_x, '  (', round(explained_variance_ratio[compo_x], digits = 2), '% var. explained)')
    title_y <- paste0(title_y, '  (', round(explained_variance_ratio[compo_y], digits = 2), '% var. explained)')
  }

  all_data <- {if(! is.null(meta_data$stimulus_id)) dplyr::inner_join(data, meta_data, by = c('donor_id','stimulus_id'))
                  else dplyr::inner_join(data, meta_data, by = 'donor_id') } %>% dplyr::mutate(size = as.numeric(1))

  pca_components <- prin_comp[["x"]] %>% as.data.frame()
  pca_components <- dplyr::bind_cols(pca_components, dplyr::select(data, dplyr::all_of(c('stimulus_id','donor_id'))))

  n_compo_x <- compo_x
  compo_x <- paste0('PC',compo_x)
  n_compo_y <- compo_y
  compo_y <- paste0('PC',compo_y)

  pca_components <- pca_components %>%
    {if(flip_x) dplyr::mutate_at(., compo_x, .funs = function(x){-x}) else .} %>%
    {if(flip_y) dplyr::mutate_at(., compo_y, .funs = function(x){-x}) else .}

  ## We produce the scatter plot
  fig <- produce_scatter_plot_colored_by_stimulus(data = pca_components,
                                                  color_scheme = color_palette,
                                                  x = compo_x,
                                                  y = compo_y,
                                                  z = {if(! is.null(compo_z)) paste0('PC',compo_z) else NULL},
                                                  marker_opacity = marker_opacity,
                                                  marker_size = size) %>%
    plotly::layout(xaxis = list(title = title_x, showticklabels=FALSE),
                   yaxis = list(title = title_y, showticklabels=FALSE))

  ## If asked we had loading plot on it
  if(loading_plot){
    fig <- fig %>% add_loading_plot(fig = .,
                                    prin_comp = prin_comp,
                                    n_vect_loading_plot = n_vect_loading_plot,
                                    fact_scale_loading_vect = fact_scale_loading_vect,
                                    flip_x = flip_x,
                                    flip_y = flip_y,
                                    compo_x = n_compo_x,
                                    compo_y = n_compo_y)
  }

  # if asked we add margin plot
  if(is.null(compo_z) && margin_plot){
    fig <- fig %>% add_margin_plot(fig = .,
                                   data = pca_components,
                                   x = compo_x,
                                   y = compo_y,
                                   color_scheme = color_palette)
  }



  return(fig)

}




produce_scatter_plot_colored_by_stimulus <- function(data,
                                                     color_scheme,
                                                     x,
                                                     y,
                                                     z = NULL,
                                                     marker_opacity = 1,
                                                     marker_size = 15){


  fig <- NULL
  for(stimulus in unique(data$stimulus_id)){

      type_graph <- {if(! is.null(z)) "scatter3d" else "scatter"}
      data_stim <- dplyr::filter(data, stimulus_id == stimulus)

      if(is.null(fig)){
         fig <- plotly::plot_ly(type = type_graph,
                                mode = "markers",
                                #x = {if(flip_x) lapply(X = data_stim$compo_x, FUN = function(x){-x}) else data_stim$compo_x},
                                #y = {if(flip_y) lapply(X = data_stim$compo_y, FUN = function(x){-x}) else data_stim$compo_y},
                                x = dplyr::pull(data_stim, x),
                                y = dplyr::pull(data_stim, y),
                                z = {if(! is.null(z)) dplyr::pull(data_stim, z)},
                                name = stimulus,
                                color = stimulus,
                                text = paste0(dplyr::pull(data_stim,'donor_id'),'_',stimulus),
                                legendgroup = stimulus,
                                size = rep(marker_size, time(nrow(data))),
                                sizes = c(marker_size, marker_size),
                                 marker = list(color = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]],
                                               opacity = marker_opacity,
                                               line = list(color = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]])))

      } else {
        fig <- fig %>% add_trace(#x = {if(flip_x) lapply(X = data_stim$compo_x, FUN = function(x){-x}) else data_stim$compo_x},
                                 #y = {if(flip_y) lapply(X = data_stim$compo_y, FUN = function(x){-x}) else data_stim$compo_y},
                                 x = dplyr::pull(data_stim, x),
                                 y = dplyr::pull(data_stim, y),
                                 z = {if(! is.null(z)) dplyr::pull(data_stim, z)},
                                 name = stimulus,
                                 text = paste0(dplyr::pull(data_stim,'donor_id'),'_',stimulus),
                                 type = type_graph,
                                 mode = 'markers',
                                 color = stimulus,
                                 size = rep(marker_size, time(nrow(data))),
                                 sizes = c(marker_size,marker_size),
                                 opacity = marker_opacity,
                                 legendgroup = stimulus,
                                 marker = list(color = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]],
                                               opacity = marker_opacity,
                                               line = list(color = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]])))
    }
  }

  return(fig)
}


#' @param fig the plotly scatterplot figure
#' @param data the point coordinate
#'
add_margin_plot <- function(fig,
                            data,
                            x,
                            y,
                            flip_x = FALSE,
                            flip_y = FALSE,
                            color_scheme){

    top_subplot <- NULL
    side_subplot <- NULL

    for(stimulus in unique(data$stimulus_id)){

      data_stim <- dplyr::filter(data, stimulus_id == stimulus) %>%
      {if(flip_x) dplyr::mutate_at(., x, .funs = function(x){-x}) else .} %>%
      {if(flip_y) dplyr::mutate_at(., y, .funs = function(x){-x}) else .}

      if(is.null(top_subplot) || is.null(side_subplot)){
        top_subplot <- plot_ly(x = dplyr::pull(.data = data_stim, var = x),
                               type = "box",
                               colors = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]],
                               color = stimulus, boxpoints = FALSE, legendgroup = stimulus) %>%
          style(showlegend = FALSE)
        side_subplot <- plot_ly(y = dplyr::pull(.data = data_stim, var = y),
                                type = "box",
                                colors = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]],
                                color = stimulus,
                                boxpoints = FALSE, legendgroup = stimulus) %>%
          style(showlegend = FALSE)

      } else {

        top_subplot <-  top_subplot %>%
          add_trace(x =  dplyr::pull(.data = data_stim, var = x),
                    type = "box",
                    colors = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]],
                    color = stimulus,
                    boxpoints = FALSE,
                    legendgroup = stimulus) %>%
          style(showlegend = FALSE)

        side_subplot <-  side_subplot %>%
          add_trace(y = dplyr::pull(.data = data_stim, var = y),
                    type = "box",
                    colors = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]],
                    color = stimulus,
                    boxpoints = FALSE,
                    legendgroup = stimulus) %>%
          style(showlegend = FALSE)
      }

    }

  # merge everything together
  fig <- fig %>% style(showlegend = TRUE)
  fig <- plotly::subplot(top_subplot,
                          plotly::plotly_empty(),
                          fig,
                          side_subplot,
                          nrows = 2, heights = c(.05, .95),
                          widths = c(.95,.05),
                          margin = 0,
                          shareX = TRUE,
                          shareY = TRUE)


  return(fig)

}



#' @param fig the plotly scatterplot figure
#' @param prin_comp object generated using the prcomp function
#' @param n_vect_loading_plot n top vector to show
add_loading_plot <- function(fig,
                             prin_comp,
                             compo_x,
                             compo_y,
                             flip_x = FALSE,
                             flip_y = FALSE,
                             n_vect_loading_plot = 5,
                             fact_scale_loading_vect = 3,
                             show_vector = TRUE ){




  # helper function to compute euclidian distance
  top_euclidean_dist <-  function(data, compo_x, compo_y, n_top){
    data <- data %>% as.data.frame() %>%
      dplyr::select(dplyr::all_of(c(compo_x,compo_y)))
    data$dist <- sqrt(data[,1] ^2 + data[,2] ^2)
    data <- dplyr::top_n(data, n = n_top, wt = dist)
    return(data)
  }

  explained_variance <- summary(prin_comp)[["sdev"]]
  explained_variance <- explained_variance[c(compo_x,compo_y)]

  comp <- prin_comp[["rotation"]]
  loadings <- prin_comp[["rotation"]]

  for (i in seq(explained_variance)){
    loadings[,i] <- comp[,i] * explained_variance[i] * fact_scale_loading_vect
  }

  loadings <- loadings %>%
    {if(flip_x) dplyr::mutate_at(., compo_x, .funs = function(x){-x}) else .} %>%
    {if(flip_y) dplyr::mutate_at(., compo_y, .funs = function(x){-x}) else .}

  loadings <- top_euclidean_dist(data = loadings,
                                 compo_x = compo_x,
                                 compo_y = compo_y,
                                 n_top = min(nrow(loadings),n_vect_loading_plot))

  features <- rownames(loadings)
  for( i in 1:min(nrow(loadings),n_vect_loading_plot)){
    fig <- fig %>%
      {if(show_vector) add_segments(., x = 0, xend = loadings[i, 1], y = 0, yend = loadings[i, 2], line = list(color = 'rgba(120, 120, 145, 0.71)'),inherit = FALSE) else .} %>%
      add_annotations(x=loadings[i, 1], y=loadings[i, 2], ax = 0, ay = 0,text = features[i], xanchor = 'center', yanchor= 'bottom') %>%
          style(showlegend = FALSE)
  }

  return(fig)

}




