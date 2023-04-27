produce_boxplot_colored_by_stimulus <- function(data,
                                                color_scheme = NULL,
                                                show_outlier = TRUE,
                                                show_pvalue = TRUE,
                                                stimulus_target = NULL,
                                                fdr_adjusted = TRUE,
                                                paired_test = TRUE,
                                                log_fun = TRUE,
                                                font_size = 9,
                                                compound_order = NULL,
                                                zebra_background = FALSE,
                                                step_zebra = 1,
                                                y_title = 'analyte_concentration'){



  #### Helper functions
  boxplot_add_zerbra_background <- function(fig,
                                            data_table,
                                            step = 1.25,
                                            y_correction = -0.5,
                                            min_dat = NULL,
                                            max_dat = NULL,
                                            color_dark = 'rgba(210,210,210,0.3)',
                                            color_light = 'rgba(255,255,255,0)'){

    logp1 <- function(x){return(log(x + 1))}

    max_dat <- {if(is.null(max_dat)) data_table %>% dplyr::select(-dplyr::all_of(c("stimulus_id","donor_id"))) %>% dplyr::summarise_all(.funs = max) %>% .[1,] %>% max() %>% logp1() else max_dat}
    min_dat <- {if(is.null(min_dat)) data_table %>% dplyr::select(-dplyr::all_of(c("stimulus_id","donor_id"))) %>% dplyr::summarise_all(.funs = min) %>% .[1,] %>% min() %>% logp1() else min_dat}
    n_col_to_color <- colnames(data_table) %>% setdiff(c("stimulus_id","donor_id")) %>% length()

    list_rect <- list()
    rect_type <- list(type = 'rect',
                      xref = 'x',
                      yref='y',
                      line = list(color = 'rgba(0,0,0,0)'), layer = 'below')

    for(i in 0:(n_col_to_color - 1) ){
      tmp_rect <- rect_type
      tmp_rect$'x0' <- min_dat - 0.1#(i - 1) * step
      tmp_rect$'x1' <- max_dat + 0.1
      tmp_rect$'y0' <- i * step + y_correction
      tmp_rect$'y1' <- (i + 1) * step + y_correction
      tmp_rect$fillcolor <- {if(i %% 2 == 0) color_light else color_dark}
      list_rect[[i+1]] <- tmp_rect
    }

    #return(list_rect)
    fig <- fig %>% plotly::layout(fig, shapes= list_rect)

    return(fig)
  }


  produce_order_array <- function(list_compound, list_stimulus){

    ret <- list()
    for(comp in setdiff(list_compound, c('stimulus_id','donor_id'))){
      ret[[comp]] <- as.list(list_stimulus)
    }

    return(ret)

  }



  list_compound <- NULL
  if(! is.null(compound_order)){
    list_compound <- compound_order
  } else {
    list_compound <- setdiff(colnames(data), c('stimulus_id','donor_id'))
  }

  data_plot <- data %>% dplyr::mutate_at(.vars = dplyr::vars(dplyr::all_of(c('stimulus_id','donor_id'))), .funs = as.factor) %>%
    dplyr::mutate_at(.vars = dplyr::vars(- dplyr::all_of(c('stimulus_id','donor_id'))), .funs = as.numeric) %>%
    dplyr::filter(stimulus_id %in% stimulus_target) %>%
    reshape2::melt(id = c('donor_id','stimulus_id')) %>%
    {if(is.function(log_fun)) dplyr::mutate_at(.tbl = ., .vars = dplyr::vars('value'), .funs = log_fun)
      else if(log_fun) dplyr::mutate_at(.tbl  = ., dplyr::vars('value'), .funs = log)
        else . } %>%
    {if(show_pvalue) dplyr::mutate(., donor_id = paste0('donor_',donor_id)) %>% dplyr::mutate(suffixe = " ") else .}

  if(show_pvalue && length(unique(data_plot$stimulus_id)) == 2){

    pvalue_dat <- tribble(~compound_stat,~pvalue)
    for(compound in list_compound){

      stim_1 <- data %>% dplyr::filter(stimulus_id == stimulus_target[[1]]) %>%
        dplyr::select(dplyr::all_of(c(compound,'donor_id'))) %>% dplyr::rename(stim_1 = 1)
      stim_2 <- data %>% dplyr::filter(stimulus_id == stimulus_target[[2]]) %>%
        dplyr::select(dplyr::all_of(c(compound,'donor_id'))) %>% dplyr::rename(stim_2 = 1)

      if(paired_test){
        stim_dat <- dplyr::inner_join(stim_1, stim_2, by = 'donor_id')
        wilco <- wilcox.test(x = as.numeric(stim_dat$stim_1), y = as.numeric(stim_dat$stim_2), paired = paired_test, exact = FALSE)
        pvalue_dat <- pvalue_dat %>% dplyr::add_row(compound_stat = compound, pvalue = wilco$p.value)
      } else  {
        wilco <- wilcox.test(x = as.numeric(stim_1$stim_1), y = as.numeric(stim_2$stim_2), paired = paired_test, exact = FALSE)
        pvalue_dat <- pvalue_dat %>% dplyr::add_row(compound_stat = compound, pvalue = wilco$p.value)
      }
    }

    if(fdr_adjusted){
      pvalue_dat <-  pvalue_dat %>% dplyr::mutate(pvalue = p.adjust(pvalue, method = "fdr"))
    }

    for(compound in list_compound){
      pvalue <- pvalue_dat %>% dplyr::filter(compound_stat == compound) %>% dplyr::pull(pvalue)
      pvalue <-  pvalue[[1]]

      if(pvalue > 0.05){
        data_plot["suffixe"][data_plot['variable'] == compound] = "         "
      } else if(pvalue <= 0.05 && pvalue > 0.01 ){
        data_plot["suffixe"][data_plot['variable'] == compound] <- "     (*)"
      } else if(pvalue <= 0.01 && pvalue >= 0.001 ){
        data_plot["suffixe"][data_plot['variable'] == compound] <- "    (**)"
      } else if(pvalue <= 0.001 && pvalue >= 0.0001 ){
        data_plot["suffixe"][data_plot['variable'] == compound] <- "   (***)"
      } else {
        data_plot["suffixe"][data_plot['variable'] == compound] <- " (****)"
      }

    }
  }


  data_plot <- data_plot %>%  dplyr::mutate(variable = paste0(variable, " ", suffixe))

  fig <- NULL
  for(stimulus in unique(data_plot$stimulus_id)){

    data_plot_temp <- data_plot %>% dplyr::filter(stimulus_id == stimulus)

    if(is.null(fig)){

      fig <- plotly::plot_ly(data = data_plot_temp,
                     type = 'box',
                     y = ~variable,
                     x = ~value,
                     #color = ~stimulus_id,
                     name = ~stimulus_id,
                     marker = list(color = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]]),
                     line = list(color = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]]),
                     fillcolor = list(color = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]]),
                     #color = {if(is.null(color_scheme)) NULL else color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]]},
                     text = ~donor_id,
                     #categoryorders = produce_order_array(data_plot_temp$variable, stimulus_target),
                     #categoryorder = 'trace',
                     boxpoints = show_outlier ) %>%
        plotly::layout(boxmode = "group",
               xaxis = list( zeroline=FALSE, title = list(text = 'log(x+1)')),# categoryarray = unique(data_plot_temp$variable), categoryorder = "array"),
               yaxis = list( zeroline=FALSE, title = list(text = NULL),  categoryorder = "array", categoryarray = unique(data_plot_temp$variable))) #categoryarray = produce_order_array(data_plot_temp$variable, stimulus_target), categoryorder = "array"))

    }else {

      fig <- fig %>% plotly::add_trace(data = data_plot_temp,
                                       type = 'box',
                                       y = ~variable,
                                       x = ~value,
                                       #color = ~stimulus_id,
                                       name = ~stimulus_id,
                                       line = list(color = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]]),
                                       marker = list(color = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]]),
                                       fillcolor = list(color = color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]]),
                                       #color = {if(is.null(color_scheme)) NULL else color_scheme[color_scheme$stimulus_id == stimulus,]$color_hex[[1]]},
                                       #categoryorders = produce_order_array(data_plot_temp$variable, stimulus_target),
                                       #categoryorder = 'trace',
                                       text = ~donor_id,
                                       boxpoints = show_outlier ) %>%
        plotly::layout(boxmode = "group",
               font = list(family = "Helvetica", size = font_size),
               #legend = list(traceorder = 'normal'),
               #font = list(family = "sans serif"),
               xaxis = list( zeroline=FALSE, title = list(text = y_title), categoryarray = stimulus_target, categoryorder = "array"),
               yaxis = list( zeroline=FALSE, title = list(text = NULL), type='category', categoryorder = "array", categoryarray = unique(data_plot_temp$variable)))# categoryarray = produce_order_array(data_plot_temp$variable, stimulus_target), categoryorder = "array"))

    }
  }

  if(zebra_background){
    fig <-  fig %>% boxplot_add_zerbra_background(data_table = data, step = step_zebra, min_dat = min(data_plot$value, na.rm = TRUE), max_dat = max(data_plot$value, na.rm = TRUE))
  }

  return(fig)

}



