function (cross, phenotype.name, focal.groups = NULL, nuisance.groups = NULL, 
          genotype.names = c("AA", "AB", "BB"), xlim = NULL, ylim = NULL, 
          title = paste(phenotype.name, "by", paste(focal.groups, collapse = ", ")), 
          draw_ribbons = TRUE, se_line_size = 1, point_size = 1) 
{
  indiv.mean.estim <- indiv.mean.lb <- indiv.mean.ub <- "fake_global_for_CRAN"
  indiv.sd.estim <- indiv.sd.lb <- indiv.sd.ub <- "fake_global_for_CRAN"
  group.mean.estim <- group.mean.ub <- group.mean.lb <- "fake_global_for_CRAN"
  group.sd.estim <- group.sd.ub <- group.sd.lb <- "fake_global_for_CRAN"
  validate.mean_var_plot_model_based.input(cross = cross, phenotype.name = phenotype.name, 
                                           focal.groups = focal.groups, nuisance.groups = nuisance.groups)
  modeling.df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))
  modeling.df[[phenotype.name]] <- cross[["pheno"]][[phenotype.name]]
  marker.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.geno(cross = cross))], 
                    nuisance.groups[nuisance.groups %in% colnames(qtl::pull.geno(cross = cross))])
  phen.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.pheno(cross = cross))], 
                  nuisance.groups[nuisance.groups %in% colnames(qtl::pull.pheno(cross = cross))])
  for (marker.name in marker.names) {
    modeling.df[[marker.name]] <- factor(x = qtl::pull.geno(cross = cross)[, 
                                                                           marker.name], labels = genotype.names)
  }
  for (phen.name in phen.names) {
    modeling.df[[phen.name]] <- factor(qtl::pull.pheno(cross = cross)[[phen.name]])
  }
  modeling.df[["placeholder"]] <- NULL
  covar.form <- paste(focal.groups, collapse = "+")
  if (!is.null(nuisance.groups)) {
    covar.form <- paste(covar.form, "+", paste(nuisance.groups, 
                                               collapse = "+"))
  }
  mean.form <- paste(phenotype.name, "~", covar.form)
  var.form <- paste("~", covar.form)
  dglm.fit <- dglm::dglm(formula = stats::formula(mean.form), 
                         dformula = stats::formula(var.form), data = modeling.df)
  mean.pred <- stats::predict(dglm.fit, se.fit = TRUE)
  mean.estim <- mean.pred$fit
  mean.se <- mean.pred$se.fit
  sd.pred <- stats::predict(dglm.fit$dispersion.fit, se.fit = TRUE, 
                            dispersion = 2)
  sd.estim <- sd.pred$fit/(sd.pred$residual.scale^2)
  sd.se <- sd.pred$se.fit
  indiv.prediction.tbl <- dplyr::bind_cols(stats::na.omit(modeling.df), 
                                           dplyr::data_frame(indiv.mean.estim = mean.estim, indiv.mean.lb = mean.estim - 
                                                               mean.se, indiv.mean.ub = mean.estim + mean.se, indiv.sd.estim = exp(sd.estim), 
                                                             indiv.sd.lb = exp(sd.estim - sd.se), indiv.sd.ub = exp(sd.estim + 
                                                                                                                      sd.se)))
  group.prediction.tbl <- indiv.prediction.tbl %>% dplyr::group_by_(.dots = c(focal.groups)) %>% 
    dplyr::summarise(group.mean.estim = mean(indiv.mean.estim), 
                     group.mean.lb = mean(indiv.mean.lb), group.mean.ub = mean(indiv.mean.ub), 
                     group.sd.estim = mean(indiv.sd.estim), group.sd.lb = mean(indiv.sd.lb), 
                     group.sd.ub = mean(indiv.sd.ub))
  p <- ggplot2::ggplot(data = group.prediction.tbl, mapping = ggplot2::aes_string(color = focal.groups[1]))
  if (draw_ribbons & length(focal.groups) > 1) {
    p <- p + ggplot2::geom_path(data = group.prediction.tbl, 
                                mapping = ggplot2::aes_string(x = "group.mean.estim", 
                                                              y = "group.sd.estim", color = focal.groups[1]), 
                                size = 4, alpha = 0.3)
  }
  p <- p + ggplot2::geom_segment(mapping = ggplot2::aes(x = group.mean.lb, 
                                                        xend = group.mean.ub, y = group.sd.estim, yend = group.sd.estim), 
                                 size = se_line_size) + ggplot2::geom_segment(mapping = ggplot2::aes(x = group.mean.estim, 
                                                                                                     xend = group.mean.estim, y = group.sd.lb, yend = group.sd.ub), 
                                                                              size = se_line_size)
  p <- p + ggplot2::theme_minimal() + ggplot2::xlab("mean estimate +/- 1 SE") + 
    ggplot2::ylab("SD estimate +/- 1 SE") + ggplot2::ggtitle(title) + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
    ggplot2::guides(ggplot2::guide_legend(order = 1))
  if (length(focal.groups) == 1) {
    p <- p + ggplot2::geom_point(mapping = ggplot2::aes(x = group.mean.estim, 
                                                        y = group.sd.estim), size = point_size)
  }
  if (length(focal.groups) > 1) {
    p <- p + ggplot2::geom_point(mapping = ggplot2::aes_string(x = "group.mean.estim", 
                                                               y = "group.sd.estim", shape = focal.groups[2]), size = point_size)
  }
  if (!is.null(xlim) & !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  }
  if (!is.null(xlim) & is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }
  if (is.null(xlim) & !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = ylim)
  }
  return(p)
}
<environment: namespace:vqtl>