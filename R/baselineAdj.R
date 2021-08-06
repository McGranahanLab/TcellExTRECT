#' Calculate median values of coverage within exons in rolling windows
#'
#' @param vdj_logR_input data frame of positions and coverage values
#' @param vdj_seg Segments used for normalisation
#' @param GCcorrect Use GC corrected output or not
#' @importFrom stats sd confint
#' @import gratia
#' @return Adjusted baseline logR dataframe
#' @name baselineAdj


baselineAdj <- function(vdj_logR_input, vdj_seg, GCcorrect = TRUE){

  pos <- reads <- NULL

  ratio.col <- ifelse(GCcorrect, 'Ratio.gc.correct','Ratio')
  ratio.col <- rlang::sym(ratio.col)

  adjust.baseline.value <- vdj_logR_input %>%
    dplyr::filter((pos >= vdj_seg[3,'start'] & pos <= vdj_seg[3,'end']) |
             (pos >= vdj_seg[4,'start'] & pos <= vdj_seg[4,'end'])) %>%
    dplyr::summarise(gc.adjust = mean(!!ratio.col),
              CI.95.range = 1.96*sd(!!ratio.col)/sqrt(length(!!ratio.col)))

  # Look at focal too in GAM model

  if(GCcorrect){
    adjust.model <- mgcv::gam(Ratio.gc.correct~s(pos, bs = 'cs'), data = vdj_logR_input)
  }else{
    adjust.model <- mgcv::gam(Ratio~s(pos, bs = 'cs'), data = vdj_logR_input)
  }


  adjust.fit.model <- confint(adjust.model, parm = "s(pos)",
                                           partial_match = TRUE, type = 'simultaneous',
                                           newdata = seq(vdj_seg[2,2], vdj_seg[2,3],by=100),
                                           shift = TRUE)
  fit.loc <- which(adjust.fit.model$pos > vdj_seg[2,2] & adjust.fit.model$pos < vdj_seg[2,3])
  adjust.baseline.value2 <- list(mean(adjust.fit.model$est[fit.loc]),
                                 mean(adjust.fit.model$upper[fit.loc]) - mean(adjust.fit.model$est[fit.loc]))

  adjust.value <- max(adjust.baseline.value[[1]], adjust.baseline.value2[[1]])
  ci.95.value <- ifelse(adjust.baseline.value2[[1]] > adjust.baseline.value[[1]],
                        adjust.baseline.value2[[2]], adjust.baseline.value[[2]])


  vdj_logR_input[[ratio.col]] <- vdj_logR_input[[ratio.col]] - adjust.value
  return(list(vdj_logR_input, ci.95.value))
}

