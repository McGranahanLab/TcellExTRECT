#' Method for calculating QC values to judge GAM fit
#'
#' @param tumour.logR LogR values
#' @param segs GC content of exons
#' @param ci_type FASTA file for VDJ gene e.g TCRA
#' @param GC_correct whether to use GC correction or not for output
#' @return QC values for judging fit of GAM model
#' @name calcQCvalue
#' @export

calcQCvalue <- function(tumour.logR, segs, GC_correct = TRUE, ci_type = 'simultaneous'){

  if(GC_correct){
    gam.model <- mgcv::gam(Ratio.gc.correct~s(pos, bs = 'cs'), data = tumour.logR)
  }else{
    gam.model <- mgcv::gam(Ratio~s(pos, bs = 'cs'), data = tumour.logR)
  }

  fit.model <- confint(gam.model, parm = "s(pos)",
                       partial_match = TRUE, type = ci_type,
                       data = data.frame(pos = seq(segs[1,2],segs[1,3],by=100)),
                       shift = TRUE)
  xz <- zoo::as.zoo(fit.model$est)

  rx.max <-zoo::rollapply(xz, 3, function(x) which.max(x)==2)
  rx.max.loc <- zoo::index(rx.max)[zoo::coredata(rx.max)]

  rx.min <- zoo::rollapply(xz, 3, function(x) which.min(x)==2)
  rx.min.loc <- zoo::index(rx.min)[zoo::coredata(rx.min)]

  all.minmax <- abs(c(fit.model$est[rx.min.loc],fit.model$est[rx.max.loc]))

  all.min.norm <- fit.model$est[rx.min.loc]/max(all.minmax)
  all.max.norm <- fit.model$est[rx.max.loc]/max(all.minmax)

  return(c(max(all.minmax),
           sum(abs(c(all.min.norm, all.max.norm)))))

}
