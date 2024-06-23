library(ggplot2)

#' draw the percentile plot
#' 
#' draw the percentile plot of target and background samples
#' 
#' @param target the values of target samples. -1 will be removed.
#' @param background the values of background samples. -1 will bed removed.
#' @param is.hypo If the merged methylation level of target samples is less than background samples.
#' @param p the percentile of low methylation group. If is.hypo is true, the p-th percentile is caclulated
#' @param ylab the ylab
#' @param title the title
#' @return None.
#' @examples
#' figPercentileMethy
#' @export


figPercentileMethy<-function(target, background, is.hypo=TRUE,p=0.8, ylab="",title="") {
  target<-unlist(target)
  background<-unlist(background)
  data<-data.frame(values=c(target, background),group=c(rep('Target',length(target)),rep('Backgound',length(background))))
  data<-data[data$values >=0,]
  if (is.hypo){
    percentile.low <- quantile(data$values[data$group == "Target"], p)
    percentile.high <- quantile(data$values[data$group == "Backgound"], 1-p)
    title <- paste(title, "hypo")
  }else{
    percentile.low <- quantile(data$values[data$group == "Backgound"], p)
    percentile.high <- quantile(data$values[data$group == "Target"], 1-p)
    title <- paste(title, "hyper")
  }
  set.seed(123)
  ggplot(data, aes(x = group, y = values)) +
    geom_boxplot(outlier.shape = NA, fill = NA, color = NA) +
    geom_jitter(width = 0.1, size = 2) +
    geom_hline(yintercept = percentile.high, linetype = "dashed", color = "red3") +
    annotate("text", x = 2.5, y = percentile.high, label = percentile.high, color = "red3") +
    geom_hline(yintercept = percentile.low, linetype = "dashed", color = "green3") +
    annotate("text", x = 2.5, y = percentile.low, label = percentile.low, color = "green3") +
    coord_flip() +
    labs(x = ylab, y = "Mehtylation Level", title = title) +
    theme_classic() +
    theme(axis.text.y = element_text(angle = 90))
}
