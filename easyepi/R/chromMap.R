library(ggplot2)
library(dplyr)
source('arrayTobed.R')
# chromMap----------------------------------------------------------------------
#' Chromosome Mapping and Correlation Plotting Function
#'
#' This function takes a data frame, an value vector, and a specified format and chip type, then performs 
#' various transformations and calculations to map chromosome data, calculate correlations, and create 
#' a color-coded scatter plot of the correlations.
#'
#' @param df A data frame containing the data to be processed. It should be bed format, which have columns including 'chr' 
#' (chromosome information), 'start' (start position) and 'end' (end position). Input df in bed format. If the rowname of df 
#' is the probe ID, the format needs to be set to 'array'
#' @param value A numeric vector representing the value data for correlation calculation.
#' @param format A character string indicating the format of the input data. Supported formats are 'bed' and 'array'.
#' @param chip A character string specifying the methylation chip type. Supported chips are '40K', '450K', and '850K'.
#' @threshold A threshold of value, if set, will indicate the corresponding label on the graph, which is the row name of df.
#' @xlab The xlab of x axis
#' @ylab The ylab of y axis
#' @export

chromMap <- function(df, value, format = 'bed', chip, threshold, xlab = 'Chromosome', ylab = 'value'){
  if(format == 'array'){
    if (missing(chip) || is.null(chip)) {
      stop("pleas provide the chip: '40K', '450K' or '850K'.")
    }
    df <- arrayTobed(df, chip)
    df <- df[!grepl("chrX", df$chr), ]
    df <- df[!grepl("chrY", df$chr), ]
  }
  colnames(df)[1] <- 'chr'
  spec_value <- data.frame(lable = rownames(df), value = value, 
                        chr = gsub("[^0-9]", "", (df$chr)) %>% as.numeric(),
                        mapinfo = df$start+1)
  spec_value <- spec_value[order(spec_value$chr, spec_value$mapinfo),]
  
  spec_value$pos <- NA
  spec_value$index <- NA
  
  ind = 0
  
  for (i in unique(spec_value$chr)){
    ind =  ind +1 
    spec_value[spec_value$chr ==i, ]$index = ind
  }
  
  nchr <- length(unique(spec_value$chr))
  lastbase = 0
  ticks = NULL
  
  for (i in unique(spec_value$index)) {
    if (i == 1) {
      spec_value[spec_value$index == i, ]$pos = spec_value[spec_value$index == i, ]$mapinfo
    } else {
      lastbase = lastbase+tail(subset(spec_value, index == i-1)$mapinfo, 1)
      spec_value[spec_value$index == i, ]$pos = spec_value[spec_value$index == i, ]$mapinfo +lastbase
    }
    ticks <- c(ticks, (min(spec_value[spec_value$index == i,]$pos)+max(spec_value[spec_value$index == i, ]$pos))/2+1)
  }
  
  # plot color chrom map
  color_func <- colorRampPalette(c("blue",'green','orange', "red"))
  colors <- color_func(length(unique(spec_value$chr)))[spec_value$index]
  spec_value$color <- colors
  
  ymax <- round(max(abs(spec_value$value))+0.07, 1)

  p <- ggplot(spec_value, aes(x = pos, y = value, color = color)) +
    geom_point(size = 5, shape = 20) +
    scale_color_identity() +
    scale_x_continuous(
      breaks = ticks,           # 设置刻度位置
      labels = unique(spec_value$chr),            # 设置刻度标签
      minor_breaks = NULL      # 可选：设置次要刻度位置
    ) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, vjust = 7),
          axis.text.y = element_text(size = 15),
          # panel.grid.major.x = element_line(color = "grey", linetype = "solid"),
          panel.grid=element_blank()) +
    ylim(-ymax, ymax) +
    labs(x = xlab, y = ylab) +
    geom_hline(yintercept = 0, color = "black", linewidth = 2) +
    geom_hline(yintercept = c(round(ymax/3,1), -round(ymax/3,1), round(ymax*2/3,1), -round(ymax*2/3,1), ymax, -ymax), 
               color = "grey", linewidth = c(1, 1, 1, 1, 1, 1))
  if(!missing(threshold) && !is.null(threshold)){
    p <- p+
      ggrepel::geom_text_repel(data = subset(spec_value, abs(value) > threshold),
                               aes(label = lable), 
                               size = 6,
                               color = subset(spec_value, abs(value) > threshold)$color,
                               hjust = 1, vjust = 0)
  }
  print(p)
}
