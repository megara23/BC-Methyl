##################################################
# helper functions
##################################################

########################################
# summarize a single result, with
#     1 = FC > 1, significant
#     2 = FC > 1, ns
#     3 = FC < 1, ns
#     4 = FC < 1, significant
########################################
summarizeOne <- function(x) {
  if (is.null(x)) {
    return(NULL)
  } else if (is.na(x$FC) | is.na(x$FDR)) {
    return(NULL)
  }
  if (x$FC > 1 & x$FDR <= 0.05) {
    return(1)
  } else if (x$FC > 1 & x$FDR > 0.05) {
    return(2)
  } else if (x$FC < 1 & x$FDR > 0.05) {
    return(3)
  } else if (x$FC < 1 & x$FDR <= 0.05) {
    return(4)
  } 
  return(NULL)
}

####################################
# get group labels from data frame
####################################
label.value <-function(x){
  r = x$group
  r[x$value==0]=NA
  r = gsub(",", ",\n", r)
  r
}

##########################################################
# Main plotting function -- plots summary based on 
#   the list X that includes FC and FDR for each dataset
##########################################################

plotSummary <- function(X) {

    # get summary vector
  s = sapply(X, summarizeOne)
  
  # levels corresponding to 1-4
  levels = c("methylated, (adjusted) P<0.05)", "methylated, (adjusted) P>0.05)",
           "demethylated, (adjusted) P>0.05)", "demethylated, (adjusted) P<0.05)")

  colors = c("#FF0000", "#FFC0CB", "#00FFFF", "#0000FF")
  names(colors) = levels

  f = factor(s, levels = 1:4)
  t = table(f)
  df = data.frame(group = factor(levels, levels = levels), value = as.integer(t))
  df = df[4:1,]

  # Barplot
  bp<- ggplot(df, aes(x="", y=value, fill=group))+
    ggtitle("Results Summary")+
    theme(axis.text.x=element_blank(), 
        legend.position = "none",
        axis.title.y = element_text(size = 15, face = "bold", vjust = 2),
        plot.title = element_text(size = 18, face = "bold", vjust = 2),
        plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
    xlab("") + ylab("# results") +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values=colors) + ylim(0,4) +
    geom_label(aes(y = cumsum(df$value) - diff(c(0,cumsum(df$value))) / 2,
        label = label.value(df)), size=4, fill = "lightgray") 

    print(bp)
}