rm(list=ls())
library(survival)
library(ggplot2)
library(R.utils);

#' My default preferences for ggplots
#' @export
#'
#' @description BW / remove major/minor panels, bold lines, etc
#'
myggplotdefaults <- function() {
  p <- theme_bw() + theme(	panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.background = element_blank(),
                           panel.border = element_blank(),
                           axis.line.x = element_line(colour = "black", size=1.0),
                           axis.line.y = element_line(colour = "black", size=1.0),
                           title=element_text(face="bold", size=20),
                           axis.title=element_text(face="bold", size=18),
                           axis.title.x=element_text(vjust=0.3),
                           axis.text=element_text(face="bold", size=14),
                           legend.text=element_text(size=14)
  );
  return(p)
}

#' Survival Curves by a Variable. displays p-values by log-rank & cox model
#' @export
#'
#' @param event_time   Time to Event
#' @param event        Event occured or not (0 or 1)
#' @param Max          Maximum time to display on graph
#' @param Var          Variable to split Survival Curves by (should be categorical or ordinal)
#' @param lgnd         Legend text (optional)
#' @param tit         Title of graph  (default = 'Survival')
#'
#' @return P-value from log-rank test
#'
#' @seealso \code{\link{qkm}}, \code{\link{ggkm}}, \code{\link{qkm_var_med}}, \code{\link{qkm_var_anal}}
#'

qkm2 <- function(event_time, event, max, var, lgnd, tit="Survival") {
  lbls = c("1","2","3","4","5","6","7","8","9");
  so <- Surv(event_time, event)~var
  colors<-seq(1,9,1);
  mfit<-survfit(so)
  #print(mfit)
  #print(survdiff(so))
  plot(mfit, xlim=c(0,max), conf.int=FALSE, xlab="Years", ylab="Fraction", col=colors)
  if ( missing(lgnd) ) {
    lnum = length(unique(na.omit(var)));
    ltxt = lbls[1:lnum];
    legend("bottomright", ltxt, lty=1, col=colors)
  } else {
    legend("bottomright", lgnd, lty=1, col=colors)

  }
  m <- coxph(so)
  #print(summary(m))
  title(tit)

  # Get annotations to add to graph for survival curve
  annot <- km_get_pval(so);
  if ( is.null(annot) ) { return (NA ) }
  stext(annot$str1, pos=0, side=3);
  stext(annot$str2, pos=1,side=3);

  return(annot$pv);
}


# Utility funciton to get strings to annotate graphs
# returns pvalue, and a couple of string annotations for graphs (ggkm, qkm2)
# input is a survival function object
km_get_pval <- function(so) {
  tmp <- tryCatch({survdiff(so)},
                  error = function(e) {
                    warning("Error computing survival difference for group");
                  }
  );

  # some problem computing survival difference
  if ( is.list(tmp) == FALSE ) {
    return(NULL);
  }
  # add log-rank p-value to top of graph
  df <- length(tmp$obs)-1;
  pv <- 1 - pchisq(tmp$chisq, df);
  if ( pv < 0.01) {
    str <- sprintf("log-rank p=%.2e", pv);
  } else {
    str <- sprintf("log-rank p=%.2f", pv);
  }


  # add N in each group to bottom of graph
  str2 <- "";
  for (j in 1:length(tmp$n)) {
    str2 <- paste(str2, "N",j,"= ",tmp$n[j], " ",sep="");
  }

  annot <- list();
  annot$str1 <- str;
  annot$str2 <- str2
  annot$pv <- pv;
  return(annot);
}

ovdat <- read.csv("Supplementary_Files/Ovary_OS_KM_matrix.txt", sep="\t")

# Create a group bi-allelic pathogenic vs. all else
# also create an ID column
ovdat$group2 <- 0
ovdat$group2[ovdat$group == "bi_patho"] <- 1
ovdat$id2 <- paste("id",ovdat$ids, sep="")
rownames(ovdat) <- ovdat$id2


# qkm2(ovdat$death/365, ovdat$status, 5, ovdat$group)
pdf("Results_Figures_and_P_Values/Supplementary_Fig3_OS_Ovary.pdf",width=7,height=5)
qkm2(ovdat$death/365, ovdat$status, 5, ovdat$group2, c("Wild-type","Bi-allelic Pathogenic"))
dev.off()

# Read in ovarian data
clindat <- read.csv("Supplementary_Files/ovary.csv", stringsAsFactors = FALSE)
clindat$samp_type <- substr(clindat$SAMPLE.ID, 14,15)
# remove smaples that are not 01
clindat <- subset(clindat, samp_type=="01")
clindat$id2 <- paste("id", substr(clindat$PATIENT.ID, 9,12), sep="")
rownames(clindat) <- clindat$id2


combdat <- merge(ovdat, clindat, by=0)
combdat$tmpStage <- combdat$Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Group.Stage
combdat$new_stage <- NA
i <- which(combdat$tmpStage %in% c("Stage IA", "Stage IB", "Stage IC"))
combdat$new_stage[i] <- 1
i <- which(combdat$tmpStage %in% c("Stage IIA", "Stage IIB", "Stage IIC"))
combdat$new_stage[i] <- 2
i <- which(combdat$tmpStage %in% c("Stage IIIA", "Stage IIIB", "Stage IIIC"))
combdat$new_stage[i] <- 3
i <- which(combdat$tmpStage %in% c("Stage IV"))
combdat$new_stage[i] <- 4
	
# Group Stage I-II vs III/IV
combdat$new_stage2 <- 1
i <- which(combdat$new_stage > 2)
combdat$new_stage2[i] <- 2

# qkm2(combdat$death/365, combdat$status, 5, combdat$group2)
# qkm2(combdat$death/365, combdat$status, 5, combdat$Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Group.Stage)
# qkm2(combdat$death/365, combdat$status, 5, combdat$new_stage)
# qkm2(combdat$death/365, combdat$status, 5, combdat$new_stage2)

m<-coxph(Surv(death, status)~group2+new_stage+Age, data=combdat)
vars <- c("Genotype", "Stage","Age")
HRs <- c(0.46, 1.66, 1.01)
lowerCI <- c(0.28,1.25,1.001)
upperCI <- c(0.75,2.21,1.029)
# forestPlot(vars, HRs, lowerCI, upperCI)
  
d <- data.frame(x = vars, y = HRs, ylo = lowerCI, yhi = upperCI)
xlb <- c(0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 10)

# pdf("Results_Figures_and_P_Values/Supplementary_Fig3_OS_Ovary.pdf",width=7,height=5)
# p <- ggplot(d, aes(x = x, y = y, ymin = ylo, ymax = yhi)) + geom_pointrange() + coord_flip()
# p <- p + scale_y_log10(labels = xlb, breaks = xlb) 
# p <- p + geom_hline(yintercept = 1, lty = 2)
# p <- p + xlab("") + ylab("Hazard Ratio")
# p <- p + myggplotdefaults() 
# + theme(axis.line.y = element_blank())
# dev.off()
