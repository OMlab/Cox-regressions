setwd("H:/MDC/sample_sheets/")	
library(openxlsx)
library(survival)
library(tidyverse)
library(ggplot2)

## This is a script with cox regressions using the survival package
## Known bug: test for proportionallity handles missing values in a weird way
## If the last metabolite has no missing vaules, the loops seems to be working
## This script has no normalization algorithm, instead, the values in the input are pre-normalized.
## The script is created with two models for different adjustments

## next thing to do is to add the covariate input to the top


# Import Data
 Import <- read.csv("MDC2018_metabolites.csv") 	# Specify data file


Data <- Import %>%
filter(vital_st != 3 & pr_af == 0)         #Filter metabolites used in analysis from excel document

start <- (which(names(Data)=="Start"))+1   #define start of metabolites  
  
filenamestart <- "coxreg_"          #defines the start of the excel file
path <- "H:/AF/preliminary/newadjust/test/"  #Enter path to export data

DV <- "inc_af"			                       #Enter dependent variable
follow_up <- "fuaf"                        #enter time varible

create_excel = FALSE                       
create_excel = TRUE                       # Run this line if you want an excel file to be created
create_ggplot = 2 #                    # 1 = classic ggplot, 2 = metabolites on x-axis, 0 = no ggplot

## Remakes the input from DV and follo_up into numeric
Endpoint <- as.numeric(t(Data[paste(DV)]))		
Time <- as.numeric(t(Data[paste(follow_up)]))

## for specific heart failure analysis, remakes the variable to 1/0 instead of 1/2/0
# Data$inc_hf[Data$inc_hf == 2 ] <- 1
# Data$inc_hf[Data$pr_hf == 2 ] <- 1


coxreg=NULL


a <- seq(from=start,to=ncol(Data),by=1)
for (i in 1:length(a)){
  coxreg[[i]] <- coxph(Surv(Time, Endpoint) ~ Data[,a[i]] 
                       
                        +age+female
                      # + SBP+current_smoker+ALKO_B+pr_dm+AHT_B+pr_hf+pr_mc
                      # +LNNtBNP
                       , data = Data
  )
}


coxreg_adjusted <- NULL
for (i in 1:length(a)){
coxreg_adjusted[[i]] <- coxph(Surv(Time, Endpoint) ~ Data[,a[i]]
                     
                     +age+female
                      + SBP+current_smoker+ALKO_B+pr_dm+AHT_B+pr_hf+pr_mc
                     # +LNNtBNP
                     , data = Data
)
}
##### creation of output functions#####
e <- seq(from=1,to=length(a),by=1) # length of output
HR <- function(n){summary(coxreg[[n]])$coefficients[,2][1]}
p <- function(n){summary(coxreg[[n]])$coefficients[,5][1]}
lowint <- function(n){exp(confint(coxreg[[n]])[,1][1])}
highint <- function(n){exp(confint(coxreg[[n]])[,2][1])}
p_pHz_mod <- function(n){cox.zph(coxreg[[n]])$table[,"p"][length(cox.zph(coxreg[[n]])$table[,"p"])]}
p_pHz_met <- function(n){cox.zph(coxreg[[n]])$table[,"p"][1]}

HR_adjusted <- function(n){summary(coxreg_adjusted[[n]])$coefficients[,2][1]}
p_adjusted <- function(n){summary(coxreg_adjusted[[n]])$coefficients[,5][1]}
lowint_adjusted <- function(n){exp(confint(coxreg_adjusted[[n]])[,1][1])}
highint_adjusted <- function(n){exp(confint(coxreg_adjusted[[n]])[,2][1])}
p_pHz_mod_adjusted <- function(n){cox.zph(coxreg_adjusted[[n]])$table[,"p"][length(cox.zph(coxreg_adjusted[[n]])$table[,"p"])]}
p_pHz_met_adjusted <- function(n){cox.zph(coxreg_adjusted[[n]])$table[,"p"][1]}
#####ouput for the first model #####
output <- data.frame(matrix(nrow=length(e),ncol=1))
output["metabolites"] <- c(names(Data[start:ncol(Data)]))
output["p"] <- sapply(e,"p")
output["p_fdr"] <- p.adjust(output["p"][,1],method="fdr",n=nrow(output["p"]))
output["p_unique"] <- sapply(e,"p") # this variable is for sorting both data-sets according to the p-value of model 1
output["p_selection"] <- p.adjust(output["p"][,1],method="fdr",n=nrow(output["p"])) # this variable is for selection to the ggplot
output["HR"] <- sapply(e,"HR")
output["lowint"] <- sapply(e,"lowint")
output["highint"] <- sapply(e,"highint")
output["p_pHz_met"] <- vapply(e,"p_pHz_met", numeric(1))
output["p_pHz_mod"] <- sapply(e,"p_pHz_mod")
output["model"] <- as.factor(1) #using as factor makes the ggplot behave

#####output for the second model#####
output2 <- data.frame(matrix(nrow=length(e),ncol=1))
output2["metabolites"] <- c(names(Data[start:ncol(Data)]))
output2["p"] <- sapply(e,"p_adjusted")
output2["p_fdr"] <- p.adjust(output2["p"][,1],method="fdr",n=nrow(output2["p"]))
output2["p_unique"] <- sapply(e,"p") # this variable is for sorting both data-sets according to the p-value of model 1
output2["p_selection"] <- p.adjust(output["p"][,1],method="fdr",n=nrow(output["p"])) # this variable is for selection to the ggplot
output2["HR"] <- sapply(e,"HR_adjusted")
output2["lowint"] <- sapply(e,"lowint_adjusted")
output2["highint"] <- sapply(e,"highint_adjusted")
output2["p_pHz_met"] <- vapply(e,"p_pHz_met_adjusted", numeric(1))
output2["p_pHz_mod"] <- sapply(e,"p_pHz_mod_adjusted")
output2["model"] <- as.factor(2)

##### cleaning up the output, orders them by p-value and combines them #####
drops <- c("matrix.nrow...length.e...ncol...1.")
output <- output[ , !(names(output) %in% drops)]
output2 <- output2[ , !(names(output2) %in% drops)]
output <- output[order(output$p,decreasing=FALSE),]
output2 <- output2[order(output2$p_unique,decreasing=FALSE),]
output_both <-rbind(output, output2)

filetype <- c(".xlsx",".pdf",".eps")
##decides the selection for the graph. "p_selection" is the fdr-p-calue of model 1

output_graph <- output_both %>%
  group_by(model) %>%
  filter(p_selection < 0.05) #p selection is the fdr-value of model one in both dataframes, therfore, it selects all fdr-significant in model 1 to be analyzed in model 2 aswell
 

##### Create ggplots #####
p1 <- ggplot(output_graph, aes(x=reorder(metabolites,-p), y=HR, group=model, color=model, ymin=lowint,ymax=highint)) +
  geom_pointrange(position = position_dodge(width =-0.5))+
  theme(panel.background = element_rect(fill="white"),axis.text.y = element_text(size =7),axis.ticks.y = element_line(size=0.1))+
  geom_point(size=3, shape=21,position = position_dodge(width =-0.5)) +
  coord_flip() +
  geom_hline(yintercept = 1,linetype=2)+
  xlab("Metabolites")+
  ylab("HR")

p2 <- ggplot(output_graph, aes(x=reorder(metabolites,p_unique), y=HR, color=model
                       , ymin=lowint,ymax=highint)) + 
  geom_pointrange(position = position_dodge(width =0.5))+ 
  scale_y_discrete(limits = c(0,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,2))+
  theme(panel.background = element_rect(fill="white"),
        axis.text.y = element_text(size =7),axis.ticks.y = element_line(size=0.1),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank())+
  geom_point(size=1, position = position_dodge(width =0.5)) +
  labs(colour = "Model") +
  geom_hline(yintercept = 1,linetype=3)
p3 <- p2 + geom_text(data=output_graph, 
                     aes( x = metabolites, y= 1.3, label = ifelse(model==1,metabolites,""),angle=45 ), 
                     size = 2, show.legend = FALSE, colour = "black")



##### creates a new data-frame and exports it to excel#####
output$model <- 1
output2$model <- 2
output_excel <- cbind(output, output2)

output_excel <- output_excel[c(-4,-5,-15,-16)] #removes p_selection and p_unique
if(create_excel == TRUE){
  write.xlsx(output_excel,paste(path,filenamestart,DV,filetype[1],sep=""))
  
  paste("Excel file created ",path,filenamestart,DV,filetype[1],sep="")
} else {
  paste("Excel file is not created as per instruction")
}
##### exports ggplot  #####

## here is a selection between two different ggplots depening on the create_ggplot_variable
if(create_ggplot == 1){
  pdf(paste(path,filenamestart,DV,filetype[2],sep=""))
  print(p1) 
   dev.off()
  
  paste("Classic ggplot created ",path,filenamestart,DV,filetype[2],sep="")
} else if(create_ggplot == 2){ 
  pdf(paste(path,filenamestart,DV,filetype[2],sep=""))
  print(p3) 
  dev.off()
  
  paste("Modern ggplot created ",path,filenamestart,DV,filetype[2],sep="")
}else{
  paste("ggplot not created as per instruction")
}

