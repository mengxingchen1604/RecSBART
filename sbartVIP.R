library(MASS)
library(msm)
library(glmnet)
library(SoftBart)
library(diversitree)
library(extraDistr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(patchwork)
library(writexl)
#competing methods
library(recforest)
source("Variable Inclusion proportion for RecSBART.R")
#dataset
library(frailtypack)
data("readmission")

#set up dataset
readmission$gender<-as.numeric(readmission$sex=="Male")
readmission$group<-as.numeric(readmission$chemo=="Treated")
readmission$dukes1<-as.numeric(readmission$dukes)

readmission$charlson_base <- as.numeric(ave(readmission$charlson, 
                                            readmission$id, 
                                            FUN = function(x) rep(x[1], length(x))))
#1-2 = 1.5
readmission$charlson_base<-sapply(readmission$charlson_base,function(x)if(x==2){x=1.5} else{x})
#A-B = 0.5
readmission$dukes1<-sapply(readmission$dukes1,function(x)if(x==1){x=0.5} else{x})

#convert dataset

#X is total for bayesian methods
X <- matrix(NA, nrow = 403, ncol = 4)
for (i in 1:403) {
  m <- subset(readmission, id == i)
  X[i, ] <- as.numeric(m[1, 12:15])
}

ID<-sort(unique(readmission$id))
readmission$terminal<-ave(readmission$t.stop,readmission$id,FUN = max)
t<- subset(readmission, t.stop==terminal)
terminal<-as.numeric(t$terminal)


recurrent_event <- lapply(split(readmission, readmission$id), function(x) {
  as.numeric(subset(x, event == 1)$t.stop)
})

recurrent_event<-lapply(recurrent_event,function(i) i/max(terminal))
terminal<- terminal/max(terminal)

readmission_softbart<-fit_RecSBART_VIP(X_train = X,
                                              X_test = X,
                                              recurrent_train = recurrent_event,
                                              recurrent_test = recurrent_event,
                                              terminal_train = terminal,
                                              terminal_test = terminal,
                                              num_burn = 2500,
                                              num_thin = 1,
                                              num_save = 2500)
#VIP for four covariates
readmission_softbart<-readmission_softbart[2:5,]
sum_dom<-colSums(readmission_softbart)
#gender

vip_gender<-mean(readmission_softbart[1,]/sum_dom)

#treatment
sum_dom_trt<-colSums(readmission_softbart)
vip_trt<-mean(readmission_softbart[2,]/sum_dom)

#dukes
sum_dom_duke<-colSums(readmission_softbart)
vip_duke<-mean(readmission_softbart[3,]/sum_dom)

#charlson
sum_dom_charlson<-colSums(readmission_softbart)
vip_charlson <- mean(readmission_softbart[4,]/sum_dom)