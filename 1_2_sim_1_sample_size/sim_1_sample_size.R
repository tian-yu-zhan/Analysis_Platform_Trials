
rm(list = ls())

library(weights)
library(doParallel)

setwd("~/Dropbox/Research/AbbVie/Ratio_master_protocol/code/v2/1_2_sim_1_sample_size/")

n.itt = 10^6
n.ind = 7
results.mat = matrix(NA, nrow = n.ind, ncol = 8)
n.cluster = 8


for (ind in c(1:7)){

print(ind)
  
if (ind==1){mu_t = 0.5; mu_add = 0; ratio_2 = 1; sd_p1 = sd_p2 = 2; 
sd_t1 = sd_t2 = 2} 
if (ind==2){mu_t = 0.5; mu_add = 0; ratio_2 = 0.5; sd_p1 = sd_p2 = 2; 
sd_t1 = sd_t2 = 2} 
if (ind==3){mu_t = 0.5; mu_add = 0.3; ratio_2 = 1; sd_p1 = sd_p2 = 2; 
sd_t1 = sd_t2 = 2}   
if (ind==4){mu_t = 0.5; mu_add = 0.3; ratio_2 = 0.5; sd_p1 = sd_p2 = 2; 
sd_t1 = sd_t2 = 2}  
if (ind==5){mu_t = 0.6; mu_add = 0.3; ratio_2 = 0.5; sd_p1 = sd_p2 = 1; 
sd_t1 = sd_t2 = 4}  
if (ind==6){mu_t = 0.7; mu_add = 0.3; ratio_2 = 0.5; sd_p1 = sd_p2 = 4; 
sd_t1 = sd_t2 = 1}  
if (ind==7){mu_t = 0.5; mu_add = 0.3; ratio_2 = 0.5; sd_p1 = 1; sd_p2 = 2; 
sd_t1 = 2; sd_t2 = 3}    
  
mu_p = 0
ratio_1 = 1
n_p = 120
alpha = 0.05
n_t1 = n_p*ratio_1
n_t2 = n_p*ratio_2
mu_diff = mu_t - mu_p
# out.mat = matrix(NA, n.itt, 10)

w1.e1 = ((sd_t2^2/ratio_2+sd_p2^2)/
           (sd_t1^2/ratio_1+sd_t2^2/ratio_2+sd_p1^2+sd_p2^2))

w2.e1 = 1-w1.e1
w.e1 = mu_diff

pow.pred = 1-pnorm(qnorm(1-alpha), mean = 
                     w.e1/sqrt(w1.e1^2*sd_t1^2/n_p/ratio_1+
                                 w1.e1^2*sd_p1^2/n_p+
                                 w2.e1^2*sd_t2^2/n_p/ratio_2+
                                 w2.e1^2*sd_p2^2/n_p))


results.mat[ind, ] = c(
  mu_t, mu_add, ratio_2, sd_p1, sd_p2, sd_t1, sd_t2, 
  pow.pred
                )

}

colnames(results.mat) = c(
  "mu_t", "mu_add", "ratio_2", "sd_p1", "sd_p2", "sd_t1", "sd_t2",
  "pred_power"
                          )

print(results.mat)

write.csv(results.mat, "sim1_pow_results.csv")

##########################################################
## latex table
library(xtable)
results.mat = data.frame(results.mat)

results.pow.mat = data.frame(read.csv("sim1_results.csv")[,-1])

data.out = data.frame(
  "scen" = rep(c("S1", "S2", "S3", "S4", "S5", "S6", "S7")), 
  "mu_t" = results.mat$mu_t,
  "pred_pow" = paste0(sprintf("%.2f", results.mat$pred_power*100), "%"),
  "opt" = paste0(sprintf("%.2f", results.pow.mat$power_m2_e1[8:14]*100), "%"),
  "hat" = paste0(sprintf("%.2f", results.pow.mat$power_m2_e2[8:14]*100), "%")
)

print(xtable(data.out), include.rownames = FALSE)







































