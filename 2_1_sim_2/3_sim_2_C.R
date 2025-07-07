
rm(list = ls())

library(weights)
library(doParallel)

setwd("~/Dropbox/Research/AbbVie/Ratio_master_protocol/code/v2/2_1_sim_2/")

n.itt = 10^6
n.ind = 6
results.mat = matrix(NA, nrow = n.ind, ncol = 9+13)
n.cluster = 8

m2.func = function(w.vec.in){
  w1.in = w.vec.in[1]
  w2.in = w.vec.in[2]
  
  m2.est.in = w1.in*(mean(c(x_t2)) - mean(c(x_p2))) + 
    w2.in*(mean(c(x_t3)) - mean(c(x_p3)))
  
  m2.p.val.in = 1-pnorm(m2.est.in/sqrt(w1.in^2*var(x_t2)/n_t2+
                                         w1.in^2*var(x_p2)/n_p2+
                                         w2.in^2*var(x_t3)/n_t3+
                                         w2.in^2*var(x_p3)/n_p3
                                         ))
  
  return(c(m2.est.in, m2.p.val.in))
}

for (ind in c(1:6)){

  sd_t1 = 1.4; sd_t2 = 2.7; sd_t3 = 2; sd_t4 = 3.3  
  sd_p1 = 2; sd_p2 = 1.2; sd_p3 = 3.5; sd_p4 = 2.9
  
  print(ind)
  if (ind==1){mu_p1 = 11.08; mu_p2 = 10.3; mu_p3 = 10.45; mu_p4 = 11.6;
  mu_diff = 0}  
  if (ind==2){mu_p1 = 11.08; mu_p2 = 12.3; mu_p3 = 12.45; mu_p4 = 11.6;
  mu_diff = 0}    
  if (ind==3){mu_p1 = 11.6; mu_p2 = 12.45; mu_p3 = 10.3; mu_p4 = 11.08;
  mu_diff = 0}  
  
  if (ind==4){mu_p1 = 11.08; mu_p2 = 10.3; mu_p3 = 10.45; mu_p4 = 11.6;
  mu_diff = 0.45}  
  if (ind==5){mu_p1 = 11.08; mu_p2 = 12.3; mu_p3 = 12.45; mu_p4 = 11.6;
  mu_diff = 0.45}    
  if (ind==6){mu_p1 = 11.6; mu_p2 = 12.45; mu_p3 = 10.3; mu_p4 = 11.08;
  mu_diff = 0.45}  

  
alpha = 0.05

## total sample size per stage
n_p1 = round(200/(1+sqrt(2)))
n_t1 = round(200/(1+sqrt(2))/sqrt(2))

n_p2 = round(600/(1+sqrt(3)))
n_t2 = round(600/(1+sqrt(3))/sqrt(3))

n_p3 = round(300/(1+sqrt(3)))
n_t3 = round(300/(1+sqrt(3))/sqrt(3))

n_p4 = round(200/(1+sqrt(2)))
n_t4 = round(200/(1+sqrt(2))/sqrt(2))

ratio_1 = n_t1/n_p1
ratio_2 = n_t2/n_p2
ratio_3 = n_t3/n_p3
ratio_4 = n_t4/n_p4

cl = makeCluster(n.cluster)
registerDoParallel(cl)
out.vec = foreach(itt=1:n.itt) %dopar% {  

# for (itt in 1:n.itt){
  set.seed(itt*ind+n.itt*n.ind)
  x_p1 = rnorm(n_p1, mu_p1, sd_p1)
  x_t1 = rnorm(n_t1, mu_p1+mu_diff, sd_t1)
  
  x_p2 = rnorm(n_p2, mu_p2, sd_p2)
  x_t2 = rnorm(n_t2, mu_p2+mu_diff, sd_t2)
  
  x_p3 = rnorm(n_p3, mu_p3, sd_p3)
  x_t3 = rnorm(n_t3, mu_p3+mu_diff, sd_t3)

  x_p4 = rnorm(n_p4, mu_p4, sd_p4)
  x_t4 = rnorm(n_t4, mu_p4+mu_diff, sd_t4)
  
  ## regression
  time.1 = Sys.time()
  data.reg = data.frame("y" = c(x_p2, x_p3, x_t2, x_t3),
                        "grp" = c(rep(0, n_p2+n_p3), 
                                  rep(1, n_t2+n_t3)),
                        "drift_1" = c(rep(0, n_p2), rep(1, n_p3),
                                      rep(0, n_t2), rep(1, n_t3)),
                        "true_w" = c(rep(1/sd_p2^2, n_p2), 
                                     rep(1/sd_p3^2, n_p3),
                                     rep(1/sd_t2^2, n_t2), 
                                     rep(1/sd_t3^2, n_t3))
  )
  
  lm.fit.2 = lm(y~grp+drift_1, data = data.reg, weights = true_w)
  
  #######################################################################
  ## method 1: direct method
  # m1.p.val = t.test(c(x_p1, x_p2), c(x_t1, x_t2),
  #                   var.equal = FALSE)
  
  m1.est = mean(c(x_t2, x_t3)) - 
    mean(c(x_p2, x_p3))
  m1.se = sqrt((var(c(x_t2, x_t3))/(n_t2+n_t3))+
                 (var(c(x_p2, x_p3))/(n_p2+n_p3)))
  m1.p.val = 1-pnorm(m1.est / m1.se)
 
  ########################################################################
  #### method 2: stratified
  ## weight with known variance, optimal unknown weight
  w.true.1 = 1/(sd_t2^2/n_t2+sd_p2^2/n_p2)
  w.true.2 = 1/(sd_t3^2/n_t3+sd_p3^2/n_p3)
  
  w.true.vec = c(w.true.1/(w.true.1+w.true.2),
                 w.true.2/(w.true.1+w.true.2))
  w1.e1.fit = m2.func(w.true.vec)
  
  ## weight that is estimated by data
  w.est.1 = 1/(var(x_t2)/n_t2+var(x_p2)/n_p2)
  w.est.2 = 1/(var(x_t3)/n_t3+var(x_p3)/n_p3)
  
  w.est.vec = c(w.est.1/(w.est.1+w.est.2),
                 w.est.2/(w.est.1+w.est.2))
  w1.e2.fit = m2.func(w.est.vec)
   
  ## weight that is identical to IPTW
  w.IPTW.vec = c((n_p2+n_t2)/(n_p2+n_t2+n_p3+n_t3),
                 (n_p3+n_t3)/(n_p2+n_t2+n_p3+n_t3)
  )
  
  w1.e3.fit = m2.func(w.IPTW.vec)
 
  return( c(
    m1.est - mu_diff,
    w1.e1.fit[1] - mu_diff,
    w1.e2.fit[1] - mu_diff,
    w1.e3.fit[1] - mu_diff,
    summary(lm.fit.2)$coefficients[2,1] - mu_diff,
    m1.p.val<=alpha, 
    w1.e1.fit[2]<=alpha,
    w1.e2.fit[2]<=alpha,
    w1.e3.fit[2]<=alpha,
    (1-pt(summary(lm.fit.2)$coefficients[2,3], lm.fit.2$df.residual)) <= alpha
                     ))
}

out.mat = matrix(unlist(out.vec),nrow = n.itt, ncol=10, byrow = TRUE)

results.mat[ind, ] = c(
  mu_p1, mu_p2, mu_p3, mu_p4, sd_p1, sd_t1, mu_diff, 
  apply(out.mat, 2, mean),
  apply(out.mat, 2, var)[1:5]+(apply(out.mat, 2, mean)[1:5])^2
                )

}

colnames(results.mat) = c(
  "mu_p1", "mu_p2", "mu_p3", "mu_p4", "sd_p", "sd_t", "mu_diff", 
  "bias_m1",
                          "bias_m2_e1",
                          "bias_m2_e2",
  "bias_m2_e3",
  "bias_WLS",
                          "power_m1",
                          "power_m2_e1",
                          "power_m2_e2",
  "power_m2_e3",
  "power_WLS",
                          "MSE_m1",
                          "MSE_m2_e1",
                          "MSE_m2_e2",
  "MSE_m2_e3",
  "MSE_WLS"
                          )

print(results.mat)

##########################################################
## latex table
library(xtable)
results.mat = data.frame(results.mat)

data.out = data.frame(
  "scen" = rep(c("S1", "S2", "S3"), 2), 
  "mu_t" = results.mat$mu_diff,
  "m1" = paste0(sprintf("%.2f", results.mat$bias_m1*100), " (", 
                sprintf("%.2f", results.mat$MSE_m1*100), ")"),
  "m2_IPTW" = paste0(sprintf("%.2f", results.mat$bias_m2_e3*100), " (", 
                     sprintf("%.2f", results.mat$MSE_m2_e3*100), ")"),
  "m2_hat" = paste0(sprintf("%.2f", results.mat$bias_m2_e2*100), " (", 
                    sprintf("%.2f", results.mat$MSE_m2_e2*100), ")"),
  "m2_opt" = paste0(sprintf("%.2f", results.mat$bias_m2_e1*100), " (", 
                sprintf("%.2f", results.mat$MSE_m2_e1*100), ")"),
  "power_m1" = paste0(sprintf("%.2f", results.mat$power_m1*100), "%"),
  "power_m2_IPTW" = paste0(sprintf("%.2f", results.mat$power_m2_e3*100), "%"),
  "power_m2_hat" = paste0(sprintf("%.2f", results.mat$power_m2_e2*100), "%"),
  "power_m2_opt" = paste0(sprintf("%.2f", results.mat$power_m2_e1*100), "%")
)

print(xtable(data.out), include.rownames = FALSE)








































