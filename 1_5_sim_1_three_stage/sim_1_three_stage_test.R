
rm(list = ls())

library(weights)
library(doParallel)

setwd("~/Dropbox/Research/AbbVie/Ratio_master_protocol/code/v2/1_5_sim_1_three_stage/")

n.itt = 10^6
n.ind = 8
results.mat = matrix(NA, nrow = n.ind, ncol = 8+13)
n.cluster = 8

m2.func = function(w.vec.in){
  w1.in = w.vec.in[1]
  w2.in = w.vec.in[2]
  w3.in = 1-w.vec.in[1]-w.vec.in[2]
  
  m2.est.in = w1.in*(mean(c(x_t1)) - mean(c(x_p1))) + 
    w2.in*(mean(c(x_t2)) - mean(c(x_p2)))+ 
    w3.in*(mean(c(x_t3)) - mean(c(x_p3)))
  
  m2.p.val.in = 1-pnorm(m2.est.in/sqrt(w1.in^2*var(x_t1)/n_t1+
                                         w1.in^2*var(x_p1)/n_p1+
                                         w2.in^2*var(x_t2)/n_t2+
                                         w2.in^2*var(x_p2)/n_p2+
                                         w3.in^2*var(x_t3)/n_t3+
                                         w3.in^2*var(x_p3)/n_p3
                                         ))
  
  return(c(m2.est.in, m2.p.val.in))
}

for (ind in c(1:8)){

print(ind)
if (ind==1){mu_p1 = 0; mu_p2 = 0.7; mu_p3 = 1; 
sd_p = 1; sd_t = 4; mu_diff = 0} 
if (ind==2){mu_p1 = 0; mu_p2 = 0.7; mu_p3 = 1; 
sd_p = 4; sd_t = 1; mu_diff = 0} 
if (ind==3){mu_p1 = 0; mu_p2 = 1.5; mu_p3 = 0.5; 
sd_p = 1; sd_t = 4; mu_diff = 0} 
if (ind==4){mu_p1 = 0; mu_p2 = 1.5; mu_p3 = 0.5; 
sd_p = 4; sd_t = 1; mu_diff = 0}    

if (ind==5){mu_p1 = 0; mu_p2 = 0.7; mu_p3 = 1; 
sd_p = 1; sd_t = 4; mu_diff = 0.5} 
if (ind==6){mu_p1 = 0; mu_p2 = 0.7; mu_p3 = 1; 
sd_p = 4; sd_t = 1; mu_diff = 0.5} 
if (ind==7){mu_p1 = 0; mu_p2 = 1.5; mu_p3 = 0.5; 
sd_p = 1; sd_t = 4; mu_diff = 0.5} 
if (ind==8){mu_p1 = 0; mu_p2 = 1.5; mu_p3 = 0.5; 
sd_p = 4; sd_t = 1; mu_diff = 0.5}   
  
alpha = 0.05

## total sample size per stage
ratio_1 = 0.8
n_p1 = 200
n_t1 = n_p1*ratio_1

ratio_2 = 0.5
n_p2 = 100
n_t2 = n_p2*ratio_2

ratio_3 = 0.7
n_p3 = 300
n_t3 = n_p3*ratio_3

cl = makeCluster(n.cluster)
registerDoParallel(cl)
out.vec = foreach(itt=1:n.itt) %dopar% {  

# for (itt in 1:n.itt){
  set.seed(itt*ind+n.itt*n.ind)
  x_p1 = rnorm(n_p1, mu_p1, sd_p)
  x_t1 = rnorm(n_t1, mu_p1+mu_diff, sd_t)
  
  x_p2 = rnorm(n_p2, mu_p2, sd_p)
  x_t2 = rnorm(n_t2, mu_p2+mu_diff, sd_t)
  
  x_p3 = rnorm(n_p3, mu_p3, sd_p)
  x_t3 = rnorm(n_t3, mu_p3+mu_diff, sd_t)
   
  ## regression
  time.1 = Sys.time()
  data.reg = data.frame("y" = c(x_p1, x_p2, x_p3, x_t1, x_t2, x_t3),
                        "grp" = c(rep(0, n_p1+n_p2+n_p3), rep(1, n_t1+n_t2+n_t3)),
                        "drift_1" = c(rep(0, n_p1), rep(1, n_p2), rep(0, n_p3),
                                    rep(0, n_t1), rep(1, n_t2), rep(0, n_t3)),
                        "drift_2" = c(rep(0, n_p1), rep(0, n_p2), rep(1, n_p3),
                                      rep(0, n_t1), rep(0, n_t2), rep(1, n_t3)),
                        "true_w" = c(rep(1/sd_p^2, n_p1), 
                                     rep(1/sd_p^2, n_p2),
                                     rep(1/sd_p^2, n_p3),
                                     rep(1/sd_t^2, n_t1), 
                                     rep(1/sd_t^2, n_t2),
                                     rep(1/sd_t^2, n_t3))
  )
 
  lm.fit.2 = lm(y~grp+drift_1 + drift_2, data = data.reg, weights = true_w)
  
  #######################################################################
  ## method 1: direct method
  # m1.p.val = t.test(c(x_p1, x_p2), c(x_t1, x_t2),
  #                   var.equal = FALSE)
  
  m1.est = mean(c(x_t1, x_t2, x_t3)) - 
    mean(c(x_p1, x_p2, x_p3))
  m1.se = sqrt((var(c(x_t1, x_t2, x_t3))/(n_t1+n_t2+n_t3))+
                 (var(c(x_p1, x_p2, x_p3))/(n_p1+n_p2+n_p3)))
  m1.p.val = 1-pnorm(m1.est / m1.se)
 
  ########################################################################
  #### method 2: stratified
  ## weight with known variance, optimal unknown weight
  w.true.1 = 1/(sd_t^2/n_t1+sd_p^2/n_p1)
  w.true.2 = 1/(sd_t^2/n_t2+sd_p^2/n_p2)
  w.true.3 = 1/(sd_t^2/n_t3+sd_p^2/n_p3)
  
  w.true.vec = c(w.true.1/(w.true.1+w.true.2+w.true.3),
                 w.true.2/(w.true.1+w.true.2+w.true.3),
                 w.true.3/(w.true.1+w.true.2+w.true.3))
  w1.e1.fit = m2.func(w.true.vec)
  
  ## weight that is estimated by data
  w.est.1 = 1/(var(x_t1)/n_t1+var(x_p1)/n_p1)
  w.est.2 = 1/(var(x_t2)/n_t2+var(x_p2)/n_p2)
  w.est.3 = 1/(var(x_t3)/n_t3+var(x_p3)/n_p3)
  
  w.est.vec = c(w.est.1/(w.est.1+w.est.2+w.est.3),
                 w.est.2/(w.est.1+w.est.2+w.est.3),
                 w.est.3/(w.est.1+w.est.2+w.est.3))
  w1.e2.fit = m2.func(w.est.vec)
  
  ## weight that is identical to IPTW
  w.IPTW.vec = c((ratio_1+1)/(ratio_1 + ratio_2 + ratio_3 + 3),
                 (ratio_2+1)/(ratio_1 + ratio_2 + ratio_3 + 3),
                 (ratio_3+1)/(ratio_1 + ratio_2 + ratio_3 + 3)
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

out.mat = matrix(unlist(out.vec), nrow = n.itt, ncol=10, byrow = TRUE)

results.mat[ind, ] = c(
  mu_p1, mu_p2, mu_p3, sd_p, sd_t, mu_diff, 
  apply(out.mat, 2, mean),
  apply(out.mat, 2, var)[1:5]+(apply(out.mat, 2, mean)[1:5])^2
  # apply(out.mat, 2, var)[1:5]
                )

}

colnames(results.mat) = c(
  "mu_p1", "mu_p2", "mu_p3", "sd_p", "sd_t", "mu_diff", 
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

write.csv(results.mat, "sim1_three_stage.csv")

##########################################################
## latex table
library(xtable)
results.mat = data.frame(results.mat)

data.out = data.frame(
  "scen" = rep(c("S12", "S13", "S14", "S15"), 2), 
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
  "power_m3_IPTW" = paste0(sprintf("%.2f", results.mat$power_m2_e3*100), "%"),
  "power_m2_hat" = paste0(sprintf("%.2f", results.mat$power_m2_e2*100), "%"),
  "power_m2_opt" = paste0(sprintf("%.2f", results.mat$power_m2_e1*100), "%")
)

print(xtable(data.out), include.rownames = FALSE)

stop("a")

results.mat = data.frame(results.mat)

data.out = data.frame(
  "scen" = rep(c("S12", "S13", "S14", "S15"), 2), 
  "mu_t" = results.mat$mu_diff,
  "m2_opt" = paste0(sprintf("%.6f", results.mat$bias_m2_e1*100), " (", 
                    sprintf("%.6f", results.mat$MSE_m2_e1*100), ")"),
  "WLS" = paste0(sprintf("%.6f", results.mat$bias_WLS*100), " (", 
                    sprintf("%.6f", results.mat$MSE_WLS*100), ")")
)

print(xtable(data.out), include.rownames = FALSE)




































