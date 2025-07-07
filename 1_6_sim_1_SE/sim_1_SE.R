
rm(list = ls())

library(weights)
library(doParallel)

setwd("~/Dropbox/Research/AbbVie/Ratio_master_protocol/code/v2/1_6_sim_1_SE/")

n.itt = 10^6
n.ind = 10
results.mat = matrix(NA, nrow = n.ind, ncol = 15+7+9+2)
n.cluster = 8

m2.func = function(w1.in){
  w2.in = 1-w1.in
  m2.est.in = w1.in*(mean(c(x_t1)) - mean(c(x_p1))) + 
    w2.in*(mean(c(x_t2)) - mean(c(x_p2)))
  
  m2.p.val.in = 1-pnorm(m2.est.in/sqrt(w1.in^2*var(x_t1)/n_p1/ratio_1+
                                         w1.in^2*var(x_p1)/n_p1+
                                         w2.in^2*var(x_t2)/n_p2/ratio_2+
                                         w2.in^2*var(x_p2)/n_p2))
  
  return(c(m2.est.in, m2.p.val.in))
}

for (ind in c(1:10)){

print(ind)
mu_t = 0; mu_add = 0; ratio_2 = 1; sd_p1 = sd_p2 = 2; sd_t1 = sd_t2 = 2  
   
mu_p = 0
ratio_1 = 1
n_p1 = 120
n_p2 = 120
alpha = 0.05
n_t1 = n_p1*ratio_1
n_t2 = n_p2*ratio_2
mu_diff = mu_t - mu_p
# out.mat = matrix(NA, n.itt, 10)

cl = makeCluster(n.cluster)
registerDoParallel(cl)
out.vec = foreach(itt=1:n.itt) %dopar% {  

# for (itt in 1:n.itt){
  set.seed(itt+n.itt*ind*10)
  x_p1 = rnorm(n_p1, mu_p, sd_p1)
  x_p2 = rnorm(n_p2, mu_p+mu_add, sd_p2)
  
  x_t1 = rnorm(n_t1, mu_p + mu_t, sd_t1)
  x_t2 = rnorm(n_t2, mu_p + mu_t+mu_add, sd_t2)
  
  ## regression
  time.1 = Sys.time()
  data.reg = data.frame("y" = c(x_p1, x_p2, x_t1, x_t2),
                        "grp" = c(rep(0, n_p1+n_p2), rep(1, n_t1+n_t2)),
                        "drift" = c(rep(0, n_p1), rep(1, n_p2),
                                    rep(0, n_t1), rep(1, n_t2)),
                        "true_w" = c(rep(1/sd_p1^2, n_p1), 
                          rep(1/sd_p2^2, n_p2),
                          rep(1/sd_t1^2, n_t1), 
                          rep(1/sd_t2^2, n_t2))
  )
  int.lm.fit = lm(y~grp+drift, data = data.reg)
  
  data.reg$w = c(rep(1/mean(int.lm.fit$residuals[data.reg$grp==0&data.reg$drift==0]^2), n_p1), 
                 rep(1/mean(int.lm.fit$residuals[data.reg$grp==0&data.reg$drift==1]^2), n_p2),
                 rep(1/mean(int.lm.fit$residuals[data.reg$grp==1&data.reg$drift==0]^2), n_t1), 
                 rep(1/mean(int.lm.fit$residuals[data.reg$grp==1&data.reg$drift==1]^2), n_t2))
  
  lm.fit.3 = lm(y~grp+drift, data = data.reg, weights = w)
  
  time.1.diff = (Sys.time()-time.1)
  
  lm.fit.1 = int.lm.fit
  lm.fit.2 = lm(y~grp+drift, data = data.reg, weights = true_w)
  
  #######################################################################
  ## method 1: direct method
  # m1.p.val = t.test(c(x_p1, x_p2), c(x_t1, x_t2),
  #                   var.equal = FALSE)
  
  m1.est = mean(c(x_t1, x_t2)) - mean(c(x_p1, x_p2))
  m1.se = sqrt((var(c(x_t1, x_t2))/(n_p1*ratio_1 + n_p2*ratio_2))+
                 (var(c(x_p1, x_p2))/(n_p1 + n_p2)))
  m1.p.val = 1-pnorm(m1.est / m1.se)
 
  ########################################################################
  #### method 2: stratified
  # s1.p.val = t.test(c(x_p1), c(x_t1),
  #                   var.equal = FALSE)$p.value
  # s2.p.val = t.test(c(x_p2), c(x_t2),
  #                   var.equal = FALSE)$p.value
  # m2.p.val = 1 - pnorm(qnorm(1-s1.p.val)*w1+qnorm(1-s2.p.val)*w2)
  # m2.p.val = 0
  
  ## weight with known variance, optimal unknown weight
  w1.e1 = ((sd_t2^2/ratio_2+sd_p2^2)/
             (sd_t1^2/ratio_1+sd_t2^2/ratio_2+sd_p1^2+sd_p2^2))
  w1.e1.fit = m2.func(w1.e1)
  
  ## weight that is estimated by data
  time.2 = Sys.time()
  w1.e2 = ((1/ratio_2*var(x_t2)+var(x_p2))/
          (1/ratio_1*var(x_t1)+1/ratio_2*var(x_t2)+var(x_p1)+var(x_p2)))
  w1.e2.fit = m2.func(w1.e2)
  
  time.2.diff = (Sys.time()-time.2)
   
  ## weight that is identical to IPTW
  w1.e3 = ((ratio_1+1)/
          (ratio_1+ratio_2+2))
  w1.e3.fit = m2.func(w1.e3)
  
  ## weight that is estimated by equal-variance assumption
  w1.e4 = ((1/ratio_2+1)/
             (1/ratio_1+1/ratio_2+2))
  w1.e4.fit = m2.func(w1.e4)
  
  ## method 3: IPTW
  # weight.t.vec = (c(rep(ratio_1+1, n_p),
  #                  rep(ratio_2+1, n_p),
  #                  rep((ratio_1+1)/ratio_1, (ratio_1)*n_p),
  #                  rep((ratio_2+1)/ratio_2, (ratio_2)*n_p)))
  # weight.t.vec = weight.t.vec
  # 
  # # weight.t.vec = (c(rep((1+1)/1, n_p),
  # #                   rep((1+1)/1, n_p),
  # #                   rep((ratio_1+ratio_2)/ratio_1, (ratio_1)*n_p),
  # #                   rep((ratio_2+ratio_2)/ratio_2, (ratio_2)*n_p)))
  # 
  # # weight.t.vec = weight.t.vec/sum(weight.t.vec)
  # 
  # data.in = data.frame("x" = c(c(x_p1, x_p2), c(x_t1, x_t2)),
  #                      "s" = c(rep(1, n_p), rep(2, n_p), 
  #                              rep(1, n_p*ratio_1), rep(2, n_p*ratio_2)),
  #                      "grp" = c(rep(0, (1+1)*n_p),
  #                                rep(1, (ratio_2+ratio_1)*n_p)),
  #                      "w" = weight.t.vec)
  # # m3.p.val = coef(summary(lm(x~grp, data = data.in, weights = w)))[2, 4]
  # m3.p.val = 0
  # 
  # m3.est = sum((data.in$x*data.in$w)[data.in$grp==1])/n_p/(2+ratio_1+ratio_2)-
  #   sum((data.in$x*data.in$w)[data.in$grp==0])/n_p/(2+ratio_1+ratio_2)

  return( c(
    m1.est - mu_diff,
    w1.e1.fit[1] - mu_diff,
    w1.e2.fit[1] - mu_diff,
    w1.e3.fit[1] - mu_diff,
    w1.e4.fit[1] - mu_diff,
    summary(lm.fit.1)$coefficients[2,1] - mu_diff,
    summary(lm.fit.2)$coefficients[2,1] - mu_diff,
    summary(lm.fit.3)$coefficients[2,1] - mu_diff,
    m1.p.val<=alpha, 
    w1.e1.fit[2]<=alpha,
    w1.e2.fit[2]<=(alpha),
    w1.e3.fit[2]<=alpha,
    w1.e4.fit[2]<=alpha,
    (1-pt(summary(lm.fit.1)$coefficients[2,3], lm.fit.1$df.residual)) <= alpha,
    (1-pt(summary(lm.fit.2)$coefficients[2,3], lm.fit.2$df.residual)) <= alpha,
    (1-pt(summary(lm.fit.3)$coefficients[2,3], lm.fit.3$df.residual)) <= alpha,
    time.2.diff, 
    time.1.diff
                     ))
}

out.mat = matrix(unlist(out.vec),nrow = n.itt, ncol=18, byrow = TRUE)

results.mat[ind, ] = c(
  mu_t, mu_add, ratio_2, sd_p1, sd_p2, sd_t1, sd_t2, 
  apply(out.mat, 2, mean),
  apply(out.mat, 2, var)[1:8]+(apply(out.mat, 2, mean)[1:8])^2
                )

}

colnames(results.mat) = c(
  "mu_t", "mu_add", "ratio_2", "sd_p1", "sd_p2", "sd_t1", "sd_t2",
  "bias_m1",
                          "bias_m2_e1",
                          "bias_m2_e2",
                          "bias_m2_e3",
                          "bias_m2_e4",
                          "bias_lm",
  "bias_wlm_true",
  "bias_wlm",
                          "power_m1",
                          "power_m2_e1",
                          "power_m2_e2",
                          "power_m2_e3",
                          "power_m2_e4",
                          "power_lm",
  "power_wlm_true",
  "power_wlm",
  "time_m2_e2",
  "time_wlm",
                          "MSE_m1",
                          "MSE_m2_e1",
                          "MSE_m2_e2",
                          "MSE_m2_e3",
                          "MSE_m2_e4",
                          "MSE_lm",
  "MSE_wlm_true",
  "MSE_wlm"
                          )

print(results.mat)

write.csv(results.mat, "sim1_results.csv")

results.mat = data.frame(results.mat)
