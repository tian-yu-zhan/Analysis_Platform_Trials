
rm(list = ls())

library(weights)
library(doParallel)
library(ggplot2)

setwd("~/Dropbox/Research/AbbVie/Ratio_master_protocol/code/v2/1_3_sim_1_diff/")

n.itt = 5*10^3
n.ind = 1
results.mat = matrix(NA, nrow = n.ind, ncol = 15+5)
n.cluster = 8

m2.func = function(w1.in){
  w2.in = 1-w1.in
  m2.est.in = w1.in*(mean(c(x_t1)) - mean(c(x_p1))) + 
    w2.in*(mean(c(x_t2)) - mean(c(x_p2)))
  m2.se.in = sqrt(w1.in^2*var(x_t1)/n_p/ratio_1+
                    w1.in^2*var(x_p1)/n_p+
                    w2.in^2*var(x_t2)/n_p/ratio_2+
                    w2.in^2*var(x_p2)/n_p)
  
  m2.p.val.in = 1-pnorm(m2.est.in/m2.se.in)
  
  return(c(m2.est.in, m2.p.val.in, m2.est.in - m2.se.in*qnorm(1-alpha)))
}

for (ind in c(1)){

print(ind)

mu_t = 0.5; mu_add = 0.3; ratio_2 = 0.5; sd_p = 1; sd_t = 4
 
mu_p = 0
ratio_1 = 1
n_p = 120
alpha = 0.05
n_t1 = n_p*ratio_1
n_t2 = n_p*ratio_2
mu_diff = mu_t - mu_p
# out.mat = matrix(NA, n.itt, 10)

cl = makeCluster(n.cluster)
registerDoParallel(cl)
out.vec = foreach(itt=1:n.itt) %dopar% {  

# for (itt in 1:n.itt){
  set.seed(itt+n.itt*ind)
  x_p1 = rnorm(n_p, mu_p, sd_p)
  x_p2 = rnorm(n_p, mu_p+mu_add, sd_p)
  
  x_t1 = rnorm(n_t1, mu_t, sd_t)
  x_t2 = rnorm(n_t2, mu_t+mu_add, sd_t)
  
  ########################################################################
  #### method 1: stratified
  ## weight that is estimated by equal-variance assumption
  w1.e4 = ((1/ratio_2+1)/
             (1/ratio_1+1/ratio_2+2))
  w1.e4.fit = m2.func(w1.e4)

  #######################################################################
  #### method 2: combination test approach
  p1 = 1-pnorm((mean(x_t1) - mean(x_p1))/
                 sqrt(var(x_t1)/n_p/ratio_1+var(x_p1)/n_p))
  # print(p1)
  # t.test(x_t1, x_p1,alternative = "greater")
  
  p2 = 1-pnorm((mean(x_t2) - mean(x_p2))/
                 sqrt(var(x_t2)/n_p/ratio_2+var(x_p2)/n_p))
  
  p.com = 1-pnorm(qnorm(1-p1)*sqrt(w1.e4) + qnorm(1-p2)*sqrt(1-w1.e4))
  
  return( c(
    mean(x_t1) - mean(x_p1),
    mean(x_t2) - mean(x_p2),
    w1.e4.fit[2]<=alpha,
    w1.e4.fit[3]>=0,
    p.com<=alpha

                     ))
}

out.mat = matrix(unlist(out.vec),nrow = n.itt, ncol=5, byrow = TRUE)

}

###################################################################
## plot 1
out.mat = data.frame(out.mat)
out.mat$cat = 1
out.mat$cat[out.mat$X4==1&out.mat$X5==1] = 2
out.mat$cat[out.mat$X4==1&out.mat$X5==0] = 3
out.mat$cat[out.mat$X4==0&out.mat$X5==1] = 4

out.mat$cat = factor(out.mat$cat)

png("sim_1_diff_1.png", width = 1800, height = 1800)
ggplot.fit = ggplot(out.mat, aes(x = X1, y = X2)) +
  geom_point(size = 8, shape = 16, aes(alpha = cat, color = cat))+
  # scale_size_manual(values=c(10, 4, 4, 3))+
  scale_alpha_manual(values=c(0.1, 0.1, 1, 1))+
  # scale_shape_manual(values=c(15, 15, 16, 17))+
  scale_color_manual(values=c("black", "black", "#E69F00", "#56B4E9"))+
  # geom_point(aes(x=knots,y=0), color = "blue",
  #            data=data.frame("knots"=knots[c(1:6,12:16)]),
  #            size = 16)+
  # geom_point(aes(x=knots,y=0), color = "blue",
  #            data=data.frame("knots"=knots[7:11]),
  #            size = 16)+
  # scale_color_manual(values=c("blue", "red"))+
  # geom_line(aes(color = Method, linetype=Method),
  #           size = 2)
  # scale_linetype_manual(values=c("solid", "longdash", "solid", "longdash", "solid", "solid")) +
  # # scale_alpha_manual(values=c(1, 1, rep(0.7, 5)))+
  # scale_color_manual(values=
  #                      c("black", "#56B4E9", "#0072B2", "#E69F00", "#D55E00", "#009E73"))+
  scale_x_continuous(breaks = c(-0.7, 0, 0.7, 1.4, 2.1), limits = c(-0.7, 2.1)) +
  # # scale_y_continuous(sec.axis = sec_axis(~., name = "Type I error")) +
  scale_y_continuous(breaks = c(-1.5, -0.5, 0.5, 1.5, 2.6), limits = c(-1.5, 2.6))+
  labs(title = "") +
  ylab ("Estimated treatment difference in stage 2") +
  xlab("Estimated treatment difference in stage 1") +
  # theme_bw()+
  theme(plot.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(2,1,1,1),units="lines"),
        text = element_text(size=50),
        axis.text.x = element_text(colour="black",size=60,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=60,angle=0,hjust=1,vjust=0,face="plain"),
        axis.title.x = element_text(colour="black",size=60,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=60,angle=90,hjust=.5,vjust=.5,face="plain"),
        legend.text = element_text(colour="black", size = 37, face = "plain"),
        legend.title = element_text(colour="black", size = 42, face = "plain"),
        legend.key.size = unit(10,"line"),
        legend.position="none", plot.title = element_text(hjust = 0.5),
        legend.box="vertical",
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(10, "lines"))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE),
         linetype=guide_legend(nrow=1,byrow=TRUE),
         color=guide_legend(nrow=1,byrow=TRUE))

print(ggplot.fit)
dev.off()

##############################################
## plot 2
out.mat = data.frame(out.mat)
out.mat$cat = 1
out.mat$cat[out.mat$X3==1&out.mat$X4==1] = 2
out.mat$cat[out.mat$X3==1&out.mat$X4==0] = 3
out.mat$cat[out.mat$X3==0&out.mat$X4==1] = 4

out.mat$cat = factor(out.mat$cat)

png("sim_1_diff_2.png", width = 1800, height = 1800)
ggplot.fit = ggplot(out.mat, aes(x = X1, y = X2)) +
  geom_point(size = 8, shape = 16, aes(alpha = cat, color = cat))+
  # scale_size_manual(values=c(10, 4, 4, 3))+
  scale_alpha_manual(values=c(0.1, 0.1, 1, 1))+
  # scale_shape_manual(values=c(15, 15, 16, 17))+
  scale_color_manual(values=c("black", "black", "#E69F00", "#56B4E9"))+
  # geom_point(aes(x=knots,y=0), color = "blue",
  #            data=data.frame("knots"=knots[c(1:6,12:16)]),
  #            size = 16)+
  # geom_point(aes(x=knots,y=0), color = "blue",
  #            data=data.frame("knots"=knots[7:11]),
  #            size = 16)+
  # scale_color_manual(values=c("blue", "red"))+
  # geom_line(aes(color = Method, linetype=Method),
  #           size = 2)
  # scale_linetype_manual(values=c("solid", "longdash", "solid", "longdash", "solid", "solid")) +
  # # scale_alpha_manual(values=c(1, 1, rep(0.7, 5)))+
# scale_color_manual(values=
#                      c("black", "#56B4E9", "#0072B2", "#E69F00", "#D55E00", "#009E73"))+
scale_x_continuous(breaks = c(-0.7, 0, 0.7, 1.4, 2.1), limits = c(-0.7, 2.1)) +
  # # scale_y_continuous(sec.axis = sec_axis(~., name = "Type I error")) +
  scale_y_continuous(breaks = c(-1.5, -0.5, 0.5, 1.5, 2.6), limits = c(-1.5, 2.6))+
  labs(title = "") +
  ylab ("Estimated treatment difference in stage 2") +
  xlab("Estimated treatment difference in stage 1") +
  # theme_bw()+
  theme(plot.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(2,1,1,1),units="lines"),
        text = element_text(size=50),
        axis.text.x = element_text(colour="black",size=60,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=60,angle=0,hjust=1,vjust=0,face="plain"),
        axis.title.x = element_text(colour="black",size=60,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=60,angle=90,hjust=.5,vjust=.5,face="plain"),
        legend.text = element_text(colour="black", size = 37, face = "plain"),
        legend.title = element_text(colour="black", size = 42, face = "plain"),
        legend.key.size = unit(10,"line"),
        legend.position="none", plot.title = element_text(hjust = 0.5),
        legend.box="vertical",
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(10, "lines"))+
  guides(shape=guide_legend(nrow=1,byrow=TRUE),
         linetype=guide_legend(nrow=1,byrow=TRUE),
         color=guide_legend(nrow=1,byrow=TRUE))

print(ggplot.fit)
dev.off()






















