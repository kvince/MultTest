PCKG = c("doSNOW","foreach","data.table","sandwich","tseries","StepwiseTest")
for(pckg in PCKG) require(pckg,character.only = TRUE)


set.seed(123)

load("ew_l_mf")
nx = names(ew_l_mf)[-1]
N = nrow(ew_l_mf)
M = length(nx)

B = 1000  # number of bootstrap replication

load("riskfactors")
rf = riskfactors[["rf"]] * 100
ffps = as.matrix(riskfactors[date<"2016-01-01",.(mktrf,smb,hml,umd,ps_vwf)]) * 100

calc_alpha = function(portfolio_ret,factor_ret) {
  mdl = lm(portfolio_ret~factor_ret)
  v = sqrt(vcovHAC(mdl)[1,1])
  a = mdl$coefficients[1]
  return(c(a,v))
}
calc_alpha_only = function(portfolio_ret,factor_ret) {
  # calculate alpha only
  reg_result = lm(portfolio_ret ~ factor_ret)
  alpha = reg_result$coef[1]
  return(alpha)
}

# excess returns 
ew_l_mfx = as.matrix(ew_l_mf[,nx,with=F])*100 - rf

# Fama-French 3+Momentum+Illiquidity alphas
ew_l_mf_ffps = apply(ew_l_mfx[-((N-5):N),],2,calc_alpha,ffps)

ffps_hat_ew_l  = ew_l_mf_ffps[1,]/ew_l_mf_ffps[2,]

recenter_ffps_hat_ew_l =ffps_hat_ew_l * as.numeric(ffps_hat_ew_l < -(2*log(log(N-6)))^(1/2) ) 


ran_idx2=tsbootstrap(1:(N-6),nb=B,b=1/(1-0.5),type="stationary")
cl <- makeCluster(8)
registerDoSNOW(cl)
clusterEvalQ(cl, library(data.table))
pb <- txtProgressBar(max=B, style=3)
opts <- list(progress=function(n) setTxtProgressBar(pb, n))
FFPS_STAR<-foreach(boot_i=1:B,.combine=cbind,
                   .options.snow=opts) %dopar% {
  # Fama-French 3+Momentum+illiquidity alphas  
  boot_ew_l = ew_l_mfx[ran_idx2[,boot_i] , ]  
  boot_ffps = ffps[ran_idx2[,boot_i],]
  ew_l_ffps_star = apply(boot_ew_l,2,calc_alpha_only,boot_ffps)
  
  #  Compute the bootstrap statistics #
  ffps_star_ew_l  = ew_l_ffps_star/ew_l_mf_ffps[2,]
  
  ffps_boot_ew_l  = ffps_star_ew_l  - ffps_hat_ew_l  + 
    recenter_ffps_hat_ew_l 
  
}
close(pb)
stopCluster(cl)


MT_mf_Result = FDPControl(ffps_hat_ew_l,FFPS_STAR,0.05,0.05)




