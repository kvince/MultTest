PCKG = c("doSNOW","foreach","tseries","StepwiseTest")
for(pckg in PCKG) require(pckg,character.only = TRUE)

n=500        
B=1000        
fwer_level=0.05 
k = 3
NSIM = 1000

m_plus = 20
m_o    = 70
m_minus= 30
m = m_plus+m_o+m_minus
mu_plus  =  0.3
mu_minus = -0.3

mu = c( seq(0.01,mu_plus,length.out=m_plus), 
        rep(0,m_o), 
        seq(-0.01,mu_minus,length.out=m_minus)
      )

#sig = diag( 1 + 2*sqrt(abs(mu)) );
sig = diag(m)
rho = 0.2;
omega= sig %*% ((1-rho)*diag(1,m) + rho*matrix(1,m,m)) %*% sig 
v=t(chol(omega))


false_reject_rc = power_rc = cv_rc = NULL
false_reject_rc_k = power_rc_k = cv_rc_k = NULL
false_reject_spa = power_spa = cv_spa = NULL
false_reject_spa_k = power_spa_k = cv_spa_k = NULL
for(sim in 1:NSIM) {
  error=matrix(rnorm(m*n),m,n)      
  y = mu%*%matrix(1,1,n)+ v %*% error 
  
y_mean = apply(y,1,mean)
y_sig = apply(y,1,sd)
recenter_mean = y_mean * ifelse( sqrt(n)*y_mean < -y_sig*sqrt(2*log(log(n))), 1, 0)
s = tsbootstrap(1:n,B,b=2,type="stationary")
y_mean_boot_rc = y_mean_boot_spa = matrix(NA,m,B)
for(i in 1:B) {
  y_boot = y[, s[,i]]
  y_mean_boot = apply(y_boot,1,mean)
  y_mean_boot_rc[,i] = y_mean_boot - y_mean
  y_mean_boot_spa[,i] = y_mean_boot - y_mean + recenter_mean
}

res_rc = FWERkControl(y_mean,y_mean_boot_rc,1,fwer_level)
res_rc_k = FWERkControl(y_mean,y_mean_boot_rc,k,fwer_level)

res_spa = FWERkControl(y_mean,y_mean_boot_spa,1,fwer_level)
res_spa_k = FWERkControl(y_mean,y_mean_boot_spa,k,fwer_level)

false_reject_rc[sim] = as.numeric(sum(res_rc$Reject[(m_plus+1) : m]) >= 1 )
power_rc[sim] = mean(res_rc$Reject[1:m_plus])
cv_rc[sim] = res_rc$CV

false_reject_rc_k[sim] = as.numeric(sum(res_rc_k$Reject[(m_plus+1) : m]) >= k )
power_rc_k[sim] = mean(res_rc_k$Reject[1:m_plus])
cv_rc_k[sim] = res_rc_k$CV

false_reject_spa[sim] = as.numeric(sum(res_spa$Reject[(m_plus+1) : m]) >= 1 )
power_spa[sim] = mean(res_spa$Reject[1:m_plus])
cv_spa[sim] = res_spa$CV

false_reject_spa_k[sim] = as.numeric(sum(res_spa_k$Reject[(m_plus+1) : m]) >= k)
power_spa_k[sim] = mean(res_spa_k$Reject[1:m_plus])
cv_spa_k[sim] = res_spa_k$CV
}

mean(false_reject_rc)
mean(false_reject_rc_k)
mean(false_reject_spa)
mean(false_reject_spa_k)

mean(power_rc)
mean(power_rc_k)
mean(power_spa)
mean(power_spa_k)

