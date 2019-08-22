PCKG = c("doSNOW","foreach","data.table","Rcpp","tseries","StepwiseTest")
for(pckg in PCKG) require(pckg,character.only = TRUE)
sourceCpp('BasicTTR.cpp')

load("TTR_Data")
Price = TTR_Data[,DJIA]
Rfrate = log(1+TTR_Data[,RF]/360)

# Filter rules
initiate_x = c(0.005,0.01,0.02,0.03,0.04,0.05,0.1)
liquidate_y = c(0.005,0.01,0.02,0.03,0.04,0.05)
xy_comb = merge(initiate_x,liquidate_y)
xy_comb = xy_comb[which(xy_comb[,2] < xy_comb[,1]),]
FR = mapply(Filter_Rule, x=xy_comb[,1], y=xy_comb[,2], 
             MoreArgs = list(P=Price))
colnames(FR) = apply(xy_comb,1,function(z) 
                  paste("fr","_x",z[1],"_y",z[2],sep=""))

# MA crossovers
window = c(1,5,10,20,50,100,200)
ma_comb = merge(window,window)
ma_comb = ma_comb[ma_comb[,1]<ma_comb[,2],]
MA_b0 = mapply(MA_Crossover, fast_n=ma_comb[,1], slow_n=ma_comb[,2], 
            MoreArgs = list(b=0.0,P=Price))
colnames(MA_b0) = apply(ma_comb,1,function(z) 
  paste("ma_b0","_f",z[1],"_s",z[2],sep=""))
MA_b1 = mapply(MA_Crossover, fast_n=ma_comb[,1], slow_n=ma_comb[,2], 
               MoreArgs = list(b=0.01,P=Price))
colnames(MA_b1) = apply(ma_comb,1,function(z) 
  paste("ma_b1","_f",z[1],"_s",z[2],sep=""))

# Channel breakouts
rolling_window = c(5,10,15,20) # win
channel_width = c(0.05,0.1) # x
breakout_size = c(0.001,0.005,0.01) # b
holding_period = c(10,25,50) # c
cb_comb = merge(channel_width,breakout_size)
cb_comb = merge(rolling_window,cb_comb,by=NULL) 
cb_comb = merge(holding_period,cb_comb,by=NULL)
CB = mapply(Channel_Breakout,c=cb_comb[,1],win=cb_comb[,2],
              x=cb_comb[,3],b=cb_comb[,4],MoreArgs = list(P=Price))
colnames(CB) = paste("cb","_c",cb_comb[,1],"_n",cb_comb[,2],
                         "_x",cb_comb[,3],"_b",cb_comb[,4],sep="")

# RSI
h = c(10,20,30)
v = c(20,25,30)
d = c(1,2,3)
rsi_comb = merge(h,merge(v,d),by=NULL)
RSI = mapply(RSI_Rule,h=rsi_comb[,1],v=rsi_comb[,2,],d=rsi_comb[,3],
              MoreArgs = list(P=Price) )
colnames(RSI) = paste("rsi","_h",rsi_comb[,1],"_v",rsi_comb[,2],
                        "_d",rsi_comb[,3],sep="")

AllTTR =do.call(cbind,list(FR,MA_b0,MA_b1,CB,RSI) )


N = length(Price)
exrets = log(Price[2:N]/Price[1:(N-1)]) - Rfrate[1:(N-1)]
dt=data.table(ExRets=exrets,AllTTR[1:(N-1),])
cols=colnames(AllTTR)
TTR_exrets=dt[,lapply(.SD, function(x) x*ExRets ),.SDcols=cols]
TTR_exrets_tim=dt[,lapply(.SD, function(x) (x-mean(x))*ExRets),.SDcols=cols]

mret=unlist(TTR_exrets[,lapply(.SD,mean),.SDcols=cols])
sret=unlist(TTR_exrets[,lapply(.SD,sd),.SDcols=cols])

mret_tim=unlist(TTR_exrets_tim[,lapply(.SD,mean),.SDcols=cols])
sret_tim=unlist(TTR_exrets_tim[,lapply(.SD,sd),.SDcols=cols])

SR = mret/sret*sqrt(252)
SR_tim = mret_tim/sret_tim*sqrt(252)

B=1000
ran_idx=tsbootstrap(1:(N-1),nb=B,b=1/(1-0.9),type="stationary")
cl <- makeCluster(8)
registerDoSNOW(cl)
clusterEvalQ(cl, library(data.table))
pb <- txtProgressBar(max=B, style=3)
opts <- list(progress=function(n) setTxtProgressBar(pb, n))
comb <- function(...) {
  mapply('cbind', ..., SIMPLIFY=FALSE)
}
boot_SR<-foreach(boot_i=1:B,.combine='comb',.multicombine=TRUE,
                        .options.snow=opts) %dopar% {
                          
  boot_exrets = TTR_exrets[ran_idx[,boot_i] , ]   
  boot_exrets_tim = TTR_exrets_tim[ran_idx[,boot_i] , ]   
  mret=unlist(boot_exrets[,lapply(.SD,mean),.SDcols=cols])
  sret=unlist(boot_exrets[,lapply(.SD,sd),.SDcols=cols])
  mret_tim=unlist(boot_exrets_tim[,lapply(.SD,mean),.SDcols=cols])
  sret_tim=unlist(boot_exrets_tim[,lapply(.SD,sd),.SDcols=cols])
  SR_boot = (mret/sret*sqrt(252) - SR)
  SR_tim_boot = (mret_tim/sret_tim*sqrt(252) - SR_tim)
  
  list(SR_boot,SR_tim_boot)
}
close(pb)
stopCluster(cl)


TestRes = FWERkControl(SR,boot_SR[[1]],3,0.05)
TestRes_tim = FWERkControl(SR_tim,boot_SR[[2]],3,0.05)

FDPTestRes = FDPControl(SR,boot_SR[[1]],0.5,0.1)
FDPTestRes_tim = FDPControl(SR_tim,boot_SR[[2]],0.5,0.1)
