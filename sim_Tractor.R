beta_params <- function(m, v) {
  common_factor <- (m * (1 - m)) / v - 1
  alpha <- m * common_factor
  beta <- (1 - m) * common_factor
  return(c(alpha, beta))
}

nsim=1000
N=20000

certainty=c(0.5,0.6,0.7,0.8,0.9,0.95,0.98)
aveproba_all=c(0.1,0.25,0.5)

## consider different MAF thresholds: (maf1,maf2)={(0.01,0.05),(0.05,0.1),(0.1,0.2),(0.2,0.5)}

maf1=0.01  ## 0.05 0.1 0.2
maf2=0.05  ## 0.1 0.2 0.5

trutha=samplea=bestguessa=alproba=matrix(NA,nrow=nsim,ncol=2*21)
truthb=sampleb=bestguessb=alprobb=matrix(NA,nrow=nsim,ncol=2*21)
for(ap in 1:3){
  aveproba=aveproba_all[ap]
  for(cases in 1:7){
    cat("P(a)=",aveproba,": certainty=",certainty[cases],"\n")
    for(seed in 1:nsim){
      if(seed%%10==0) print(seed)
      set.seed(seed*cases)
      
      truepopa1_temp<-sapply(1:N,function(i){sample(c(1,0),1,prob=c(aveproba,1-aveproba))})
      truepopa2_temp<-sapply(1:N,function(i){sample(c(1,0),1,prob=c(aveproba,1-aveproba))})
      
      params <- beta_params(certainty[cases], (0.05/certainty[cases])^2)
      
      e1=rbeta(N, params[1], params[2])
      e2=rbeta(N, params[1], params[2])
      
      popa1=ifelse(truepopa1_temp == 1, e1, 1-e1)
      popa2=ifelse(truepopa2_temp == 1, e2, 1-e2)
      
      truepopa1<-sapply(1:N,function(i){sample(c(1,0),1,prob=c(popa1[i],1-popa1[i]))})
      truepopa2<-sapply(1:N,function(i){sample(c(1,0),1,prob=c(popa2[i],1-popa2[i]))})
      
      fst=0.2
      f1=sample(c(runif(1, maf1, maf2), runif(1, 1 - maf2, 1 - maf1)), 1)
      f2=sample(c(runif(1, maf1, maf2), runif(1, 1 - maf2, 1 - maf1)), 1)
      truef1=rbeta(1,f1*(1-fst)/fst,(1-f1)*(1-fst)/fst)
      truef2=rbeta(1,f2*(1-fst)/fst,(1-f2)*(1-fst)/fst)
      while(!((truef1 >= maf1 & truef1 <= maf2) | (truef1 >= (1 - maf2) & truef1 <= (1 - maf1)))){
        truef1=rbeta(1,f1*(1-fst)/fst,(1-f1)*(1-fst)/fst)
      }
      while(!((truef2 >= maf1 & truef2 <= maf2) | (truef2 >= (1 - maf2) & truef2 <= (1 - maf1)))){
        truef2=rbeta(1,f2*(1-fst)/fst,(1-f2)*(1-fst)/fst)
      }
      truef=c(truef1,truef2)
      generaw1=sapply(1:N,function(i){sample(c(1,0),1,prob=c(truef[truepopa1[i]+1],1-truef[truepopa1[i]+1]))})
      generaw2=sapply(1:N,function(i){sample(c(1,0),1,prob=c(truef[truepopa2[i]+1],1-truef[truepopa2[i]+1]))})
      
      popb1=1-popa1
      popb2=1-popa2
      
      
      truepopb1=1-truepopa1
      truepopb2=1-truepopa2
      
      alleleproba=(generaw1*popa1+generaw2*popa2)
      
      alleleprobb=(generaw1*popb1+generaw2*popb2)
      
      truedosagea=generaw1*truepopa1+generaw2*truepopa2
      truedosageb=generaw1*truepopb1+generaw2*truepopb2
      
      y=truedosagea*5+truedosageb*1+rnorm(N,0,15)
      
      gene=truedosagea+truedosageb
      
      lp=(popa1+popa2)/2
      
      truemodel=summary(lm(y~truedosagea+truedosageb+lp))
      trutha[seed,c(cases+7*(ap-1),21+cases+7*(ap-1))]=truemodel$coefficients["truedosagea",c(1,2)]
      truthb[seed,c(cases+7*(ap-1),21+cases+7*(ap-1))]=truemodel$coefficients["truedosageb",c(1,2)]
      
      alprobmodel=summary(lm(y~alleleproba+alleleprobb+lp))
      alproba[seed,c(cases+7*(ap-1),21+cases+7*(ap-1))]=alprobmodel$coefficients["alleleproba",c(1,2)]
      alprobb[seed,c(cases+7*(ap-1),21+cases+7*(ap-1))]=alprobmodel$coefficients["alleleprobb",c(1,2)]
      
      ##simulation sample
      set.seed(seed)
      samplepopa1<-sapply(1:N,function(i){sample(c(1,0),1,prob=c(popa1[i],1-popa1[i]))})
      samplepopa2<-sapply(1:N,function(i){sample(c(1,0),1,prob=c(popa2[i],1-popa2[i]))})
      samplepopb1=1-samplepopa1
      samplepopb2=1-samplepopa2
      
      sampledosagea=generaw1*samplepopa1+generaw2*samplepopa2
      sampledosageb=generaw1*samplepopb1+generaw2*samplepopb2
      
      samplemodel=summary(lm(y~sampledosagea+sampledosageb+lp))
      
      samplea[seed,c(cases+7*(ap-1),21+cases+7*(ap-1))]=samplemodel$coefficients["sampledosagea",c(1,2)]
      sampleb[seed,c(cases+7*(ap-1),21+cases+7*(ap-1))]=samplemodel$coefficients["sampledosageb",c(1,2)]
      
      ##best guess
      bgpopa1<-sapply(1:N,function(i){
        if(popa1[i]>=0.5){
          1
        }else{
          0
        }
      })
      
      bgpopa2<-sapply(1:N,function(i){
        if(popa2[i]>=0.5){
          1
        }else{
          0
        }
      })
      
      bgpopb1=1-bgpopa1
      bgpopb2=1-bgpopa2
      
      bgdosagea=generaw1*bgpopa1+generaw2*bgpopa2
      bgdosageb=generaw1*bgpopb1+generaw2*bgpopb2
      
      bgmodel=summary(lm(y~bgdosagea+bgdosageb+lp))
      
      bestguessa[seed,c(cases+7*(ap-1),21+cases+7*(ap-1))]=bgmodel$coefficients["bgdosagea",c(1,2)]
      bestguessb[seed,c(cases+7*(ap-1),21+cases+7*(ap-1))]=bgmodel$coefficients["bgdosageb",c(1,2)]
    }
  }
}

all=cbind(trutha,alproba,samplea,bestguessa,truthb,alprobb,sampleb,bestguessb)
write.csv(cbind(trutha,alproba,samplea,bestguessa,truthb,alprobb,sampleb,bestguessb),paste0("Tractor_compare_sim_",maf1,"-",maf2,".csv"))




