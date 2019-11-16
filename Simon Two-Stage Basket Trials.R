# date: 10/04/2016
library(parallel) 
no.cores = 10 # uses 10 cores in parallel processing 

### specify inputs & verify defaults and grid search then run entire script ###

################################# INPUT ##################################

### target FWER & power & minimum acceptable power ###
FWER = 0.05; POWER = 0.8; min.POWER = 0.7

### number of baskets ###
No.B = 5

### scenario to calibrate for POWER ###
cali.diff.A = 0 # 0: use A = 2; 1: use A = 3 active	

### null and desired RR ###
Theta_0 = 0.15; Theta_a = 0.45

### accrual rates per basket ###
AR = rep(2,No.B)

### stage 1 sample size per basket (paper notation: n_1k) ###
N1 = 7

### stage 2 sample size per basket for homogeneous track (paper notation: N2/K) ###
N20 = 4

### futility rules for heterogeneous design (track 1) and homogeneous design(track 2) ###
Rs = 1 # track 1: min number of responders/basket required to move onto stage 2
Rc = No.B # track 2: min number of responders (across all baskets) required to move onto stage 2

############################# DEFAULTS ###################################

### FWER and power optimization margins ###
Emargin = 0.01; Pmargin = 0.02

### number of simulations ###
N.sim = 1000

### write out all results ###
out.all = 0 # 0: do not write.csv/table of all results; 1: write out all

############################# SEARCH SPACE ###############################

### stage 2 track 1 SS to search over ###
N21 =  12:22

### TOH tuning parameter, heterogeneous design path type 1 error, homogeneous design path type 1 error ###
gamma = seq(0.1,0.9,by = 0.05)
alphaS = seq(0.01,0.1,by = 0.01)
alphaC = seq(0.01,0.05, by = 0.01)

### search grid for 3 alphas: error rates ###
full.set.index = cbind(rep(1:length(gamma),each = length(alphaS) * length(alphaC)), rep(rep(1:length(alphaS),each = length(alphaC)),times = length(gamma)), rep(1:length(alphaC),times = length(gamma)*length(alphaS)))
for (it in 1:length(gamma)){
  Full.set = cbind(gamma[full.set.index[,1]],alphaS[full.set.index[,2]],alphaC[full.set.index[,3]])
}

################################# OUTPUT ##################################

### list output of full results ###
# scenario (number truly active)
# family wise error rate
# marginal rejection probabilities for each basket
# expected sample size (EN)
# proportion of trials that use heterogeneous design path (track 1)
# for each basket, proportion of trials that continue to stage 2 (in heterogeneous design path)
# proportion of trials that use homogeneous design path (track 2)
# proportion of trials that continue to stage 2 (in homogeneous design path)
# true positive rate
# false positive rate
# false negative rate
# true negative rate
# true negative rate (after stage 1)
# expected trial duration

### default output ###
# output txt file with 4 optimal design parameters: n21, gamma, alphaS, alphaC :
# e.g. "K5OptimalDesign_parameters_A_2_n1_7_N2_4.txt"
# output csv file with full results for the optimal design, assuming the 4 design parameters in the above
# e.g. "K5OptimalDesign_results_A_2_n1_7_N2_4.csv"

# for each N21, outputs full results for the best calibrated design (min S(.) - see Supplementary Materials):
# e.g. "K5CalibratedDesigns_AllResults_A_2_n1_7_N2_4_n2_12.csv"
# if no calibrated design exists for FWER within margin: keeps design with FWER closest to target FWER
# if no calibrated design exists for FWER & POWER within margin: selects design with highest power in 
# calibration alternative scenario (if multiple, uses design with smallest gamma)
# e.g. "K7MaxPowerCaliDesign_parameters_A_2_n1_14_N2_4_n2_12.csv" : outputs parameters & small summary
# if only 1 calibrated design exists for FWER & POWER outputs parameters and full results:
# e.g. "K7OnlyCaliDesign_parameters_A_2_n1_14_N2_4.txt"
# e.g. "K7OnlyCaliDesign_results_A_2_n1_14_N2_4.csv"

# output csv file with summary results for the best calibrated design for each N21
# e.g. "K5CalibratedDesigns_summary_A_2_n1_7_N2_4.csv"
# list of summary results
# columns : N21
# rows : gamma, alphaS, alphaC, family wise error rate (0 active), expected sample size (0 active),
# marginal power for Basket 1 and expected sample size for alternative scenarios, i.e., 1, 2,..., K active, 
# value of utility function U(.) - see Supplementary Materials 

### additional output ###
# if you change "out.all = 1":
# for each scenario, outputs full results for each N21 for entire grid search of 3 alphas

################################################################################################################

trial = function(it, AR, K, n1.goal, n2.1.goal, n2.0.goal, theta_0, theta_a, p.sim, gamma, alphaS, alphaC, rC, rS){ 
  
  # stage 1 min and max sample sizes before can perform TOH
  s1.max = ceiling(1.5*n1.goal) 
  no.bask.rule1 = n1.goal*K
  s1.min = floor(0.5*n1.goal)
  
  # stage 2|homogeneous design min and max sample sizes before can perform one-sample test of efficacy
  s2.max = ceiling(1.5*n2.0.goal)
  s2.min = 1
  no.bask.rule2 = n2.0.goal*K
  
  Time1 = {}
  ### generate time of enrollment for stage 1###
  Baskets = matrix(NA,nrow = s1.max*K,ncol = 2)
  for (i in 1:K){
    Baskets[(i+(i-1)*(s1.max-1)):(i*s1.max),] = cbind(rexp(s1.max,AR[i]),rep(i,s1.max))
  }
  ord = order(Baskets[,1])
  keep = Baskets[ord,][1:no.bask.rule1,]
  tab = table(factor(keep[,2],levels = 1:K))
  n1 = {}
  for (i in 1:K){
    n1[i] = as.numeric(tab)[i]
    if (n1[i] == 0) {Time1[i] = 0}
    if (n1[i] == 1) {Time1[i] = sum(keep[which(keep[,2]==i),][1])}
    if (n1[i] >= 2) {Time1[i] = sum(keep[which(keep[,2]==i),][,1])}
    if (n1[i] < s1.min){Time1[i] = Time1[i] + sum(rexp(s1.min-n1[i],AR[i])); n1[i] = s1.min}
  }
  ### generate stage 1 responses ###
  yes.s1 = rbinom(K,n1,p.sim)
  no.s1 = n1 - yes.s1
  
  ### evaluate TOH ###
  tab = matrix(c(yes.s1,no.s1),nrow = K,byrow = F)
  p.val = ifelse(sum(yes.s1)==0, fisher.test(tab)$p.value, fisher.test(tab,hybrid = T,simulate.p.value=T)$p.value) # don't simulate pval if marginal of all yes = 0
  toh = as.numeric(p.val <= gamma)
  
  yes.s2 = rep(NA,K); dec = rep(0,K); K.star = {}; stage2.t2 = 0; stage2.t1 = rep(0,K); n2.0 = rep(NA,K); n2.1 = rep(NA,K)
  Time21 = rep(0,K); Time20 = rep(0,K)
  
  if (toh == 1){ # reject TOH => heterogeneous baskets
    
    ### determine which baskets to keep ###
    K.star = which(yes.s1 >= rS) # keep baskets with min desirable RR
    stage2.t1[K.star] = 1 # keep track of baskets that go on
    
    if (length(K.star) > 0){
      n2.1[K.star] = n2.1.goal
      
      ### generate time of enrollment for stage 2|heterogeneous baskets###
      temp.Time21 = matrix(NA,nrow = n2.1.goal,ncol = K)
      for (i in 1:length(K.star)){
        temp.Time21[,K.star[i]] = rexp(n2.1[K.star[i]],AR[K.star[i]])
        Time21[K.star[i]] = sum(temp.Time21[,K.star[i]])
      }
      ### stage 2 ###
      yes.s2[K.star] = rbinom(length(K.star),n2.1[K.star],p.sim[K.star])
      no.s2 = n2.1 - yes.s2
      
      ### decision to reject in each basket ###
      for (i in 1:length(K.star)){dec[K.star[i]] = as.numeric(binom.test(yes.s1[K.star[i]]+yes.s2[K.star[i]],n1[K.star[i]]+n2.1[K.star[i]],theta_0,alternative = "greater")$p.value <= (alphaS/length(K.star)) ) }
    } # end K.star
  }else if (toh == 0) {
    
    ### generate time of enrollment for stage 2|homogeneous baskets###
    Baskets2 = matrix(NA,nrow = s2.max*K,ncol = 2)
    for (i in 1:K){
      Baskets2[(i+(i-1)*(s2.max-1)):(i*s2.max),] = cbind(rexp(s2.max,AR[i]),rep(i,s2.max))
    }
    ord = order(Baskets2[,1])
    keep2 = Baskets2[ord,][1:no.bask.rule2,]
    tab = table(factor(keep2[,2],levels = 1:K))
    
    # keeping track if baskets go on to stage 2
    stage2.t2 = ifelse(sum(yes.s1) >= rC,1,0)
    
    if (stage2.t2 == 1){
      for (i in 1:K){
        n2.0[i] = as.numeric(tab)[i]
        if (n2.0[i] >= 2) {Time20[i] = sum(keep2[which(keep2[,2]==i),][,1])}
        if (n2.0[i] == 1) {Time20[i] = sum(keep2[which(keep2[,2]==i),][1])}
        if (n2.0[i] == 0) {Time20[i] = 0}
        if (n2.0[i] < s2.min){Time20[i] =  Time20[i] + sum(rexp(s2.min-n2.0[i],AR[i])); n2.0[i] = s2.min}
      }
      ### stage 2 ###
      yes.s2 = rbinom(K,n2.0,p.sim)
      no.s2 = n2.0 - yes.s2
      
      ### decision to reject one sample ###
      dec = rep(ifelse(as.numeric(binom.test(sum(yes.s1+yes.s2),sum(n1+n2.0),theta_0,alternative = "greater")$p.value <= (alphaC) ),1,0),K)
    } # end stage2.t2 	
  } # end toh
  
  Time = max(Time1) + max(Time21,na.rm = T) + max(Time20,na.rm = T)
  
  return(c(dec, yes.s1, yes.s2, length(K.star), toh, stage2.t1,stage2.t2,n1,n2.0,Time))
} # end sim func

################################################################################################################

sim = function(sc, X, n.sim, K, n1, n2.1, n2.0, theta_0, theta_a, r.s, r.c, accrual){
  h = sc
  set.seed(4206)
  one.result = sapply(1:n.sim, trial, K = K, n1.goal = n1, n2.1.goal = n2.1, n2.0.goal = n2.0, theta_0 = theta_0, theta_a = theta_a, rS = r.s, rC = r.c, AR = accrual, gamma = X[1], alphaS = X[2], alphaC = X[3], p.sim = c(rep(theta_a,h), rep(theta_0,K-h)))
  ### Marginal Rejection Probability ###
  mrp = rowSums(one.result[(1:(K)),],na.rm = T)/n.sim
  ### Overall FWER ###
  fwer = ifelse(h < K - 1, sum(as.numeric(colSums(one.result[((1+h):(K)),],na.rm = T) >= 1))/n.sim, ifelse( h == K-1, sum(one.result[((1+h):(K)),],na.rm = T)/n.sim, NA ) ) # only one basket under the null
  ### number of trials each track ###
  toh1 = which(one.result[(3*K+2),] == 1) # Track 1
  toh0 = which(one.result[(3*K+2),] == 0) # Track 2
  results1 = one.result[,toh1]
  results0 = one.result[,toh0]
  track1 = length(toh1); track2 = length(toh0)
  ### count of stage 2 for each track ###
  stage2.t1 = rowSums(one.result[(3*K+3):(4*K+2),],na.rm = T)
  stage2.t2 = sum(one.result[(4*K+3),],na.rm = T)
  ### number of baskets that go on in Track 1 ###
  if (length(toh1) >1){results1.Ktab =table(factor(as.numeric(results1[(3*K+1),]),levels = 0:K))
  }else{results1.Ktab =table(factor(as.numeric(results1[(3*K+1)]),levels = 0:K))}
  ### calc EN ###
  n1.obs = rowSums(one.result[(4*K+4):(5*K+3),])/n.sim
  n2.0.obs = rowSums(one.result[(5*K+4):(6*K+3),],na.rm = T)/n.sim
  temp = n2.1*(stage2.t1/track1)
  EN = sum(n1.obs)+(track1/n.sim)*sum(temp)+(track2/n.sim)*( sum(n2.0.obs)*(stage2.t2/track2) )
  ### sensitivity/specificity ###
  pp = round(ifelse(h > 0,sum(one.result[1:h,]==1,na.rm = T)/(K*n.sim), NA),3)
  # truth = no; decision = reject
  np = round(ifelse(h < K,sum(one.result[(h+1):K,]==1,na.rm = T)/(K*n.sim), NA),3)
  # truth = yes; decision = dnr
  pn = round(ifelse(h > 0,sum(one.result[1:h,]==0,na.rm = T)/(K*n.sim), NA),3)
  nn = round(ifelse(h < K,sum(one.result[(h+1):K,]==0,na.rm = T)/(K*n.sim), NA),3)
  nn.1 = round(ifelse(h < K,(sum(one.result[(3*K+3+h):(4*K+2),]==0,na.rm = T)+sum(one.result[(4*K+3),]==0,na.rm = T))/(K*n.sim), NA),3)
  time = round(sum(one.result[(6*K+4),])/n.sim,1)
  return(c(h,fwer,mrp,EN,track1/n.sim,stage2.t1/n.sim,track2/n.sim,stage2.t2/n.sim,pp,np,pn,nn,nn.1,time))
}

################################################################################################################
big.list = list();trunc.ind = list(); trunc.grid = list();trunc.ind2 = list(); trunc.grid2 = list(); results.tab1 = list(); results.tab2 = list(); results.tab3 = list()
new.trunc.ind = list(); N21.keepind = {}

start = Sys.time()

### search over track 1 stage 2 sample size ###
for (i in 1:length(N21)){
  print(i)
  ### NULL SCENARIO ###
  results = mclapply(1:dim(Full.set)[1], function(j, fullset, sc, n.sim, K, n1, n2.1, n2.0, theta_0, theta_a, r.s, r.c, accrual){sim(X = fullset[j,], sc, n.sim, K, n1, n2.1, n2.0, theta_0, theta_a, r.s, r.c, accrual)}, fullset = Full.set, n.sim = N.sim, K = No.B, n1 = N1, n2.1 =  N21[i], n2.0 = N20, theta_0 = Theta_0, theta_a = Theta_a, r.s = Rs, r.c = Rc, accrual = AR, sc = 0, mc.silent = TRUE, mc.cores = no.cores, mc.preschedule = FALSE) 
  
  results.tab = matrix(data = unlist(results),ncol = dim(Full.set)[1], nrow = 12+2*No.B)
  big.list[[i]] = results.tab
  if (out.all == 1){ write.csv(results.tab,paste("K",No.B,"_Sc0_results_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,"_n2_",N21[i],".csv",sep = ""),row.names = F) }
  
  ## select alpha combos with desirable FWER ##
  a = abs(as.numeric(big.list[[i]][2,]) - FWER) 
  trunc.ind[[i]] = which( a <= Emargin)
  # make sure function returns something (if Emargin too small)
  if (length(trunc.ind[[i]])==0){trunc.ind[[i]] = which(a == min(a))}
  
  trunc.grid[[i]] = Full.set[trunc.ind[[i]],]
  if(length(trunc.ind[[i]])==1){
    upper.ind =  1
    trunc.grid[[i]] = matrix(trunc.grid[[i]],nrow = 1,ncol = 3)
  }else{ upper.ind = dim(trunc.grid[[i]])[1] } 
  
  if (out.all == 1){ write.table(trunc.ind[[i]],paste("K",No.B,"TruncatedIndex_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,"_n2_",N21[i],".txt",sep = ""))	}
  
  if (cali.diff.A == 0){
    ### SCENARIO 2 ###
    results = mclapply(1:upper.ind, function(j, fullset, sc, n.sim, K, n1, n2.1, n2.0, theta_0, theta_a, r.s, r.c, accrual){sim(X = fullset[j,], sc, n.sim, K, n1, n2.1, n2.0, theta_0, theta_a, r.s, r.c, accrual)}, fullset = trunc.grid[[i]], n.sim = N.sim, K = No.B, n1 = N1, n2.1 =  N21[i], n2.0 = N20, theta_0 = Theta_0, theta_a = Theta_a, r.s = Rs, r.c = Rc, accrual = AR, sc = 2, mc.silent = TRUE, mc.cores = no.cores, mc.preschedule = FALSE) 	
    results.tab2[[i]] = matrix(data = unlist(results),ncol = upper.ind, nrow = 12+2*No.B)
    if (out.all == 1){ write.csv(results.tab2[[i]],paste("K",No.B,"_Sc2_results_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,"_n2_",N21[i],".csv",sep = ""),row.names = F) }
  }
  
  if (cali.diff.A == 1){
    # calibrate POWER at A = 3, for larger K (more baskets)
    ### SCENARIO 3 ###
    results = mclapply(1:upper.ind, function(j, fullset, sc, n.sim, K, n1, n2.1, n2.0, theta_0, theta_a, r.s, r.c, accrual){sim(X = fullset[j,], sc, n.sim, K, n1, n2.1, n2.0, theta_0, theta_a, r.s, r.c, accrual)}, fullset = trunc.grid[[i]], n.sim = N.sim, K = No.B, n1 = N1, n2.1 =  N21[i], n2.0 = N20, theta_0 = Theta_0, theta_a = Theta_a, r.s = Rs, r.c = Rc, accrual = AR, sc = 3, mc.silent = TRUE, mc.cores = no.cores, mc.preschedule = FALSE) 	
    results.tab2[[i]] = matrix(data = unlist(results),ncol = upper.ind, nrow = 12+2*No.B)
    if (out.all == 1){ write.csv(results.tab2[[i]],paste("K",No.B,"_Sc3_results_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,"_n2_",N21[i],".csv",sep = ""),row.names = F) }
  }
  
  if ( any(results.tab2[[i]][3,] >= (POWER - Pmargin)) ){
    N21.keepind = c(N21.keepind,i)
    ### SCENARIO 1 ### 
    results = mclapply(1:upper.ind, function(j, fullset, sc, n.sim, K, n1, n2.1, n2.0, theta_0, theta_a, r.s, r.c, accrual){sim(X = fullset[j,], sc, n.sim, K, n1, n2.1, n2.0, theta_0, theta_a, r.s, r.c, accrual)}, fullset = trunc.grid[[i]], n.sim = N.sim, K = No.B, n1 = N1, n2.1 =  N21[i], n2.0 = N20, theta_0 = Theta_0, theta_a = Theta_a, r.s = Rs, r.c = Rc, accrual = AR, sc = 1, mc.silent = TRUE, mc.cores = no.cores, mc.preschedule = FALSE) 
    results.tab1[[i]] = matrix(data = unlist(results),ncol = upper.ind, nrow = 12+2*No.B)
    if (out.all == 1){ write.csv(results.tab1[[i]],paste("K",No.B,"_Sc1_results_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,"_n2_",N21[i],".csv",sep = ""),row.names = F) }
  }
  
}# end i

warnings()

Sys.time() - start

#### SELECT BEST DESIGN for each N21[i] ####
calibrate = function(des,x1,x2,x3){
  return(abs(des[4]-x1)+abs(des[6]-x2)+abs(des[8]-x3))
}

select = function(des,diff.A){
  if(diff.A == 0){x = seq(10,3+(No.B+1)*2,2)
  }else{x = seq(12,3+(No.B+1)*2,2)}
  y = seq(5,3+(No.B+1)*2,2)
  return(mean(des[x])*100-mean(des[y]))
}

des.para = matrix(NA, nrow = 3+(No.B+1)*2, ncol = length(N21))
Best.OCs = list()

for (i in 1:length(N21)){
  
  if ( i %in% N21.keepind ){
    
    B = abs(as.numeric(results.tab1[[i]][3,]) - min.POWER)
    C = abs(as.numeric(results.tab2[[i]][3, ]) - POWER)
    temp = which( B <= Pmargin & C <= Pmargin )
    
    if (length(temp) == 0) { 
      temp.ind = which(C <= Pmargin)
      
      if (length(temp.ind) > 0 ){
        B = as.numeric(results.tab1[[i]][3,]) - min.POWER
        B[B <= -Pmargin] = max(B)+1 # 1 is arbitrary - so will select a design with at least min.power within Pmargin
        temp = temp.ind[which.min( B[temp.ind] ) ] # picks first duplicate - has smallest gamma
      }else{ 
        max.pow2 = max(results.tab2[[i]][3,])
        temp.ind = which(results.tab2[[i]][3,] == max.pow2)
        min.pow1 = min(B[temp.ind])
        keep = temp.ind[which(B[temp.ind]==min.pow1)]
        min.alpha = min(Full.set[trunc.ind[[i]][keep],1]) 
        
        if (any(duplicated(Full.set[trunc.ind[[i]][keep],1] ) ) ){
          dup.alpha1 = keep[which(Full.set[trunc.ind[[i]][keep],1] == min.alpha)]
          temp = dup.alpha1[length(dup.alpha1)]
        }else{
          temp = keep[which(Full.set[trunc.ind[[i]][keep],1] == min.alpha)] # keep smallest alpha, since all have same min power in S2
        } # any dup	
      } # temp.ind
      
      designs = c(Full.set[trunc.ind[[i]][temp],],as.numeric(big.list[[i]][2,trunc.ind[[i]][temp]]),round(as.numeric(big.list[[i]][3+No.B,trunc.ind[[i]][temp] ])),as.numeric(results.tab1[[i]][3, temp ]),round(as.numeric(results.tab1[[i]][3+No.B, temp ])),as.numeric(results.tab2[[i]][3,temp]),round(as.numeric(results.tab2[[i]][3+No.B,temp])))
      if (out.all == 1){ write.csv(designs,paste("K",No.B,"TruncatedDesigns_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,"_n2_",N21[i],".csv",sep = ""),row.names = F) }
      
      i.choose.you = c(Full.set[trunc.ind[[i]][temp],],as.numeric(big.list[[i]][2,trunc.ind[[i]][temp]]),round(as.numeric(big.list[[i]][3+No.B,trunc.ind[[i]][temp] ])),as.numeric(results.tab1[[i]][3, temp ]),round(as.numeric(results.tab1[[i]][3+No.B, temp ])),as.numeric(results.tab2[[i]][3,temp]),round(as.numeric(results.tab2[[i]][3+No.B,temp])))
      
    }else if(length(temp) == 1){
      
      i.choose.you = c(Full.set[trunc.ind[[i]][temp],],as.numeric(big.list[[i]][2,trunc.ind[[i]][temp]]),round(as.numeric(big.list[[i]][3+No.B,trunc.ind[[i]][temp] ])),as.numeric(results.tab1[[i]][3, temp ]),round(as.numeric(results.tab1[[i]][3+No.B, temp ])),as.numeric(results.tab2[[i]][3,temp]),round(as.numeric(results.tab2[[i]][3+No.B,temp])))
      if (out.all == 1){ write.csv(designs,paste("K",No.B,"TruncatedDesigns_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,"_n2_",N21[i],".csv",sep = ""),row.names = F) }
      
    }else{
      designs = cbind(Full.set[trunc.ind[[i]][temp],],as.numeric(big.list[[i]][2,trunc.ind[[i]][temp]]),round(as.numeric(big.list[[i]][3+No.B,trunc.ind[[i]][temp] ])),as.numeric(results.tab1[[i]][3, temp ]),round(as.numeric(results.tab1[[i]][3+No.B, temp ])),as.numeric(results.tab2[[i]][3,temp]),round(as.numeric(results.tab2[[i]][3+No.B,temp])))
      if (out.all == 1){ write.csv(designs,paste("K",No.B,"TruncatedDesigns_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,"_n2_",N21[i],".csv",sep = ""),row.names = F) }
      
      ord = order(C[temp])
      ord.designs = designs[ord,]
      
      cali.designs = apply(ord.designs,1,calibrate,x1 = FWER, x2 = min.POWER, x3 = POWER)
      min.cali = which(cali.designs == min(cali.designs))
      
      if (length(min.cali) == 1){
        i.choose.you = ord.designs[min.cali,]
      }else{
        best.designs = ord.designs[min.cali,]
        min.alpha1 = min(best.designs[,1])
        dup.alpha1 = best.designs[which(best.designs[,1] == min.alpha1),]
        
        ## for n2 = N21[i], select best design ##
        if(sum(best.designs[,1] == min.alpha1) > 1 ){
          i.choose.you = dup.alpha1[nrow(dup.alpha1),]
        }else{i.choose.you = dup.alpha1}
      }  # min.cali
    } # length temp
    
    des.para[1:9,i] = i.choose.you
    Best.OCs[[i]] = sapply(0:No.B, function(sc, set, n.sim, K, n1, n2.1, n2.0, theta_0, theta_a, r.s, r.c, accrual){sim(sc, X = set, n.sim, K, n1, n2.1, n2.0, theta_0, theta_a, r.s, r.c, accrual)}, set = i.choose.you[1:3], n.sim = N.sim, K = No.B, n1 = N1, n2.1=  N21[i], n2.0 = N20, theta_0 = Theta_0, theta_a = Theta_a, r.s = Rs, r.c = Rc, accrual = AR)	
    Best.OCs[[i]][3+No.B,] = round(Best.OCs[[i]][3+No.B,])
    rownames(Best.OCs[[i]]) = c("No.Active","FWER",paste("P",1:No.B,sep = ""),"EN","Track1",paste("Stage2|Track1 B",1:No.B,sep = ""),"Track2","Stage2|Track2","TP","FP","FN","TN","TN(stage1)","ET")
    write.csv(Best.OCs[[i]],paste("K",No.B,"CalibratedDesigns_AllResults_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,"_n2_",N21[i],".csv",sep = ""))
    parameters = i.choose.you[1:3]
    names(parameters) = c("gamma","alphaS","alphaC")
    write.table(parameters,paste("K",No.B,"CalibratedDesigns_parameters_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,"_n2_",N21[i],".txt",sep = ""),row.names = T)
    
    if (cali.diff.A == 0){ 
      des.para[10:(3+(No.B+1)*2),i] = c(Best.OCs[[i]][c(3,3+No.B),4:(No.B+1)]) 
    }else{
      des.para[8:(3+(No.B+1)*2),i] = c(Best.OCs[[i]][c(3,3+No.B),3:(No.B+1)])
    } # cali.diff.A
    
  }else if( !(i %in% N21.keepind) ){
    
    max.pow2 = max(results.tab2[[i]][3,])
    keep = which(results.tab2[[i]][3,] == max.pow2)
    min.alpha = min(Full.set[trunc.ind[[i]][keep],1]) 
    
    if (any(duplicated(Full.set[trunc.ind[[i]][keep],1] ) ) ){
      dup.alpha1 = keep[which(Full.set[trunc.ind[[i]][keep],1] == min.alpha)]
      temp = dup.alpha1[length(dup.alpha1)]
    }else{
      temp = keep[which(Full.set[trunc.ind[[i]][keep],1] == min.alpha)] # keep smallest alpha, since all have same min power in S2
    } # any dup	
    
    TopDesign = c(Full.set[trunc.ind[[i]][temp],],as.numeric(big.list[[i]][2,trunc.ind[[i]][temp]]),round(as.numeric(big.list[[i]][3+No.B,trunc.ind[[i]][temp] ])),as.numeric(results.tab2[[i]][3,temp]),round(as.numeric(results.tab2[[i]][3+No.B,temp])))
    names(TopDesign) = c("gamma","alphaS","alphaC","FWER","EN", paste("P1 (A = ",2*(1-cali.diff.A) + 3*cali.diff.A ,")",sep = ""),"EN")
    out = cbind(TopDesign,NA)
    write.table(out,paste("K",No.B,"MaxPowerCaliDesign_parameters_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,"_n2_",N21[i],".txt",sep = ""),na = "")
    
    if (out.all == 1){ 
      Best.OCs[[i]] = sapply(0:No.B, function(sc, set, n.sim, K, n1, n2.1, n2.0, theta_0, theta_a, r.s, r.c, accrual){sim(sc, X = set, n.sim, K, n1, n2.1, n2.0, theta_0, theta_a, r.s, r.c, accrual)}, set = Full.set[trunc.ind[[i]][temp],], n.sim = N.sim, K = No.B, n1 = N1, n2.1=  N21[i], n2.0 = N20, theta_0 = Theta_0, theta_a = Theta_a, r.s = Rs, r.c = Rc, accrual = AR)	
      Best.OCs[[i]][3+No.B,] = round(Best.OCs[[i]][3+No.B,])
      rownames(Best.OCs[[i]]) = c("No.Active","FWER",paste("P",1:No.B,sep = ""),"EN","Track1",paste("Stage2|Track1 B",1:No.B,sep = ""),"Track2","Stage2|Track2","TP","FP","FN","TN","TN(stage1)","ET")
      write.csv(Best.OCs[[i]],paste("K",No.B,"MaxPowerCaliDesign_results_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,"_n2_",N21[i],".csv",sep = ""))
    } # out.all
    
  } # any i
} # i

if (length(N21.keepind) == 1){
  final.combo = c(N21[N21.keepind],des.para[1:3,N21.keepind])
  names(final.combo) = c("n2","gamma","alphaS","alphaC")
  
  write.table(final.combo,paste("K",No.B,"OnlyCaliDesign_parameters_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,".txt",sep = ""),row.names = F,col.names = F)
  write.csv(Best.OCs[[N21.keepind]],paste("K",No.B,"OnlyCaliDesign_results_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,".csv",sep = ""))
  
}else if (length(N21.keepind) > 1){
  select.design = apply(des.para,2,select, diff.A = cali.diff.A)
  final.combo = c(N21[which.max(select.design)],des.para[1:3,which.max(select.design)])
  names(final.combo) = c("n2","gamma","alphaS","alphaC")
}

if (length(N21.keepind) > 1){
  write.table(final.combo,paste("K",No.B,"OptimalDesign_parameters_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,".txt",sep = ""),row.names = T,col.names = F)
  write.csv(Best.OCs[[which.max(select.design)]],paste("K",No.B,"OptimalDesign_results_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,".csv",sep = ""))
  
  des.para = rbind(des.para, select.design)
  temp.names = {}
  for (i in 1:No.B){ temp.names = c(temp.names,paste("P1 (A = ",i,")",sep = ""),paste("EN (A = ",i,")",sep = "")) }
  rownames(des.para) = c("gamma","alphaS","alphaC","FWER","EN(A = 0)",temp.names,"U")
  colnames(des.para) = N21
  write.csv(des.para,paste("K",No.B,"CalibratedDesigns_summary_A_",2*(1-cali.diff.A)+3*cali.diff.A,"_n1_",N1,"_N2_",N20,".csv",sep = ""),row.names = T)
}

Sys.time() - start













