library(openxlsx)
library(actuar)
library(psych)
library(MASS)
library(mnormt)
library(mFilter)
library(tseries)
library(urca)
library(forecast)
library(tmvtnorm)
library(zoo)           
library(xts)            
library(plyr)
library(SDMTools) 
library(WaveletComp)


#读取数据
setwd("C:/Users/Administrator/Desktop/发表论文/2经济周期协同性-多方法综合评价")


##按指标读取

GDP=read.xlsx(xlsxFile = "GDP水平值.xlsx",sheet = 9,detectDates = TRUE, colNames = TRUE)

Time<<-439



setwd("C:/Users/Administrator/Desktop/发表论文/2经济周期协同性-多方法综合评价/程序")


source('synchronization filter.R',encoding = "UTF-8")

# x=rWishart(n = 1,df = 2,Sigma = diag(2))
# solve(matrix(x,nrow=2,ncol=2))

#先验参数

##转移概率先验参数
u_a_00=8
u_a_01=2
u_a_10=1
u_a_11=9

u_b_00=8
u_b_01=2
u_b_10=1
u_b_11=9

u_sd_00=8
u_sd_01=2
u_sd_10=1
u_sd_11=9


u_S_00=8
u_S_01=2
u_S_10=1
u_S_11=9

u_V_00=8
u_V_01=2
u_V_10=1
u_V_11=9


##均值先验参数
miu_miu=matrix(c(-1,2,-1,2),4,1)
cov_miu=diag(4)

##协方差矩阵先验参数
scale_cov=diag(2)
df_cov=0




#生成初始值
##转移概率矩阵
trans_p_a=matrix(c(0.8,0.2,0.1,0.9),2,2)
trans_p_b=matrix(c(0.8,0.2,0.1,0.9),2,2)
trans_p_sd=matrix(c(0.8,0.2,0.1,0.9),2,2)
trans_p_S=matrix(c(0.8,0.2,0.1,0.9),2,2)
trans_p_V=matrix(c(0.8,0.2,0.1,0.9),2,2)

#均值初始值
miu_vector=matrix(c(-1,2,-1,2),4,1)


#协方差初始值
cov_matrix_list=list(matrix(c(1,0,0,1),2,2),
                     matrix(c(1,0,0,1),2,2))


gibbs_number=6000

burn_in_number=1000


#储存空间
miu_a_0_saver=matrix(0,nrow=5000,ncol=1)
miu_a_1_saver=matrix(0,nrow=5000,ncol=1)
miu_b_0_saver=matrix(0,nrow=5000,ncol=1)
miu_b_1_saver=matrix(0,nrow=5000,ncol=1)

sd_0_saver=rep(list(matrix(0,nrow=2,ncol=2)),5000)
sd_1_saver=rep(list(matrix(0,nrow=2,ncol=2)),5000)

trans_p_a_saver=rep(list(matrix(0,nrow=2,ncol=2)),5000)
trans_p_b_saver=rep(list(matrix(0,nrow=2,ncol=2)),5000)
trans_p_sd_saver=rep(list(matrix(0,nrow=2,ncol=2)),5000)
trans_p_S_saver=rep(list(matrix(0,nrow=2,ncol=2)),5000)
trans_p_V_saver=rep(list(matrix(0,nrow=2,ncol=2)),5000)

state_a_saver=matrix(0,nrow=Time,ncol=5000)
state_b_saver=matrix(0,nrow=Time,ncol=5000)
state_sd_saver=matrix(0,nrow=Time,ncol=5000)
state_S_saver=matrix(0,nrow=Time,ncol=5000)
state_V_saver=matrix(0,nrow=Time,ncol=5000)


#状态变量
state_a=matrix(0,nrow=Time,ncol=1)
state_b=matrix(0,nrow=Time,ncol=1)
state_sd=matrix(0,nrow=Time,ncol=1)
state_S=matrix(0,nrow=Time,ncol=1)
state_V=matrix(0,nrow=Time,ncol=1)


#生成miu的函数4个一起同时生成
miu_vector_estimation=function(y,x,cov_vector,prior_miu,prior_cov)
{
  
  cov_parameter=matrix(0,4,4)
  
  
  for(t in 1:Time)
  {
    cov_parameter=cov_parameter+t(x[[t]])%*%solve(cov_vector[[t]])%*%x[[t]]
  }
  
  cov_miu=solve(solve(prior_cov)+cov_parameter)
  
  
  miu_paramter=matrix(0,4,1)
  
  
  for(t in 1:Time)
  {
    miu_paramter=miu_paramter+t(x[[t]])%*%solve(cov_vector[[t]])%*%matrix(y[,t],nrow=2,ncol=1)
  }
  
  miu_miu=cov_miu%*%(solve(prior_cov)%*%prior_miu+miu_paramter)
  
  factor_miu=rtmvnorm(1,mean=as.vector(miu_miu),sigma = cov_miu,lower = c(-Inf,0,-Inf,0),upper = c(Inf,Inf,Inf,Inf),algorithm = "rejection")
  
  
  return(t(factor_miu))
}


#生成协方差矩阵只生成单一状态的
#miu和state都是两个y变量合并的
cov_matrix_estimation=function(cov_time,state_matrix,y,miu_vector)
{
  df_cov=cov_time+0
  
  scale_cov_parameter=matrix(0,2,2)
  
  for(t in 1:cov_time)
  {
    scale_cov_parameter=scale_cov_parameter+
      (matrix(y[,t],nrow=2,ncol=1)-state_matrix[[t]]%*%miu_vector)%*%
      t(matrix(y[,t],nrow=2,ncol=1)-state_matrix[[t]]%*%miu_vector)
  }
  scale_cov_parameter=scale_cov_parameter+diag(2)
  
  aversed_cov=matrix(rWishart(1,df = df_cov,Sigma = solve(scale_cov_parameter)),2,2)
  
  return(solve(aversed_cov))
}


#生成转移矩阵（每次一个）
trans_p_estimation=function(state_vector)
{
  counter_00=0
  counter_01=0
  counter_10=0
  counter_11=0
  
  for(t in 1:(Time-1))
  {
    if(state_vector[t,]==0&
       state_vector[t+1,]==0)
    {
      counter_00=counter_00+1
    }else if(state_vector[t,]==0&
             state_vector[t+1,]==1)
    {
      counter_01=counter_01+1
    }else if(state_vector[t,]==1&
             state_vector[t+1,]==0)
    {
      counter_10=counter_10+1
    }else if(state_vector[t,]==1&
             state_vector[t+1,]==1)
    {
      counter_11=counter_11+1
    }
  }
  trans_p=matrix(0,2,2)
  trans_p[1,1]=rbeta(n = 1,shape1 =8+counter_00,shape2 =2+counter_01)
  trans_p[2,2]=rbeta(n = 1,shape1 =9+counter_11,shape2 =1+counter_10)
  trans_p[2,1]=1-trans_p[1,1]
  trans_p[1,2]=1-trans_p[2,2]
  
  return(trans_p)
}


for(n in 1:gibbs_number)
{
  

  #生成各个变量的状态变量
  state_list=synchronization_filter(y = t(GDP[,2:3]),
                                    miu_1_vector = t(miu_vector[1:2,]),
                                    miu_2_vector = t(miu_vector[3:4,]),
                                    cov_matrix = cov_matrix_list,
                                    trans_p_a = trans_p_a,
                                    trans_p_b = trans_p_b,
                                    trans_p_sd = trans_p_sd,
                                    trans_p_S = trans_p_S,
                                    trans_p_V = trans_p_V)
  
  state_a=state_list[[1]]
  state_b=state_list[[2]]
  state_sd=state_list[[3]]
  state_S=state_list[[4]]
  state_V=state_list[[5]]


  #估计两个序列在不同区制下的均值
  ##同时生成哑变量方便协方差矩阵的计算
  
  ##生成S_at和S_bt合并的状态变量矩阵
  combine_state_matrix=rep(list(matrix(0,2,4)),Time)
  
  for(t in 1:Time)
  {
    combine_state_matrix[[t]]=matrix(c(1,0,
                                       state_a[t,],0,
                                       0,1,
                                       0,state_b[t,]),nrow = 2,ncol=4)
  }
  
  cov_ts=rep(list(matrix(0,2,2)),Time)
  
  for(t in 1:Time)
  {
    cov_ts[[t]]=cov_matrix_list[[state_sd[t,]+1]]
  }
  
  
  miu_vector=miu_vector_estimation(y = t(GDP[,2:3]),x = combine_state_matrix,
                                 cov_vector = cov_ts,prior_miu = miu_miu,
                                 prior_cov = cov_miu)


  #估计不同区制下的协方差矩阵
  ##保证1状态的协方差矩阵一定大于0状态的（行列式）
  cov_0=list()
  cov_1=list()
  
  
  for(t in 1:Time)
  {
    if(state_sd[t,]==0)
    {
      cov_0=append(cov_0,list(combine_state_matrix[[t]]))
    }else if(state_sd[t,]==1)
    {
      cov_1=append(cov_1,list(combine_state_matrix[[t]]))
    }
  }
  
  cov_matrix_list[[1]]=cov_matrix_estimation(cov_time = length(which(state_sd==0)),
                                             state_matrix = cov_0,
                                             y = t(GDP[which(state_sd==0),2:3]),
                                             miu_vector = miu_vector)
  
  cov_matrix_list[[2]]=cov_matrix_estimation(cov_time = length(which(state_sd==1)),
                                             state_matrix = cov_1,
                                             y = t(GDP[which(state_sd==1),2:3]),
                                             miu_vector = miu_vector)



  #生成转移矩阵
  ##包含先验参数（不叠加）
  
  trans_p_a=trans_p_estimation(state_vector = state_a)
  
  trans_p_b=trans_p_estimation(state_vector = state_b)
  
  trans_p_sd=trans_p_estimation(state_vector = state_sd)
  
  trans_p_S=trans_p_estimation(state_vector = state_S)
  
  trans_p_V=trans_p_estimation(state_vector = state_V)


  #储存每次估计结果
  ##通过将每次估计的状态变量相加，方便后续估计不同时间点上取到
  ##状态的概率
  if(n>burn_in_number)
  {
    miu_a_0_saver[n- burn_in_number,]=miu_vector[1,]
    miu_a_1_saver[n- burn_in_number,]=miu_vector[2,]
    miu_b_0_saver[n- burn_in_number,]=miu_vector[3,]
    miu_b_1_saver[n- burn_in_number,]=miu_vector[4,]
    
  
    sd_0_saver[[n- burn_in_number]]=cov_matrix_list[[1]]
    sd_1_saver[[n- burn_in_number]]=cov_matrix_list[[2]]
    
    trans_p_a_saver[[n- burn_in_number]]=trans_p_a
    trans_p_b_saver[[n- burn_in_number]]=trans_p_b
    trans_p_sd_saver[[n- burn_in_number]]=trans_p_sd
    trans_p_S_saver[[n- burn_in_number]]=trans_p_S
    trans_p_V_saver[[n- burn_in_number]]=trans_p_V
    
    state_a_saver[,n- burn_in_number]=state_a
    state_b_saver[,n- burn_in_number]=state_b
    state_sd_saver[,n- burn_in_number]=state_sd
    state_S_saver[,n- burn_in_number]=state_S
    state_V_saver[,n- burn_in_number]=state_V
    
  }
  
  if(n%%100==0)
  {
    print(n)
  }
}

#参数后验概率

#后验状态概率
state_a_p=rowMeans(state_a_saver)
state_b_p=rowMeans(state_b_saver)
state_sd_p=rowMeans(state_sd_saver)
state_S_p=rowMeans(state_S_saver)
state_V_p=rowMeans(state_V_saver)


#均值
miu_a_0=mean(miu_a_0_saver)
miu_a_1=mean(miu_a_1_saver)
miu_b_0=mean(miu_b_0_saver)
miu_b_1=mean(miu_b_1_saver)

#协方差矩阵
sd_0=matrix(0,2,2)
sd_1=matrix(0,2,2)
trans_p_a=matrix(0,2,2)
trans_p_b=matrix(0,2,2)
trans_p_sd=matrix(0,2,2)
trans_p_S=matrix(0,2,2)
trans_p_V=matrix(0,2,2)



for(i in 1:(gibbs_number- burn_in_number))
{
  sd_0=sd_0+sd_0_saver[[i]]
  sd_1=sd_1+sd_1_saver[[i]]
  trans_p_a=trans_p_a+trans_p_a_saver[[i]]
  trans_p_b=trans_p_b+trans_p_b_saver[[i]]
  trans_p_sd=trans_p_sd+trans_p_sd_saver[[i]]
  trans_p_S=trans_p_S+trans_p_S_saver[[i]]
  trans_p_V=trans_p_V+trans_p_V_saver[[i]]
}


sd_0=sd_0/5000
sd_1=sd_1/5000
trans_p_a=trans_p_a/5000
trans_p_b=trans_p_b/5000
trans_p_sd=trans_p_sd/5000
trans_p_S=trans_p_S/5000
trans_p_V=trans_p_V/5000


plot(x=as.Date(GDP[,1],origin = "1900-01-01"),y=1-state_a_p,type='l')
lines(x=as.Date(GDP[,1],origin = "1900-01-01"),y=1-state_b_p,col='red')


plot(x=as.Date(GDP[,1],origin = "1900-01-01"),y=state_V_p,type='l')
plot(x=as.Date(GDP[,1],origin = "1900-01-01"),y=1-state_sd_p,type='l')
