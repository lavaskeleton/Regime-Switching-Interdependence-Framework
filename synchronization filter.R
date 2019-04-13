synchronization_filter=function(y,miu_1_vector,miu_2_vector,cov_matrix,
                                trans_p_a,trans_p_b,trans_p_sd,trans_p_S,trans_p_V)
{

  #使用二元而非一元正态分布进行
  
#边际状态概率
p_a=matrix(0,nrow=Time,ncol=2)
p_b=matrix(0,nrow=Time,ncol=2)
p_sd=matrix(0,nrow=Time,ncol=2)
p_S=matrix(0,nrow=Time,ncol=2)
p_V=matrix(0,nrow=Time,ncol=2)




#总状态概率，不包含0期
p_forecast_a=matrix(0,nrow=Time,ncol=2)
p_forecast_b=matrix(0,nrow=Time,ncol=2)
p_forecast_sd=matrix(0,nrow=Time,ncol=2)
p_forecast_S=matrix(0,nrow=Time,ncol=2)
p_forecast_V=matrix(0,nrow=Time,ncol=2)

p_y_sy=matrix(0,nrow=1,ncol=16)
p_ys=matrix(0,nrow=1,ncol = 16)

p_ya=matrix(0,nrow=1,ncol=2)
p_yb=matrix(0,nrow=1,ncol=2)
p_yssd=matrix(0,nrow=1,ncol=2)
p_yS=matrix(0,nrow=1,ncol=2)
p_yV=matrix(0,nrow=1,ncol=2)


p_y=matrix(0,nrow=1,ncol=1)

#初始值
{
  p_forecast_a[1,1]=trans_p_a[1,1]*0.5+trans_p_a[1,2]*0.5
  p_forecast_a[1,2]=trans_p_a[2,1]*0.5+trans_p_a[2,2]*0.5
  
  p_forecast_b[1,1]=trans_p_b[1,1]*0.5+trans_p_b[1,2]*0.5
  p_forecast_b[1,2]=trans_p_b[2,1]*0.5+trans_p_b[2,2]*0.5
  
  p_forecast_sd[1,1]=trans_p_sd[1,1]*0.5+trans_p_sd[1,2]*0.5
  p_forecast_sd[1,2]=trans_p_sd[2,1]*0.5+trans_p_sd[2,2]*0.5
  
  p_forecast_S[1,1]=trans_p_S[1,1]*0.5+trans_p_S[1,2]*0.5
  p_forecast_S[1,2]=trans_p_S[2,1]*0.5+trans_p_S[2,2]*0.5
  
  p_forecast_V[1,1]=trans_p_V[1,1]*0.5+trans_p_V[1,2]*0.5
  p_forecast_V[1,2]=trans_p_V[2,1]*0.5+trans_p_V[2,2]*0.5
  
  
  
  #不同区制下取到y_t的概率p(y_t|S_abstar_t,y_t-1),S_V_t不同取值公式有所不同
  #p_y_sy为过渡变量，表示当期在各种状态和y_t-1下取到y_t的概率
  p_y_sy[,1]=dmvnorm(x=y[,1],
                     mean = c(miu_1_vector[,1],miu_2_vector[,1]),
                     sigma=cov_matrix[[1]])
  
  p_y_sy[,2]=dmvnorm(x=y[,1],
                     mean = c(miu_1_vector[,1],miu_2_vector[,1]+miu_2_vector[,2]),
                     sigma=cov_matrix[[1]])
  
  p_y_sy[,3]=dmvnorm(x=y[,1],
                     mean = c(miu_1_vector[,1]+miu_1_vector[,2],miu_2_vector[,1]),
                     sigma=cov_matrix[[1]])
  
  p_y_sy[,4]=dmvnorm(x=y[,1],
                     mean = c(miu_1_vector[,1]+miu_1_vector[,2],miu_2_vector[,1]+miu_2_vector[,2]),
                     sigma=cov_matrix[[1]])
  
  p_y_sy[,5]=dmvnorm(x=y[,1],
                     mean = c(miu_1_vector[,1],miu_2_vector[,1]),
                     sigma=cov_matrix[[2]])
  
  p_y_sy[,6]=dmvnorm(x=y[,1],
                     mean = c(miu_1_vector[,1],miu_2_vector[,1]+miu_2_vector[,2]),
                     sigma=cov_matrix[[2]])
  
  p_y_sy[,7]=dmvnorm(x=y[,1],
                     mean = c(miu_1_vector[,1]+miu_1_vector[,2],miu_2_vector[,1]),
                     sigma=cov_matrix[[2]])
  
  p_y_sy[,8]=dmvnorm(x=y[,1],
                     mean = c(miu_1_vector[,1]+miu_1_vector[,2],miu_2_vector[,1]+miu_2_vector[,2]),
                     sigma=cov_matrix[[2]])
  
  
  p_y_sy[,9]=dmvnorm(x=y[,1],
                     mean = c(miu_1_vector[,1],miu_2_vector[,1]),
                     sigma=cov_matrix[[1]])
  
  p_y_sy[,10]=0
  
  p_y_sy[,11]=0
  
  p_y_sy[,12]=dmvnorm(x=y[,1],
                      mean = c(miu_1_vector[,1]+miu_1_vector[,2],miu_2_vector[,1]+miu_2_vector[,2]),
                      sigma=cov_matrix[[1]])
  
  p_y_sy[,13]=dmvnorm(x=y[,1],
                      mean = c(miu_1_vector[,1],miu_2_vector[,1]),
                      sigma=cov_matrix[[2]])
  
  p_y_sy[,14]=0
  
  p_y_sy[,15]=0
  
  p_y_sy[,16]=dmvnorm(x=y[,1],
                      mean = c(miu_1_vector[,1]+miu_1_vector[,2],miu_2_vector[,1]+miu_2_vector[,2]),
                      sigma=cov_matrix[[2]])  
  
  
  
  #当期取到y_t与S_abstar_t的联合概率
  #p(y_t,S_abstar_t|y_t-1)=p(y_t|S_abstar_t,y_t-1)*p(S_abstar_t|y_t-1)
  #p_ys为过渡变量，表示y_t-1情况下取到各个状态和y_t的概率
  
  p_ys[,1]=p_y_sy[,1]*p_forecast_a[1,1]*p_forecast_b[1,1]*p_forecast_sd[1,1]*p_forecast_V[1,1]
  p_ys[,2]=p_y_sy[,2]*p_forecast_a[1,1]*p_forecast_b[1,2]*p_forecast_sd[1,1]*p_forecast_V[1,1]
  p_ys[,3]=p_y_sy[,3]*p_forecast_a[1,2]*p_forecast_b[1,1]*p_forecast_sd[1,1]*p_forecast_V[1,1]
  p_ys[,4]=p_y_sy[,4]*p_forecast_a[1,2]*p_forecast_b[1,2]*p_forecast_sd[1,1]*p_forecast_V[1,1]
  p_ys[,5]=p_y_sy[,5]*p_forecast_a[1,1]*p_forecast_b[1,1]*p_forecast_sd[1,2]*p_forecast_V[1,1]
  p_ys[,6]=p_y_sy[,6]*p_forecast_a[1,1]*p_forecast_b[1,2]*p_forecast_sd[1,2]*p_forecast_V[1,1]
  p_ys[,7]=p_y_sy[,7]*p_forecast_a[1,2]*p_forecast_b[1,1]*p_forecast_sd[1,2]*p_forecast_V[1,1]
  p_ys[,8]=p_y_sy[,8]*p_forecast_a[1,2]*p_forecast_b[1,2]*p_forecast_sd[1,2]*p_forecast_V[1,1]
  
  p_ys[,9]=p_y_sy[,9]*p_forecast_S[1,1]*p_forecast_sd[1,1]*p_forecast_V[1,2]
  p_ys[,10]=0
  p_ys[,11]=0
  p_ys[,12]=p_y_sy[,12]*p_forecast_S[1,2]*p_forecast_sd[1,1]*p_forecast_V[1,2]
  p_ys[,13]=p_y_sy[,13]*p_forecast_S[1,1]*p_forecast_sd[1,2]*p_forecast_V[1,2]
  p_ys[,14]=0
  p_ys[,15]=0
  p_ys[,16]=p_y_sy[,16]*p_forecast_S[1,2]*p_forecast_sd[1,2]*p_forecast_V[1,2]
  
  #各个状态变量边际概率密度的p(y_t,S_*_t|y_t-1)
  #p(y_t,S_a_t|y_t-1)=sum(p(y_t,S_abstar_t|y_t-1)) for all S_b_t,S_sd_t,S_V_t
  #p(y_t,S_b_t|y_t-1)=sum(p(y_t,S_abstar_t|y_t-1)) for all S_a_t,S_sd_t,S_V_t
  #p(y_t,S_sd_t|y_t-1)=sum(p(y_t,S_abstar_t|y_t-1)) for all S_a_t,S_b_t,S_V_t
  #p(y_t,S_V_t|y_t-1)=sum(p(y_t,S_abstar_t|y_t-1)) for all S_a_t,S_b_t,S_sd_t
  
  #p(y_t,S_t|y_t-1)=sum(p(y_t,S_ab_t,S_V_t==1|y_t-1)) for all S_a_t,S_b_t,S_sd_t
  
  p_ya[,1]=p_ys[,1]+p_ys[,2]+p_ys[,5]+p_ys[,6]+p_ys[,9]+p_ys[,10]+p_ys[,13]+p_ys[,14]
  p_ya[,2]=p_ys[,3]+p_ys[,4]+p_ys[,7]+p_ys[,8]+p_ys[,11]+p_ys[,12]+p_ys[,15]+p_ys[,16]
  
  p_yb[,1]=p_ys[,1]+p_ys[,3]+p_ys[,5]+p_ys[,7]+p_ys[,9]+p_ys[,11]+p_ys[,13]+p_ys[,15]
  p_yb[,2]=p_ys[,2]+p_ys[,4]+p_ys[,6]+p_ys[,8]+p_ys[,10]+p_ys[,12]+p_ys[,14]+p_ys[,16]
  
  
  p_yssd[,1]=p_ys[,1]+p_ys[,2]+p_ys[,3]+p_ys[,4]+p_ys[,9]+p_ys[,10]+p_ys[,11]+p_ys[,12]
  p_yssd[,2]=p_ys[,5]+p_ys[,6]+p_ys[,7]+p_ys[,8]+p_ys[,13]+p_ys[,14]+p_ys[,15]+p_ys[,16]
  
  p_yS[,1]=p_ys[,1]+p_ys[,5]+p_ys[,9]+p_ys[,13]
  p_yS[,2]=p_ys[,4]+p_ys[,8]+p_ys[,12]+p_ys[,16]
  
  p_yV[,1]=p_ys[,1]+p_ys[,2]+p_ys[,3]+p_ys[,4]+p_ys[,5]+p_ys[,6]+p_ys[,7]+p_ys[,8]
  p_yV[,2]=p_ys[,9]+p_ys[,10]+p_ys[,11]+p_ys[,12]+p_ys[,13]+p_ys[,14]+p_ys[,15]+p_ys[,16]
  
  
  
  #无条件概率密度
  #p(y_t|y_t-1)=sum(p(y_t,S_abstar_t|y_t-1)) for all S_a_t,S_b_t,S_sd_t,S_V_t
  #p(y_t|y_t-1)_star=sum(p(y_t,S_t|y_t-1)) for all S_t
  
  p_y=rowSums(p_ys)
  
  
  
  
  #更新阶段p(S_*_t|y_t)
  #p(S_a_t|y_t)=p(y_t,S_a_t|y_t-1)/p(y_t|y_t-1)
  #p(S_b_t|y_t)=p(y_t,S_b_t|y_t-1)/p(y_t|y_t-1)
  #p(S_sd_t|y_t)=p(y_t,S_sd_t|y_t-1)/p(y_t|y_t-1)
  #p(S_V_t|y_t)=p(y_t,S_V_t|y_t-1)/p(y_t|y_t-1)
  #p(S_t|y_t)=p(y_t,S_t|y_t-1)/p(y_t|y_t-1)_star
  
  p_a[1,1]=p_ya[,1]/p_y
  p_a[1,2]=p_ya[,2]/p_y
  
  p_b[1,1]=p_yb[,1]/p_y
  p_b[1,2]=p_yb[,2]/p_y
  
  p_sd[1,1]=p_yssd[,1]/p_y
  p_sd[1,2]=p_yssd[,2]/p_y
  
  
  p_S[1,1]=p_yS[,1]/rowSums(p_yS)
  p_S[1,2]=p_yS[,2]/rowSums(p_yS)
  
  p_V[1,1]=p_yV[,1]/p_y
  p_V[1,2]=p_yV[,2]/p_y
}









for(t in 1:(Time-1))
{
  #t=1
  #预测阶段
  #预测当期区制概率阶段p(S_t|y_t-1)=p(S_t|S_t-1)*p(S_t-1|y_t-1)
  #前一期的后验状态概率乘以转移概率作为当期的先验状态预测概率
  #p(S_a_t|y_t-1)=sum(p(S_a_t|S_a_t-1)*p(S_a_t-1|y_t-1)) for all S_a_t-1
  #p(S_b_t|y_t-1)=sum(p(S_b_t|S_b_t-1)*p(S_b_t-1|y_t-1)) for all S_b_t-1
  #p(S_sd_t|y_t-1)=sum(p(S_sd_t|S_sd_t-1)*p(S_sd_t-1|y_t-1)) for all S_sd_t-1
  #p(S_V_t|y_t-1)=sum(p(S_V_t|S_V_t-1)*p(S_V_t-1|y_t-1)) for all S_V_t-1
  #p(S_t|y_t-1)=sum(p(S_t|S_t-1)*p(S_t-1|y_t-1)) for all S_t-1
  
  #p_a,p_b,p_sd,p_V,p_S代表当期各个状态变量处于某个状态的概率
  
  #总状态先验概率有两部分，S_V_t==0和S_V_t==1的公式有所不同
  #p(S_abstar_t|y_t-1)=p(S_V_t==0|y_t-1)*(p(S_a_t|y_t-1)*p(S_b_t|y_t-1)*p(S_sd_t|y_t-1))
  #或=p(S_V_t==1|y_t-1)*(p(S_t)*p(S_sd_t|y_t-1))
  p_forecast_a[t+1,1]=trans_p_a[1,1]*p_a[t,1]+trans_p_a[1,2]*p_a[t,2]
  p_forecast_a[t+1,2]=trans_p_a[2,1]*p_a[t,1]+trans_p_a[2,2]*p_a[t,2]
  
  p_forecast_b[t+1,1]=trans_p_b[1,1]*p_b[t,1]+trans_p_b[1,2]*p_b[t,2]
  p_forecast_b[t+1,2]=trans_p_b[2,1]*p_b[t,1]+trans_p_b[2,2]*p_b[t,2]
  
  p_forecast_sd[t+1,1]=trans_p_sd[1,1]*p_sd[t,1]+trans_p_sd[1,2]*p_sd[t,2]
  p_forecast_sd[t+1,2]=trans_p_sd[2,1]*p_sd[t,1]+trans_p_sd[2,2]*p_sd[t,2]
  
  p_forecast_S[t+1,1]=trans_p_S[1,1]*p_S[t,1]+trans_p_S[1,2]*p_S[t,2]
  p_forecast_S[t+1,2]=trans_p_S[2,1]*p_S[t,1]+trans_p_S[2,2]*p_S[t,2]
  
  p_forecast_V[t+1,1]=trans_p_V[1,1]*p_V[t,1]+trans_p_V[1,2]*p_V[t,2]
  p_forecast_V[t+1,2]=trans_p_V[2,1]*p_V[t,1]+trans_p_V[2,2]*p_V[t,2]
  
  
  
  #不同区制下取到y_t的概率p(y_t|S_abstar_t,y_t-1),S_V_t不同取值公式有所不同
  #p_y_sy为过渡变量，表示当期在各种状态和y_t-1下取到y_t的概率
  p_y_sy[,1]=dmvnorm(x=y[,t+1],
                     mean = c(miu_1_vector[,1],miu_2_vector[,1]),
                     sigma= cov_matrix[[1]])
  
  p_y_sy[,2]=dmvnorm(x=y[,t+1],
                     mean = c(miu_1_vector[,1],miu_2_vector[,1]+miu_2_vector[,2]),
                     sigma= cov_matrix[[1]])
  
  p_y_sy[,3]=dmvnorm(x=y[,t+1],
                     mean = c(miu_1_vector[,1]+miu_1_vector[,2],miu_2_vector[,1]),
                     sigma= cov_matrix[[1]])
  
  p_y_sy[,4]=dmvnorm(x=y[,t+1],
                     mean = c(miu_1_vector[,1]+miu_1_vector[,2],miu_2_vector[,1]+miu_2_vector[,2]),
                     sigma= cov_matrix[[1]])
  
  p_y_sy[,5]=dmvnorm(x=y[,t+1],
                     mean = c(miu_1_vector[,1],miu_2_vector[,1]),
                     sigma= cov_matrix[[2]])
  
  p_y_sy[,6]=dmvnorm(x=y[,t+1],
                     mean = c(miu_1_vector[,1],miu_2_vector[,1]+miu_2_vector[,2]),
                     sigma= cov_matrix[[2]])
  
  p_y_sy[,7]=dmvnorm(x=y[,t+1],
                     mean = c(miu_1_vector[,1]+miu_1_vector[,2],miu_2_vector[,1]),
                     sigma= cov_matrix[[2]])
  
  p_y_sy[,8]=dmvnorm(x=y[,t+1],
                     mean = c(miu_1_vector[,1]+miu_1_vector[,2],miu_2_vector[,1]+miu_2_vector[,2]),
                     sigma= cov_matrix[[2]])
  
  
  p_y_sy[,9]=dmvnorm(x=y[,t+1],
                     mean = c(miu_1_vector[,1],miu_2_vector[,1]),
                     sigma= cov_matrix[[1]])
  
  p_y_sy[,10]=0
  
  p_y_sy[,11]=0
  
  p_y_sy[,12]=dmvnorm(x=y[,t+1],
                      mean = c(miu_1_vector[,1]+miu_1_vector[,2],miu_2_vector[,1]+miu_2_vector[,2]),
                      sigma= cov_matrix[[1]])
  
  p_y_sy[,13]=dmvnorm(x=y[,t+1],
                      mean = c(miu_1_vector[,1],miu_2_vector[,1]),
                      sigma= cov_matrix[[2]])
  
  p_y_sy[,14]=0
  
  p_y_sy[,15]=0
  
  p_y_sy[,16]=dmvnorm(x=y[,t+1],
                      mean = c(miu_1_vector[,1]+miu_1_vector[,2],miu_2_vector[,1]+miu_2_vector[,2]),
                      sigma= cov_matrix[[2]])  
  
  
  
  #当期取到y_t与S_abstar_t的联合概率
  #p(y_t,S_abstar_t|y_t-1)=p(y_t|S_abstar_t,y_t-1)*p(S_abstar_t|y_t-1)
  #p_ys为过渡变量，表示y_t-1情况下取到各个状态和y_t的概率
  
  p_ys[,1]=p_y_sy[,1]*p_forecast_a[t+1,1]*p_forecast_b[t+1,1]*p_forecast_sd[t+1,1]*p_forecast_V[t+1,1]
  p_ys[,2]=p_y_sy[,2]*p_forecast_a[t+1,1]*p_forecast_b[t+1,2]*p_forecast_sd[t+1,1]*p_forecast_V[t+1,1]
  p_ys[,3]=p_y_sy[,3]*p_forecast_a[t+1,2]*p_forecast_b[t+1,1]*p_forecast_sd[t+1,1]*p_forecast_V[t+1,1]
  p_ys[,4]=p_y_sy[,4]*p_forecast_a[t+1,2]*p_forecast_b[t+1,2]*p_forecast_sd[t+1,1]*p_forecast_V[t+1,1]
  p_ys[,5]=p_y_sy[,5]*p_forecast_a[t+1,1]*p_forecast_b[t+1,1]*p_forecast_sd[t+1,2]*p_forecast_V[t+1,1]
  p_ys[,6]=p_y_sy[,6]*p_forecast_a[t+1,1]*p_forecast_b[t+1,2]*p_forecast_sd[t+1,2]*p_forecast_V[t+1,1]
  p_ys[,7]=p_y_sy[,7]*p_forecast_a[t+1,2]*p_forecast_b[t+1,1]*p_forecast_sd[t+1,2]*p_forecast_V[t+1,1]
  p_ys[,8]=p_y_sy[,8]*p_forecast_a[t+1,2]*p_forecast_b[t+1,2]*p_forecast_sd[t+1,2]*p_forecast_V[t+1,1]
  
  p_ys[,9]=p_y_sy[,9]*p_forecast_S[t+1,1]*p_forecast_sd[t+1,1]*p_forecast_V[t+1,2]
  p_ys[,10]=0
  p_ys[,11]=0
  p_ys[,12]=p_y_sy[,12]*p_forecast_S[t+1,2]*p_forecast_sd[t+1,1]*p_forecast_V[t+1,2]
  p_ys[,13]=p_y_sy[,13]*p_forecast_S[t+1,1]*p_forecast_sd[t+1,2]*p_forecast_V[t+1,2]
  p_ys[,14]=0
  p_ys[,15]=0
  p_ys[,16]=p_y_sy[,16]*p_forecast_S[t+1,2]*p_forecast_sd[t+1,2]*p_forecast_V[t+1,2]
  
  #各个状态变量边际概率密度的p(y_t,S_*_t|y_t-1)
  #p(y_t,S_a_t|y_t-1)=sum(p(y_t,S_abstar_t|y_t-1)) for all S_b_t,S_sd_t,S_V_t
  #p(y_t,S_b_t|y_t-1)=sum(p(y_t,S_abstar_t|y_t-1)) for all S_a_t,S_sd_t,S_V_t
  #p(y_t,S_sd_t|y_t-1)=sum(p(y_t,S_abstar_t|y_t-1)) for all S_a_t,S_b_t,S_V_t
  #p(y_t,S_V_t|y_t-1)=sum(p(y_t,S_abstar_t|y_t-1)) for all S_a_t,S_b_t,S_sd_t
  
  #p(y_t,S_t|y_t-1)=sum(p(y_t,S_ab_t,S_V_t==1|y_t-1)) for all S_a_t,S_b_t,S_sd_t
  
  p_ya[,1]=p_ys[,1]+p_ys[,2]+p_ys[,5]+p_ys[,6]+p_ys[,9]+p_ys[,10]+p_ys[,13]+p_ys[,14]
  p_ya[,2]=p_ys[,3]+p_ys[,4]+p_ys[,7]+p_ys[,8]+p_ys[,11]+p_ys[,12]+p_ys[,15]+p_ys[,16]
  
  p_yb[,1]=p_ys[,1]+p_ys[,3]+p_ys[,5]+p_ys[,7]+p_ys[,9]+p_ys[,11]+p_ys[,13]+p_ys[,15]
  p_yb[,2]=p_ys[,2]+p_ys[,4]+p_ys[,6]+p_ys[,8]+p_ys[,10]+p_ys[,12]+p_ys[,14]+p_ys[,16]
  
  
  p_yssd[,1]=p_ys[,1]+p_ys[,2]+p_ys[,3]+p_ys[,4]+p_ys[,9]+p_ys[,10]+p_ys[,11]+p_ys[,12]
  p_yssd[,2]=p_ys[,5]+p_ys[,6]+p_ys[,7]+p_ys[,8]+p_ys[,13]+p_ys[,14]+p_ys[,15]+p_ys[,16]
  
  p_yS[,1]=p_ys[,1]+p_ys[,5]+p_ys[,9]+p_ys[,13]
  p_yS[,2]=p_ys[,4]+p_ys[,8]+p_ys[,12]+p_ys[,16]
  
  p_yV[,1]=p_ys[,1]+p_ys[,2]+p_ys[,3]+p_ys[,4]+p_ys[,5]+p_ys[,6]+p_ys[,7]+p_ys[,8]
  p_yV[,2]=p_ys[,9]+p_ys[,10]+p_ys[,11]+p_ys[,12]+p_ys[,13]+p_ys[,14]+p_ys[,15]+p_ys[,16]
  
  
  
  #无条件概率密度
  #p(y_t|y_t-1)=sum(p(y_t,S_abstar_t|y_t-1)) for all S_a_t,S_b_t,S_sd_t,S_V_t
  #p(y_t|y_t-1)_star=sum(p(y_t,S_t|y_t-1)) for all S_t
  
  p_y=rowSums(p_ys)
  
  
  
  
  #更新阶段p(S_*_t|y_t)
  #p(S_a_t|y_t)=p(y_t,S_a_t|y_t-1)/p(y_t|y_t-1)
  #p(S_b_t|y_t)=p(y_t,S_b_t|y_t-1)/p(y_t|y_t-1)
  #p(S_sd_t|y_t)=p(y_t,S_sd_t|y_t-1)/p(y_t|y_t-1)
  #p(S_V_t|y_t)=p(y_t,S_V_t|y_t-1)/p(y_t|y_t-1)
  #p(S_t|y_t)=p(y_t,S_t|y_t-1)/p(y_t|y_t-1)_star
  
  p_a[t+1,1]=p_ya[,1]/p_y
  p_a[t+1,2]=p_ya[,2]/p_y
  
  p_b[t+1,1]=p_yb[,1]/p_y
  p_b[t+1,2]=p_yb[,2]/p_y
  
  p_sd[t+1,1]=p_yssd[,1]/p_y
  p_sd[t+1,2]=p_yssd[,2]/p_y
  
  
  p_S[t+1,1]=p_yS[,1]/rowSums(p_yS)
  p_S[t+1,2]=p_yS[,2]/rowSums(p_yS)
  
  p_V[t+1,1]=p_yV[,1]/p_y
  p_V[t+1,2]=p_yV[,2]/p_y
  
  #t=t+1
}

#平滑
# #状态变量
state_a=matrix(0,nrow=Time,ncol=1)
state_b=matrix(0,nrow=Time,ncol=1)
state_sd=matrix(0,nrow=Time,ncol=1)
state_S=matrix(0,nrow=Time,ncol=1)
state_V=matrix(0,nrow=Time,ncol=1)


# ###生成最后一期状态变量
state_a[Time,]=sample(x=c(0,1),size = 1,prob = p_a[Time,])
state_b[Time,]=sample(x=c(0,1),size = 1,prob = p_b[Time,])
state_sd[Time,]=sample(x=c(0,1),size = 1,prob = p_sd[Time,])
state_S[Time,]=sample(x=c(0,1),size = 1,prob = p_S[Time,])
state_V[Time,]=sample(x=c(0,1),size = 1,prob = p_V[Time,])


# ##反向更新生成状态变量
# #p(s_t|y_t)=p(S_t+1|S_t)*p(S_t|y_t)/p(S_t+1|y_t)
# ###最终状态变量概率
p_a_final=matrix(0,nrow=Time,ncol=2)
p_b_final=matrix(0,nrow=Time,ncol=2)
p_sd_final=matrix(0,nrow=Time,ncol=2)
p_S_final=matrix(0,nrow=Time,ncol=2)
p_V_final=matrix(0,nrow=Time,ncol=2)


  for(t in (Time-1):1)
  {

    p_a_final[t,1]=trans_p_a[state_a[t+1,]+1,1]*p_a[t,1]/p_forecast_a[t+1,state_a[t+1,]+1]
    p_a_final[t,2]=trans_p_a[state_a[t+1,]+1,2]*p_a[t,2]/p_forecast_a[t+1,state_a[t+1,]+1]

    p_b_final[t,1]=trans_p_b[state_b[t+1,]+1,1]*p_b[t,1]/p_forecast_b[t+1,state_b[t+1,]+1]
    p_b_final[t,2]=trans_p_b[state_b[t+1,]+1,2]*p_b[t,2]/p_forecast_b[t+1,state_b[t+1,]+1]


    p_sd_final[t,1]=trans_p_sd[state_sd[t+1,]+1,1]*p_sd[t,1]/p_forecast_sd[t+1,state_sd[t+1,]+1]
    p_sd_final[t,2]=trans_p_sd[state_sd[t+1,]+1,2]*p_sd[t,2]/p_forecast_sd[t+1,state_sd[t+1,]+1]


    p_S_final[t,1]=trans_p_S[state_S[t+1,]+1,1]*p_S[t,1]/p_forecast_S[t+1,state_S[t+1,]+1]
    p_S_final[t,2]=trans_p_S[state_S[t+1,]+1,2]*p_S[t,2]/p_forecast_S[t+1,state_S[t+1,]+1]


    p_V_final[t,1]=trans_p_V[state_V[t+1,]+1,1]*p_V[t,1]/p_forecast_V[t+1,state_V[t+1,]+1]
    p_V_final[t,2]=trans_p_V[state_V[t+1,]+1,2]*p_V[t,2]/p_forecast_V[t+1,state_V[t+1,]+1]


    state_a[t,]=sample(x=c(0,1),size = 1,prob = p_a_final[t,])
    state_b[t,]=sample(x=c(0,1),size = 1,prob = p_b_final[t,])
    state_sd[t,]=sample(x=c(0,1),size = 1,prob = p_sd_final[t,])
    state_S[t,]=sample(x=c(0,1),size = 1,prob = p_S_final[t,])
    state_V[t,]=sample(x=c(0,1),size = 1,prob = p_V_final[t,])

  }
  return(list(state_a,state_b,state_sd,state_S,state_V))
}

