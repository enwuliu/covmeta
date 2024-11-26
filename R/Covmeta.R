#' overall covariance of random effect meta-analysis in presence of interaction
#' @param b1 the coefficient of main effect from the regression
#' @param v1 the variance of main effect from the regression
#' @param b2 the coefficient of the interaction effect from the regression
#' @param v2 the variance of interaction effect from the regression
#' @param cov_b1b2 the covariance of main effect and interaction effect from the regression
#' @param sample_size the sample size of the study
#' @param effect fixed or random
#' @return overall main effect, overall interaction effect and overall covariance
#' @examples
#' #read in datatable
#' studies<-read.csv('https://raw.githubusercontent.com/enwuliu/meta-analysis/main/random_effect_meta_sim.csv',header=T)
#' #function arguments
#' b1<-studies$b1; #beta1 of the 10 cohorts
#' v1<-studies$var_b1; #variances of the coefficients
#' b2<-studies$b2; #beta1 of the 10 cohorts
#' v2<-studies$var_b2; #variances of the coefficients
#' cov_b1b2<-studies$cov_b1b2; #covariance of b1 and b2
#' sample_size<-studies$sample_size; #sample size
#' cov_meta(b1,v1,b2,v2, cov_b1b2, sample_size,'fixed')
#' cov_meta(b1,v1,b2,v2, cov_b1b2, sample_size,'random')
#' @export

cov_meta<- function(b1=NULL,
                       v1=NULL,
                        b2=NULL,
                        v2=NULL,
                        cov_b1b2=NULL,
                        sample_size=NULL,
                        effect='random',
                        ...){

  if(!is.numeric(b1)){
    stop("The coefficients for the main effect  must be numeric.", call. = FALSE)
  }

  if(!is.numeric(v1)){
    stop("The variance for the main effect must be numeric.", call. = FALSE)
  }


   if(!is.numeric(b2)){
     stop("The coefficients for the interaction effect must be numeric.", call. = FALSE)
   }

   if(!is.numeric(v2)){
     stop("The variance for the interaction effect must be numeric.", call. = FALSE)
   }

  if(!is.numeric(cov_b1b2)){
     stop("The covvariance for the main and interaction effect must be numeric.", call. = FALSE)
  }

 if(!is.numeric(sample_size)){
     stop("The sample size for the studies must be numeric.", call. = FALSE)
 }

  
############functions######################  
  
  fixed_effect_meta<-function(B,V){
    W=1/V
    Beta<-sum(W*B)/sum(W)
    Var=1/sum(W)
    resultlist<-list('overall beta'=Beta,'overall variance'=Var)
    return(resultlist)
  }

  random_effect_meta<-function(B,V){
    W<-1/V                                #Inverse variance weight
    Q<-sum(W*B^2)-(sum(W*B))^2/sum(W)     #total variance or Q statistics
    df<-length(B)-1                       #degree of freedom, number of studies minus 1
    c<-sum(W)-sum(W^2)/sum(W)             #c is the scaling factor

    tau_square<-ifelse(Q>df,(Q-df)/c,0)   #between study variance

    V_star=V+tau_square
    W_star=1/V_star                        #new weight
    Beta_star<-sum(W_star*B)/sum(W_star)   #overall random effect coefficients
    Var_star=1/sum(W_star)                 #ovall variance

    resultlist<-list('overall beta'=Beta_star,'overall variance'=Var_star)
    return(resultlist)

  }

  fixed_effect_meta_r<-function(v1,v2,cov12,sample_size,n_studies){
    r<-cov12/sqrt(v1*v2)                   #correlation coefficient
    z<-0.5*log((1+r)/(1-r))                #Fisher z transformation
    v<-1/(sample_size-3)                   #variance for z
    W<-1/v                                 #weight
    z_overall<-sum(W*z)/sum(W)             #overall random effect of z
    r_overall<-(exp(2*z_overall)-1)/(exp(2*z_overall)+1) #transform back to r
    return(r_overall)
  } 
  
  
  random_effect_meta_r<-function(v1,v2,cov_b1b2,sample_size){
    r<-cov_b1b2/sqrt(v1*v2)           #correlation coefficient
    z<-0.5*log((1+r)/(1-r))          #Fisher z transformation
    V<-1/(sample_size-3)                  #variance for z
    W<-1/V                            #weight
    Q<-sum(W*z^2)-(sum(W*z))^2/sum(W) #total variance or Q statistics
    df<-length(v1)-1                   #degree of freedom, number of studies minus 1
    c<-sum(W)-sum(W^2)/sum(W)         #c is the scaling factor
    
    tau_square<-ifelse(Q>df,(Q-df)/c,0) #between study variance
    
    V_star=V+tau_square
    W_star=1/V_star                     #new weight
    z_star<-sum(W_star*z)/sum(W_star)   #overall random effect of z
    r_overall<-(exp(2*z_star)-1)/(exp(2*z_star)+1) #transform overall z* back to overall correlation coefficient
    return(r_overall)
  }
  
 ####################################################### 
  
  ####fixed effect meta-analysis results
  
  fixed_b1_results<-fixed_effect_meta(b1,v1)
  fixed_b2_results<-fixed_effect_meta(b2,v2)
  fixed_r_results<-fixed_effect_meta_r(v1,v2,cov_b1b2,sample_size)
  fixed_overall_covariance<-fixed_r_results*sqrt(fixed_b1_results$`overall variance`)*sqrt(fixed_b2_results$`overall variance`)
  
  fixed_b1_results$`overall beta`
  
  
  ####random effect meta-analysis results

  random_b1_results<-random_effect_meta(b1,v1)
  random_b2_results<-random_effect_meta(b2,v2)
  random_r_results<-random_effect_meta_r(v1,v2,cov_b1b2,sample_size)
  random_overall_covariance<-random_r_results*sqrt(random_b1_results$`overall variance`)*sqrt(random_b2_results$`overall variance`)
  


if(effect=='fixed'){
   resultlist<-list(
     'fixed overall b1'=  fixed_b1_results$`overall beta`,
     'fixed overall v1'=fixed_b1_results$`overall variance`, 
     'fixed overall b2'=fixed_b2_results$`overall beta`,
     'fixed overall v2'=fixed_b2_results$`overall variance`,
     'fixed overall covariance'= fixed_overall_covariance)
   return(resultlist)
 }
 

if(effect=='random'){
  resultlist<-list(
    'random overall b1'=random_b1_results$`overall beta`,
    'random overall v1'=random_b1_results$`overall variance`, 
    'random overall b2'=random_b2_results$`overall beta`,
    'random overall v2'=random_b2_results$`overall variance`,
    'random overall covariance'= random_overall_covariance)
  return(resultlist)
}

 
}









