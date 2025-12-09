fit_RecSBART_VIP<-function(X_train,
                                  X_test,
                                  recurrent_train,
                                  recurrent_test,
                                  terminal_train,
                                  terminal_test,
                                  num_burn=2000,
                                  num_save=500,
                                  num_thin=1){
  
  #original recurrent events
  n<-nrow(X_train)
  p<-ncol(X_train)
  num_recurrent<-sapply(recurrent_train,length)
  N<-sum(num_recurrent)
  
  #scale time in [0,1]
  recurrent_train<-lapply(recurrent_train,function(i){i/max(terminal_train)})
  
  terminal_train<-terminal_train/max(terminal_train)
  
  #estimate the "empirical" intensity function
  safe_max <- function(x) {
    if (length(x) == 0) return(0)
    return(max(x, na.rm = TRUE))
  }
  
  lambda<-mean(num_recurrent)/mean(sapply(recurrent_train,safe_max))
  
  
  shape_lambda_prior<-1
  rate_lambda_prior <- 1/lambda
  
  #hypers for tree structures
  hypers      <- Hypers(cbind(sapply(recurrent_train,safe_max), X_train), 
                        num_recurrent/terminal_train,sigma_hat = 1, num_tree = 50)
  opts        <- Opts(update_sigma = TRUE, update_s = FALSE, 
                      update_alpha = FALSE, update_sigma_mu = FALSE)
  my_forest   <- MakeForest(hypers, opts,warn=FALSE)
  
  #prior distributions
  shape_W_prior<-40
  W <- rgamma(n, shape=shape_W_prior,rate=shape_W_prior)
  
  #log-likelihood function
  loglik_shape_W_prior <- function(shape_W_prior){
    logl <- sum(dgamma(W, shape = shape_W_prior, rate = shape_W_prior, log = T))
    return(logl)
  }
  
  #MCMC iterations
  num_iter <- num_burn+num_save*num_thin
  
  #posterior
  save_idx <- 0
  baseline <- array(NA,num_save)
  shape_baseline_post <- array(NA, num_save)
  rate_baseline_post <- array(NA, num_save)
  
  frailty <- array(NA, dim = c(num_save,n))
  frailty_par <-array(NA, num_save)
  
  #data augmentation
  q <- array(NA,n)
  m <- array(NA,n)
  
  #For computation, assume the recurrence time is after terminal time for 0 recurrent event subject
  empty_indices <- which(num_recurrent == 0)
  
  #estimate the intensity function
  grid <- seq(0.05,1,by=0.05)
  # intensity_estimate <- array(NA, dim = c(num_iter,n,length(grid)))
  
  #estimate the cumulative intensity function
  cumulative_intensity_estimate <- array(NA, dim = c(num_save,n,length(grid)))
  
  #estimate the martingale residual
  sort_time <- sort(unique(unlist(recurrent_train)))
  partition<-c(0,sort_time)
  delta <- diff(partition)
  
  
  #store the variable inclusion proportions
  tree_infor<-matrix(NA,nrow=p+1,ncol=num_save)
  
  start.time = Sys.time()
  for(iter in 1:num_iter){
    z<-NULL
    X_i<-NULL
    X_store<-NULL
    g<-NULL
    time_points <- NULL
    
    for(i in 1:n){
      g<-NULL
      #==================Poisson augmentation=========
      q[i] <- rpois(1, 2*lambda*W[i]*terminal_train[i])
      
      if(q[i] == 0){
        m[i] =0
      }
      else{
        c<-runif(q[i],min = 0, max = 2*lambda*W[i]*terminal_train[i])
        a<-c/(2*lambda*W[i])
        a_x_a<-cbind(a,matrix(rep(X_train[i,],length(a)),ncol = p,nrow = length(a),byrow=T))
        u<-runif(q[i])
        l<-my_forest$do_predict(a_x_a)
        g<-a[which(u<1-pnorm(l))]
        m[i]<- length(g)
      }
      #=======================normal augmentation========
      if(num_recurrent[i]==0 & m[i]==0){
        next
      }
      if(num_recurrent[i]==0 & m[i]!=0){
        z_g<-array(0,m[i])
        for(j in 1:m[i]){
          z_g[j]<-msm::rtnorm(1,mean= my_forest$do_predict(cbind(g[j], t(X_train[i,]))), sd=1, lower=-Inf, upper=0)
        }
        
        z <-c(z,z_g)
        X_g <- matrix(rep(X_train[i,], m[i]), nrow = m[i], ncol = p, byrow = TRUE)
        X_store <- rbind(X_store, X_g)
        time_points <- c(time_points, g)
      }
      if(num_recurrent[i]!=0 & m[i]==0){
        z_t<-array(0,num_recurrent[i])
        for(j in 1:num_recurrent[i]){
          z_t[j]<-msm::rtnorm(1,mean = my_forest$do_predict(cbind(recurrent_train[[i]][j], t(X_train[i,]))), sd=1, lower=0, upper=Inf)
        }
        X_t <- matrix(rep(X_train[i,], num_recurrent[i]), 
                      nrow = num_recurrent[i], ncol = p, byrow = TRUE)
        
        z <- c(z, z_t)
        X_store <- rbind(X_store, X_t)
        time_points <- c(time_points, recurrent_train[[i]])
        
      }
      if(num_recurrent[i]!=0 & m[i]!=0){
        z_t <- array(0,num_recurrent[i])
        for(j in 1:num_recurrent[i]){
          z_t[j]<-msm::rtnorm(1,mean = my_forest$do_predict(cbind(recurrent_train[[i]][j], t(X_train[i,]))), sd=1, lower=0, upper=Inf)
        }
        
        z_g<-array(0,m[i])
        for(j in 1:m[i]){
          z_g[j]<-msm::rtnorm(1,mean= my_forest$do_predict(cbind(g[j], t(X_train[i,]))), sd=1, lower=-Inf, upper=0)
        }
        
        z <- c(z,z_t,z_g)
        X_tg<-matrix(rep(X_train[i,],m[i]+num_recurrent[i]), nrow=m[i]+num_recurrent[i], ncol=p, byrow=T)
        X_store<-rbind(X_store,X_tg)
        time_points <- c(time_points,recurrent_train[[i]],g)
        
      }
      
    }
    #=============================update=====================
    #update baseline
    sum_points<-N+sum(m)
    lambda<-rgamma(1,shape=(shape_lambda_prior+sum_points),rate=(rate_lambda_prior+2*sum(W*terminal_train)))
    
    #update W
    for(i in 1:n){
      W[i]<-rgamma(1,shape=(shape_W_prior+m[i]+num_recurrent[i]),
                   rate=(shape_W_prior+lambda*terminal_train[i]))
    }
    
    
    #update tree
    mu_hat<-my_forest$do_gibbs(cbind(time_points,X_store),z,cbind(time_points,X_store),1)
    
    
    
    #update the eta
    shape_W_prior <-( diversitree::mcmc(lik=loglik_shape_W_prior,
                                        nsteps=1, w=1, x.init=c(shape_W_prior),
                                        prior = function(shape_W_prior){
                                          dgamma(shape_W_prior, shape = 40, rate=0.5, log = TRUE)} 
                                        ,lower=20, upper=Inf) )$pars
    
    
    if (iter > num_burn && ((iter - num_burn - 1) %% num_thin == 0)) {
      save_idx <- (iter - num_burn) / num_thin
    tree_infor[,save_idx]<-my_forest$get_counts()}
      
      
  }
  
 
  
    
  
  
  
  
  
  
  return(tree_infor)
  
}