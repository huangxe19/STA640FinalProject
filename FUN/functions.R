# Generate synthetic data
datagen <- function(N, k, theta, sparse, random_d, pd) {
  N = N
  k = k
  b = 1 / (1:k)
  var=1
  
  # Generate covariance matrix of X
  sigma <- genPositiveDefMat(k, "unifcorrmat")$Sigma
  sigma <- cov2cor(sigma)
  
  # Generate covariates X
  X <- rmvnorm(N, sigma = sigma)
  
  if (sparse==T) {
    for (i in 1:k){
      sparse_id <- sample_n(data.frame(id=1:N), 0.6*N)[,1] %>% as.vector()
      X[sparse_id, i] <- 0
    }
  }
  
  # Options for D (m_(X))
  if (random_d == T) {
    d <- rbinom(n=N, 1, p=pd)
  } else {
    d_prop <- pnorm(X %*% b) # D is dependent on Za
    d <- as.numeric(rbinom(N, prob = d_prop, size = 1))
  }
  
  # Options for theta
  if (theta == "het") {
    theta_s <- 1 + 2*X[,2]*X[,5] + X[,10]
    theta <-
      (theta_s - min(theta_s)) * (1 - (-1)) / (max(theta_s) - min(theta_s)) - 1
  } else {
    theta == theta
  }
  
  g <- as.vector(cos(X %*% b) ^ 2)
  
  y <- theta * d + g + rnorm(N,0,var)
  
  data <- as.data.frame(y)
  data <- cbind(data, theta, d, X)
  colnames(data) <- c("Y", "theta", "Z", c(paste0("V", 1:k)))
  
  return(data)
}


# DML to estimate ATE
my_dml_cfit <- function(DF, folds){
  res_dml <- rep(NA, 2)
  DF <- DF %>% select(-ID) %>% mutate(Z = as.factor(Z))
  for(i in 1:2){
    # sample splitting for honesty criterion
    samp1 = DF[-folds[[i]],]
    samp2 = DF[folds[[i]],]
    
    # tuning mtry for random forest model E(Y|X)
    tuneGrid <- data.frame(.mtry = c(1,2))
    
    rf_u <- train(Y~.,
                  data = samp1 %>% dplyr::select(-Z),
                  method = 'rf',
                  tuneGrid = tuneGrid,
                  trControl = trainControl(method = "cv", 
                                           number = 5, 
                                           verboseIter = FALSE))
    # predict on the untouched half
    pred_u <- predict(rf_u, 
                      newdata= samp2 %>% select(-Z)) %>% as.numeric()
    
    # get the residual
    u_res <- ifelse(samp2$Y==TRUE, 1, 0) - pred_u
    
    # tuning mtry for random forest model E(Z|X)
    rf_v <- train(Z~.,
                  data = samp1 %>% dplyr::select(-Y),
                  method = 'rf',
                  tuneGrid = tuneGrid,
                  trControl = trainControl(method = "cv", 
                                           number = 5, 
                                           verboseIter = FALSE))
    
    # predict on the untouched half
    pred_v <- predict(rf_v, 
                      newdata= samp2 %>% dplyr::select(-Y),
                      type = "prob")[, 2]
    
    # get the residual
    v_res <- ifelse(samp2$Z==TRUE, 1, 0) - pred_v
    
    lm_fit <- lm(u_res~v_res, data = data.frame(cbind(u_res, v_res)))
    res_dml[i] <- as.numeric(lm_fit$coefficients[2])
  }
  return(mean(res_dml))
}


# DML to estimate ATE without kfold
my_dml_naive <- function(DF){
    DF <- DF %>% select(-ID) %>% mutate(Z = as.factor(Z))
    # tuning mtry for random forest model E(Y|X)
    tuneGrid <- data.frame(.mtry = c(1,2))
    
    rf_u <- train(Y~.,
                  data = DF %>% dplyr::select(-Z),
                  method = 'rf',
                  tuneGrid = tuneGrid,
                  trControl = trainControl(method = "cv", 
                                           number = 5, 
                                           verboseIter = FALSE))
    
    # predict on the untouched half
    pred_u <- predict(rf_u, 
                      newdata= DF %>% select(-Z)) %>% as.numeric()
    
    # get the residual
    u_res <- ifelse(DF$Y==TRUE, 1, 0) - pred_u
    
    # tuning mtry for random forest model E(Z|X)
    rf_v <- train(Z~.,
                  data = DF %>% dplyr::select(-Y),
                  method = 'rf',
                  tuneGrid = tuneGrid,
                  trControl = trainControl(method = "cv", 
                                           number = 5, 
                                           verboseIter = FALSE))
    
  # predict on the untouched half
  pred_v <- predict(rf_v, 
                    newdata= DF %>% dplyr::select(-Y),
                    type = "prob")[, 2]
    
  # get the residual
  v_res <- ifelse(DF$Z==TRUE, 1, 0) - pred_v
    
  lm_fit <- lm(u_res~v_res, data = data.frame(cbind(u_res, v_res)))
  res_dml <- as.numeric(lm_fit$coefficients[2])
    
  return(res_dml)
}


# BART naive
bart_naive <- function(DF) {
  bart_out <- bartc(confounders = DF[, 4:ncol(DF)],
                     treatment = DF$Z,
                     response = DF$Y,
                     method.rsp='bart',
                     method.trt='none',
                     estimand='ate',
                     p.scoreAsCovariate = F,
                     keepTrees=T,
                    verbose=F)
  bart_post <- predict(bart_out,
                       newdata = DF,
                       type = 'icate')
  ite_bart <- data.frame(x = colMeans(bart_post))
  return(ite_bart)
}


# BART sample-splitting
bart_cfit <- function(DF, M=15) {
  DF <- DF %>% mutate(row_id = c(1:nrow(DF)))
  res_mat <- data.frame(row_id = c(1:nrow(DF)), ITE = rep(NA, nrow(DF)))
  
  if(nrow(DF) == 2000){
    nsub = 1900
  }
  else{
    nsub=450
  }
  for (i in 1:M) {
    # select nsub IDs from original nrow IDs
    samp_curr_id <- sample(c(1:nrow(DF)), size = nsub, replace = F)
    
    folds1 <- sample(samp_curr_id, size = length(samp_curr_id)/2, replace = F)
    folds2 <- setdiff(samp_curr_id, folds1)
    
    samp1 <- DF %>% filter(row_id %in% folds1)
    samp2 <- DF %>% filter(row_id %in% folds2)
    
    bart_out <- bartc(confounders = samp1[, 4:(ncol(samp1)-1)],
                      treatment = samp1$Z,
                      response = samp1$Y,
                      method.rsp ='bart',
                      method.trt ='none',
                      estimand ='ate',
                      p.scoreAsCovariate = F,
                      keepTrees = T,
                      verbose=F)
    
    bart_post <- predict(bart_out,
                         newdata = samp2,
                         type = 'icate')
    
    pred_curr <- data.frame(row_id = folds2, 
                            res = colMeans(bart_post))
    
    new_res <- merge(x = res_mat, 
                     y = pred_curr, 
                     by='row_id', 
                     all.x=T)
    meanpred <- rowMeans(new_res[,2:3], na.rm = T)
    res_mat$ITE <- meanpred
  }
  return(res_mat)
}


# BART-ps naive
bartps_naive <- function(DF) {
  bartps_out <- bartc(confounders = DF[, 4:ncol(DF)],
                      treatment = DF$Z,
                      response = DF$Y,
                      method.rsp='bart',
                      method.trt='glm',
                      estimand='ate',
                      keepTrees=T,
                      verbose=F)
  bartps_post <- predict(bartps_out,
                       newdata = DF,
                       type = 'icate')
  ite_bartps <- data.frame(x = colMeans(bartps_post))
  return(ite_bartps)
}


# BART-ps sample-splitting
bartps_cfit <- function(DF, M=15) {
  DF <- DF %>% mutate(row_id = c(1:nrow(DF)))
  res_mat <- data.frame(row_id = c(1:nrow(DF)), ITE = rep(NA, nrow(DF)))
  
  if(nrow(DF) == 2000){
    nsub = 1920
  }
  else{
    nsub=450
  }
  for (i in 1:M) {
    # select nsub IDs from original nrow IDs
    samp_curr_id <- sample(c(1:nrow(DF)), size = nsub, replace = F)

    folds1 <- sample(samp_curr_id, size = length(samp_curr_id)/3, replace = F)
    folds_remain <- setdiff(samp_curr_id, folds1)
    folds2 <- sample(folds_remain, size = length(folds_remain)/2, replace = F)
    folds3 <- setdiff(folds_remain, folds2)
    
    samp1 <- DF %>% filter(row_id %in% folds1)
    samp2 <- DF %>% filter(row_id %in% folds2)
    samp3 <- DF %>% filter(row_id %in% folds3)
    
    # train ps model on sample 1
    ps1_mod <- glm(Z~.-theta-Y-row_id, data = samp1, family='binomial')
    
    # predict ps on sample 2
    ps1 <- predict(ps1_mod, newdata = samp2, type='response')
    samp2$ps_est <- ps1
    ps2 <- predict(ps1_mod, newdata = samp3, type='response')
    samp3$ps_est <- ps2
    
    # train outcome on sample 2
    bart_out <- bartc(confounders = samp2[, 4:(ncol(samp2)-1)],
                      treatment = samp2$Z,
                      response = samp2$Y,
                      method.rsp ='bart',
                      method.trt ='none',
                      estimand ='ate',
                      p.scoreAsCovariate = F,
                      keepTrees = T,
                      verbose=F)
    
    bart_post <- predict(bart_out,
                         newdata = samp3,
                         type = 'icate')
    
    pred_curr <- data.frame(row_id = folds3, 
                            res = colMeans(bart_post))
    
    new_res <- merge(x = res_mat, 
                     y = pred_curr, 
                     by='row_id', 
                     all.x=T)
    meanpred <- rowMeans(new_res[,2:3], na.rm = T)
    res_mat$ITE <- meanpred
  }
  return(res_mat)
}


# BCF naive
bcf_naive <- function(DF) {
  ps1_mod <- glm(Z~.-theta-Y, data=DF, family='binomial')
  ps1 <- predict(ps1_mod, newdata=DF, type='response')
  
  x_design <- makeModelMatrixFromDataFrame(DF[, 4:ncol(DF)])
  
  bcf_out <- bcf(y = DF$Y, 
                 z = DF$Z,
                 x_control = x_design, 
                 x_moderate = x_design,
                 pihat = ps1, 
                 nburn =1000, 
                 nsim =3000)
  ite_bcf <- colMeans(bcf_out$tau)
  return(ite_bcf)
}


# BCF sample-splitting
bcf_cfit <- function(DF, M=15) {
  DF <- DF %>% mutate(row_id = c(1:nrow(DF)))
  res_mat <- data.frame(row_id = c(1:nrow(DF)), ITE = rep(NA, nrow(DF)))
  
  if(nrow(DF) == 2000){
    nsub = 1920
  }
  else{
    nsub=450
  }
  for (i in 1:M) {
    # select nsub IDs from original nrow IDs
    samp_curr_id <- sample(c(1:nrow(DF)), size = nsub, replace = F)

    folds1 <- sample(samp_curr_id, size = length(samp_curr_id)/3, replace = F)
    folds_remain <- setdiff(samp_curr_id, folds1)
    folds2 <- sample(folds_remain, size = length(folds_remain)/2, replace = F)
    folds3 <- setdiff(folds_remain, folds2)
    
    samp1 <- DF %>% filter(row_id %in% folds1)
    samp2 <- DF %>% filter(row_id %in% folds2)
    samp3 <- DF %>% filter(row_id %in% folds3)
    
    # train ps model on sample 1
    ps1_mod <- glm(Z~.-theta-Y-row_id, data = samp1, family='binomial')
    
    # predict ps on sample 2 and 3
    ps1 <- predict(ps1_mod, newdata = samp2, type='response')
    ps2 <- predict(ps1_mod, newdata = samp3, type='response')

    # train outcome on sample 2
    samp2_design <- makeModelMatrixFromDataFrame(samp2[, 4:(ncol(samp2)-1)])
    
    bcf_out <- bcf(y = samp2$Y, 
                   z = samp2$Z,
                   x_control = samp2_design, 
                   x_moderate = samp2_design,
                   pihat = ps1, 
                   nburn =1000, 
                   nsim =3000)
    
    samp3_design <- makeModelMatrixFromDataFrame(samp3[, 4:(ncol(samp2)-1)])
    pred_curr <- predict(bcf_out,
                         x_predict_control = samp3_design,
                         x_predict_moderate = samp3_design,
                         pi_pred = ps2,
                         z_pred = samp3$Z,
                         save_tree_directory = '..')
    
    new_res <- merge(x = res_mat, 
                     y = pred_curr, 
                     by='row_id', 
                     all.x=T)
    meanpred <- rowMeans(new_res[,2:3], na.rm = T)
    res_mat$ITE <- meanpred
  }
  return(res_mat)
}