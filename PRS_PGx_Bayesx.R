# PRS-PGx-Bayesx: Extend PRS-PGx-Bayes to cross-population

# PGx_GWAS: a numeric list containing PGx GWAS summary statistics of K populations
# G_reference: a numeric list containing the individual-level genotype matrix from the reference panel (e.g., 1KG) of K populations
# n.itr: a numeric value indicating the total number of MCMC iteration
# n.burnin: a numeric value indicating the number of burn in
# n.gap: a numeric value indicating the MCMC gap
# paras: a numeric vector containg hyper-parameters v (the degree of freedom) and phi (the global shrinkage parameter)
# standardize: a logical flag indicating should phenotype and genotype be standardized
# N: a numeric vector indicating the sample size in K PGx GWAS summary statistics
# K: a numeric value indicating the number of populations
# vali_Y: a numeric vector indicating the outcome in validation data
# vali_T: a numeric vector indicating the treatment assignment in validation data
# vali_G: a numeric matrix indicating the genotype in validation data
# tar_G: a numeric matrix indicating the genotype in target data

require(dplyr)
require(data.table)
require(bigsnpr)
require(Matrix)
require(GIGrvg)
require(MCMCpack)
require(bdsmatrix)
require(bigsparser)
require(lmtest)
require(mvtnorm)
require(propagate)
require(bigparallelr)
require(methods)
require(bigstatsr)
require(Rfast)
require(matrixcalc)

PRS_PGx_Bayesx <- function(PGx_GWAS, G_reference, n.itr = 1000, n.burnin = 500, n.gap = 10, paras, N, K, standardize = TRUE, vali_Y, vali_T, vali_G, tar_G){
  re <- PRS_PGx_Bayesx_sep(PGx_GWAS, G_reference, n.itr = 1000, n.burnin = 500, n.gap = 10, paras, N, K, standardize = TRUE)
  beta <- re$coef.G; alpha <- re$coef.TG
  
  PRS_G <- vali_G %*% beta
  PRS_GT <- vali_G %*% alpha
  
  df <- cbind.data.frame(Y=vali_Y, Tr=vali_T, PRS_G, vali_T*PRS_GT)
  colnames(df) <- c("Y","T",paste0("Pop",1:K),paste0("T_Pop",1:K))
  fit <- lm(Y ~ ., data=df)
  w <- fit$coefficients; w <- w[-c(1,2)]
  w_G <- w[c(1:K)]; w_GT <- w[c((K+1):(2*K))]
  
  PRS_G <- tar_G %*% beta %*% w_G %>% as.vector()
  PRS_GT <- tar_G %*% alpha %*% w_GT %>% as.vector()
  
  re <- list(PRS_G=PRS_G, PRS_GT=PRS_GT, beta=beta, alpha=alpha, w_G=w_G, w_GT=w_GT)
  re
}

PRS_PGx_Bayesx_sep <- function(PGx_GWAS, G_reference, n.itr = 1000, n.burnin = 500, n.gap = 10, paras, N, K, standardize = TRUE){
  # prepare scaling parameters for standarization and unstandardization
  sd.y.vec <- sd.G.mat <- sd.TG.mat <- double()
  for(k in 1:K){
    PGx_GWAS_k <- PGx_GWAS[[k]]
    MAF = PGx_GWAS_k$G_INFO$MAF
    sd.y = PGx_GWAS_k$sd.y; sd.y.vec <- c(sd.y.vec, sd.y)
    meanT = PGx_GWAS_k$meanT; varT = PGx_GWAS_k$varT
    mu.G = 2*MAF; var.G = 2*MAF*(1-MAF); sd.G = sqrt(var.G); sd.G.mat <- cbind(sd.G.mat, sd.G)
    var.TG = meanT^2*var.G + varT*mu.G^2 + varT*var.G; sd.TG = sqrt(var.TG); sd.TG.mat <- cbind(sd.TG.mat, sd.TG)
  }
  
  # prepare standardized b_hat = (beta_hat, alpha_hat)
  beta_hat <- alpha_hat <- double()
  for (k in 1:K) {
    PGx_GWAS_k <- PGx_GWAS[[k]]
    beta_hat_k <- PGx_GWAS_k$G_INFO$beta_hat; alpha_hat_k <- PGx_GWAS_k$G_INFO$alpha_hat
    beta_hat <- cbind(beta_hat, beta_hat_k)
    alpha_hat <- cbind(alpha_hat, alpha_hat_k)
  }
  if(standardize == TRUE){
    for (k in 1:K) {
      beta_hat[,k] <- beta_hat[,k]*sd.G.mat[,k]/sd.y.vec[k]; alpha_hat[,k] <- alpha_hat[,k]*sd.TG.mat[,k]/sd.y.vec[k]
    }
  }
  b_hat <- rbind(beta_hat, alpha_hat)
  
  # main Bayes
  D <- calculate_D(PGx_GWAS, G_reference, K)

  mod <- one_block_Bayes(D, b_hat, N, n.itr, n.burnin, n.gap, paras)
  coef.G <- mod$coef.G; coef.TG <- mod$coef.TG

  if(standardize == TRUE){
    for (k in 1:K) {
      coef.G[,k] <- coef.G[,k]*sd.y.vec[k]/sd.G.mat[,k]; coef.TG[,k] <- coef.TG[,k]*sd.y.vec[k]/sd.TG.mat[,k]
    }
  }
  rownames(coef.G) <- rownames(coef.TG) <- PGx_GWAS[[1]]$G_INFO$ID

  re <- list(coef.G = coef.G, coef.TG = coef.TG)
  return(re)
}

calculate_D <- function(PGx_GWAS, G_reference, K){
  D <- list()
  for (k in 1:K) {
    PGx_GWAS_k <- PGx_GWAS[[k]]; G_reference_k <- G_reference[[k]]
    
    M <- ncol(G_reference_k); G_reference_k <- as.matrix(G_reference_k)
    meanT = PGx_GWAS_k$meanT; varT = PGx_GWAS_k$varT
    MAF <- PGx_GWAS_k$G_INFO$MAF
    
    lefttop <- cora(G_reference_k); Sigma <- cova(G_reference_k)
    righttop <- matrix(NA, ncol = M, nrow = M)
    for (i in 1:M) {
      for (j in 1:M) {
        righttop[i,j] <- meanT*Sigma[i,j]/sqrt(Sigma[i,i]*(meanT^2*Sigma[j,j]+4*MAF[j]^2*varT+Sigma[j,j]*varT))
      }
    }
    
    leftbottom <- t(righttop)
    
    rightbottom <- matrix(NA, ncol = M, nrow = M)
    for (i in 1:M) {
      for (j in 1:M) {
        up <- (meanT^2+varT)*Sigma[i,j]+4*varT*MAF[i]*MAF[j]
        down <- sqrt((meanT^2*Sigma[i,i]+4*MAF[i]^2*varT+Sigma[i,i]*varT)*(meanT^2*Sigma[j,j]+4*MAF[j]^2*varT+Sigma[j,j]*varT))
        rightbottom[i,j] <- up/down
      }
    }
    Dk <- rbind(cbind(lefttop, righttop), cbind(leftbottom, rightbottom))
    
    D[[k]] <- Dk
  }
  return(D)
}

one_block_Bayes <- function(D, b_hat, N, n.itr, n.burnin, n.gap, paras){
  # initialization
  M <- nrow(b_hat)/2; K <- ncol(b_hat)
  b0 <- matrix(rep(0, 2*M*K),ncol = K); sigma20 <- rep(1,K)
  psi0 <- rep(1,M); xi0 <- rep(1,M); rho0 <- rep(0,M)
  cov0 <- rho0*sqrt(psi0*xi0); Omega0 <- rbind(cbind(diag(psi0),diag(cov0)),cbind(diag(cov0),diag(xi0)))
  delta0 <- rep(1, M); lambda0 <- rep(1, M)
  b1 <- b2 <- 1/2

  # Gibbs sampler
  n.pst <- (n.itr-n.burnin)/n.gap

  b.old <- b0; sigma2.old <- sigma20
  psi.old <- psi0; xi.old <- xi0; rho.old <- rho0; Omega.old <- Omega0
  delta.old <- delta0; lambda.old <- lambda0

  b.mat <- list()

  for (i in 1:n.itr) {
    stats <- update.once(b.old = b.old, sigma2.old = sigma2.old, psi.old = psi.old, xi.old = xi.old, rho.old = rho.old, Omega.old = Omega.old, delta.old = delta.old, lambda.old = lambda.old, D = D, b_hat = b_hat, paras = paras, N=N, M=M)

    b.old <- b.mat[[i]] <- stats$b.new

    sigma2.old <- stats$sigma2.new
    psi.old <- stats$psi.new; xi.old <- stats$xi.new; rho.old <- stats$rho.new
    Omega.old <- stats$Omega.new

    delta.old <- stats$delta.new; lambda.old <- stats$lambda.new
  }
  
  b.est <- double()
  for (k in 1:K) {
    b.est.k <- rep(0,2*M)
    for (i in 1:n.itr) {
      if(i > n.burnin && ((i-n.burnin) %% n.gap == 0)){
        b.est.k <- b.est.k + b.mat[[i]][,k]/n.pst
      }
    }
    b.est <- cbind(b.est, b.est.k)
  }

  coef.G <- b.est[c(1:M),]; coef.TG <- b.est[c((M+1):(2*M)),]

  re <- list(coef.G = coef.G, coef.TG = coef.TG, b.mat=b.mat)
  return(re)
}

update.once <- function(b.old, sigma2.old, psi.old, xi.old, rho.old, Omega.old, delta.old, lambda.old, D, b_hat, paras, N, M){
  # extract "paras"
  v <- paras[1]; phi <- paras[2]; b1 = b2 = 1/2

  # update b = (beta, alpha)
  b.new <- double()
  for (k in 1:K) {
    Psi.old.inv <- Psi.Inv(psi.old, xi.old, rho.old, M)
    D.Psi.old.inv <- solve(D[[k]]+Psi.old.inv, tol=1e-40)
    
    mu <- D.Psi.old.inv%*%b_hat[,k]; mu <- as.vector(mu)
    Sigma0 <- sigma2.old[k]/N[k]*D.Psi.old.inv; Sigma0 <- round(Sigma0, 5)
    if(is.positive.definite(Sigma0) == FALSE){Sigma0 <- nearPD(Sigma0, corr = FALSE); Sigma0 <- as.matrix(Sigma0$mat)}
    Sigma0_chol <- chol(Sigma0)
    b.new.k <- mu + as.vector(rnorm(length(mu))%*%Sigma0_chol)
    
    b.new <- cbind(b.new, b.new.k)
  }

  # update sigma2
  sigma2.new <- double()
  for (k in 1:K) {
    input.shape <- M + N[k]/2
    
    input.scale1 <- N[k]/2*(t(b.new[,k])%*%(D[[k]]+Psi.old.inv)%*%b.new[,k]+1-2*t(b_hat[,k])%*%b.new[,k])
    input.scale1 <- as.vector(input.scale1)
    input.scale2 <- N[k]/2*(t(b.new[,k])%*%Psi.old.inv%*%b.new[,k])
    input.scale2 <- as.vector(input.scale2)
    input.scale0 <- max(input.scale1, input.scale2)
    
    sigma2.new.k <- rinvgamma(1, shape = input.shape, scale = input.scale0)
    sigma2.new <- c(sigma2.new, sigma2.new.k)
  }

  # update Omega
  beta.new <- b.new[c(1:M),]; alpha.new <- b.new[c((M+1):(2*M)),]

  psi.new <- double(); xi.new <- double(); rho.new <- double()

  Comb.new <- sapply(1:M, generate.M, N, sigma2.new, beta.new, alpha.new, delta.old, lambda.old, paras)

  psi.new <- Comb.new[1,]; xi.new <- Comb.new[2,]; rho.new <- Comb.new[3,]

  corr <- rho.new*sqrt(psi.new*xi.new)
  Omega.new <- rbind(cbind(diag(psi.new),diag(corr)),cbind(diag(corr),diag(xi.new)))

  # update delta and lambda
  delta.new <- sapply(1:M, generate.delta, psi.new, rho.new, paras)
  lambda.new <- sapply(1:M, generate.lambda, xi.new, rho.new, paras)

  # posterior
  re <- list(b.new = b.new,
             sigma2.new = sigma2.new,
             psi.new = psi.new,
             xi.new = xi.new,
             rho.new = rho.new,
             Omega.new = Omega.new,
             delta.new = delta.new,
             lambda.new = lambda.new)
  return(re)
}

Psi.Inv <- function(psi,xi,rho,n){
  n <- length(psi)
  Psi_inv = matrix(0, ncol=2*n, nrow=2*n)

  rho11 = 1/(psi*(1-rho^2)); rho22 = 1/(xi*(1-rho^2))
  rho12 = -rho/(1-rho^2)/(sqrt(psi)*sqrt(xi))

  Psi_inv[row(Psi_inv)==(col(Psi_inv)-n)] = Psi_inv[col(Psi_inv)==(row(Psi_inv)-n)] = rho12
  diag(Psi_inv)[1:n] = rho11; diag(Psi_inv)[(n+1):(2*n)] = rho22

  return(Psi_inv)
}

generate.delta <- function(i, psi, rho, paras){
  v <- paras[1]; phi <- paras[2]; b1 = b2 = 1/2

  rate <- phi+(2*v)/(psi[i]*(1-rho[i]^2))
  re <- rgamma(1, shape = (v+b1+1/2), rate = rate)
  return(re)
}

generate.lambda <- function(i, xi, rho, paras){
  v <- paras[1]; phi <- paras[2]; b1 = b2 = 1/2

  rate <- phi+(2*v)/(xi[i]*(1-rho[i]^2))
  re <- rgamma(1, shape = (v+b2+1/2), rate = rate)
  return(re)
}

generate.M <- function(j, N, sigma2.new, beta.new, alpha.new, delta.old, lambda.old, paras){
  v <- paras[1]; phi <- paras[2]
  
  A <- matrix(c(0,0,0,0),nrow = 2)
  for (k in 1:K) {
    Ak <- N[k]/sigma2.new[k]*matrix(c(beta.new[j,k]^2, alpha.new[j,k]*beta.new[j,k], alpha.new[j,k]*beta.new[j,k], alpha.new[j,k]^2),ncol=2)
    A <- A + Ak
  }
  B <- diag(c(4*v*delta.old[j],4*v*lambda.old[j]))
  S <- A+B; S <- nearPD(S, corr = FALSE); S <- as.matrix(S$mat)

  Mj <- riwish(v = 2*v+K+1, S = S)

  psi.new <- ifelse(Mj[1,1] > 1/phi, 1/phi, Mj[1,1])
  xi.new <- ifelse(Mj[2,2] > 1/phi, 1/phi, Mj[2,2])
  rho.new <- Mj[1,2]/sqrt(Mj[1,1]*Mj[2,2])

  re <- c(psi.new, xi.new, rho.new); names(re) <- c("psi","xi","rho")
  return(re)
}




