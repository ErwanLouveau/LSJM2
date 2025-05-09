#' @export

summary.lsjm_classicSingle <- function(object,...)
{
  x <- object
  if(!inherits(x, "lsjm_classicSingle")) stop("use only \"lsjm_classicSingle\" objects")

  if(x$result_step1$istop != 1){
    cat("--------------------------------------------------------------------------------------------------------- \n")
    cat("WARNING", "\n")
    cat("The first step estimation didn't reach convergence because :")
    if(x$result_step1$istop==2) cat("Maximum number of iteration reached without convergence.")
    if(x$result_step1$istop==4) {cat("The program stopped abnormally. No results can be displayed. \n")}
    cat("\n")
    cat("We recommend to not interpret the following results and to try to reach convergence for the first step. \n")

    cat("--------------------------------------------------------------------------------------------------------- \n")
    stop("\n")
  }



  cat("Joint model for a single event fitted by maximum likelihood method", "\n")

  #ajouter le code d'appelle à la fonction
  cat("\n")
  cat("Statistical Model:", "\n")
  cat(paste("    Number of subjects:", x$control$Ind),"\n")
  cat(paste("    Number of observations:", nrow(x$control$data.long)),"\n")

  cat("\n")
  cat("Iteration process:", "\n")

  if(!is.null(x$info_conv_step2)){
    if(x$info_conv_step2$conv==1) cat("    Convergence criteria satisfied")
    if(x$info_conv_step2$conv==2) cat("    Maximum number of iteration reached without convergence")
    if(x$info_conv_step2$conv==4) cat("    The program stopped abnormally. No results can be displayed. \n")
  }
  cat("\n")
  cat(paste("     Number of iterations: "), "\n")
  cat(paste("          Step 1: ",x$result_step1$ni), "\n")
  cat(paste("          Step 2: ",x$info_conv_step2$niter), "\n")
  cat(paste("     Convergence criteria (Step1): parameters =" ,signif(x$info_conv_step1$convcrit[1],3)), "\n")
  cat(paste("                                 : likelihood =" ,signif(x$info_conv_step1$convcrit[2],3)), "\n")
  cat(paste("                                 : second derivatives =" ,signif(x$info_conv_step1$convcrit[3],3)), "\n")
  cat(paste("     Time of computation :" , format(x$time.computation)),  "\n")

  cat("\n")
  cat("Goodness-of-fit statistics:")
  cat("\n")
  cat(paste("    Likelihood: ", round(x$result_step1$fn.value,3)),"\n")
  cat(paste("    AIC: ", round(2*length(x$result_step1$b) - 2* x$result_step1$fn.value,3)),"\n")

  cat("\n")
  cat("Maximum Likelihood Estimates:")
  #Manage parameters
  curseur <- 1
  param <- x$table.res$Estimation
  param.se <- x$table.res$SE
  param.names <- rownames(x$table.res)

  #Manage parameter
  curseur <- 1
  ## Risque 01
  ### Hazard baseline
  if(x$control$hazard_baseline_01 == "Weibull"){
    shape_01 <- param[curseur]
    shape_01.se <- param.se[curseur]
    shape_01.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if(x$control$hazard_baseline_01 == "Gompertz"){
    Gompertz.1_01 <- param[curseur]
    Gompertz.2_01 <- param[curseur+1]
    Gompertz.1_01.se <- param.se[curseur]
    Gompertz.2_01.se <- param.se[curseur+1]
    Gompertz.1_01.name <- param.names[curseur]
    Gompertz.2_01.name <- param.names[curseur+1]
    curseur <- curseur + 2
  }
  if(x$control$hazard_baseline_01 == "Splines"){
    gamma_01 <- param[(curseur):(curseur+x$control$nb.knots.splines[1]+2+1)]
    gamma_01.se <- param.se[(curseur):(curseur+x$control$nb.knots.splines[1]+2+1)]
    gamma_01.name <- param.names[(curseur):(curseur+x$control$nb.knots.splines[1]+2+1)]
    curseur <- curseur + x$control$nb.knots.splines[1]+2 + 2
  }
  ### Covariables :
  nb.alpha_01 <- x$control$nb.alpha[1]
  if(nb.alpha_01 >=1){
    alpha_01 <-  param[(curseur):(curseur+nb.alpha_01-1)]
    alpha_01.se <-  param.se[(curseur):(curseur+nb.alpha_01-1)]
    alpha_01.name <-  param.names[(curseur):(curseur+nb.alpha_01-1)]
    curseur <- curseur+nb.alpha_01
  }
  ### Association
  if("random effects" %in% x$control$sharedtype_01){
    alpha.re_01 <- param[curseur:(curseur+x$control$Objectlsmm$control$nb.e.a-1)]
    alpha.re_01.se <- param.se[curseur:(curseur+x$control$Objectlsmm$control$nb.e.a-1)]
    alpha.re_01.name <- param.names[curseur:(curseur+x$control$Objectlsmm$control$nb.e.a-1)]
    curseur <- curseur + x$control$Objectlsmm$control$nb.e.a
  }
  if("value" %in% x$control$sharedtype_01){
    alpha.current_01 <-  param[curseur]
    alpha.current_01.se <-  param.se[curseur]
    alpha.current_01.name <-  param.names[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% x$control$sharedtype_01){
    alpha.slope_01 <- param[curseur]
    alpha.slope_01.se <- param.se[curseur]
    alpha.slope_01.name <- param.names[curseur]
    curseur <- curseur + 1
  }


  ## Effets fixes trend :
  beta <- param[curseur:(curseur+x$control$Objectlsmm$control$nb.beta-1)]
  beta.se <- param.se[curseur:(curseur+x$control$Objectlsmm$control$nb.beta-1)]
  beta.name <- param.names[curseur:(curseur+x$control$Objectlsmm$control$nb.beta-1)]
  curseur <- curseur+x$control$Objectlsmm$control$nb.beta
  sigma_epsilon <- param[curseur]
  sigma_epsilon.se <- param.se[curseur]
  sigma_epsilon.name <- param.names[curseur]
  curseur <- curseur +1

  borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a, k = 2) + x$control$Objectlsmm$control$nb.e.a - 1
  C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a)**2),nrow=x$control$Objectlsmm$control$nb.e.a,ncol=x$control$Objectlsmm$control$nb.e.a)
  C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
  MatCov <- C1
  MatCov <- as.matrix(MatCov)
  borne2 <- borne1 + choose(n = x$control$Objectlsmm$control$nb.e.a, k = 2) + x$control$Objectlsmm$control$nb.e.a
  Matcov.name <- unique( unique(gsub("\\*.*", "", gsub("__", "_", param.names[(borne1+1):borne2]))))



  cat("\n")
  cat("Longitudinal model:")
  cat("\n")

  cat("      Fixed effects:")

  betas_tab <- matrix(nrow = length(beta), ncol = 4)
  betas_tab[,1] <- beta
  betas_tab[,2] <- beta.se
  betas_tab[,3] <- betas_tab[,1]/betas_tab[,2]
  betas_tab[,4] <- 1 - pchisq(betas_tab[,3]**2,1)
  betas_tab <- as.data.frame(betas_tab)
  rownames(betas_tab) <- gsub("_Y", "",beta.name)
  colnames(betas_tab) <- c("Coeff", "SE", "Wald", "Pvalue")
  betas_tab <- round(betas_tab, 4)
  betas_tab$Pvalue <- ifelse(betas_tab$Pvalue < 0.001, "<0.001", round(betas_tab$Pvalue,3))
  cat("\n")
  print(betas_tab)

  cat("\n")
  cat("      Residual variability:")
  sigma_eps_tab <- matrix(nrow = 1, ncol = 4)
  sigma_eps_tab[,1] <- sigma_epsilon
  sigma_eps_tab[,2] <- sigma_epsilon.se
  sigma_eps_tab[,3] <- sigma_eps_tab[,1]/sigma_eps_tab[,2]
  sigma_eps_tab[,4] <- 1 - pchisq(sigma_eps_tab[,3]**2,1)
  sigma_eps_tab <- as.data.frame(sigma_eps_tab)
  colnames(sigma_eps_tab) <- c("Coeff", "SE", "Wald", "Pvalue")
  sigma_eps_tab <- round(sigma_eps_tab, 4)
  sigma_eps_tab$Pvalue <- ifelse(sigma_eps_tab$Pvalue < 0.001, "<0.001", round(sigma_eps_tab$Pvalue,3))
  cat("\n")
  print(sigma_eps_tab)


  cat("\n")

  cat("     Covariance matrix of the random effects:")
  cat("\n")
  Cov <- MatCov%*%t(MatCov)
  colnames(Cov) <- Matcov.name
  rownames(Cov) <- Matcov.name
  print(Cov)
  cat("\n")

  cat("Survival model:")
  cat("\n")
  #browser()
  e1_var_tab <- NULL
  e1_share_random_tab <- NULL
  e1_share_current_tab <- NULL
  e1_share_slope_tab <- NULL
  e1_alpha_tab <- NULL
  e1_names_tab <- c()

  if(c("random effects") %in% x$control$sharedtype_01){
    e1_share_random_tab <- matrix(nrow = 1, ncol = 4)
    e1_share_random_tab[,1] <- alpha.re_01
    e1_share_random_tab[,2] <- alpha.re_01.se
    e1_share_random_tab[,3] <- e1_share_random_tab[,1]/e1_share_random_tab[,2]
    e1_share_random_tab[,4] <- 1 - pchisq(e1_share_random_tab[,3]**2,1)
    e1_names_tab <- c(e1_names_tab, alpha.re_01.name)
  }
  if(c("value") %in% x$control$sharedtype_01){
    e1_share_current_tab <- matrix(nrow = 1, ncol = 4)
    e1_share_current_tab[,1] <- alpha.current_01
    e1_share_current_tab[,2] <- alpha.current_01.se
    e1_share_current_tab[,3] <- e1_share_current_tab[,1]/e1_share_current_tab[,2]
    e1_share_current_tab[,4] <- 1 - pchisq(e1_share_current_tab[,3]**2,1)
    e1_names_tab <- c(e1_names_tab, alpha.current_01.name)
  }
  if(c("slope") %in% x$control$sharedtype_01){
    e1_share_slope_tab <- matrix(nrow = 1, ncol = 4)
    e1_share_slope_tab[,1] <- c(alpha.slope_01)
    e1_share_slope_tab[,2] <- c(alpha.slope_01.se)
    e1_share_slope_tab[,3] <- e1_share_slope_tab[,1]/e1_share_slope_tab[,2]
    e1_share_slope_tab[,4] <- 1 - pchisq(e1_share_slope_tab[,3]**2,1)
    e1_names_tab <- c(e1_names_tab, alpha.slope_01.name)
  }
  if(x$control$hazard_baseline_01 == "Splines"){
    if(x$control$nb.alpha[1] >=1){
      e1_alpha_tab <- matrix(nrow = length(alpha_01), ncol = 4)
      e1_alpha_tab[,1] <- alpha_01
      e1_alpha_tab[,2] <- alpha_01.se
      e1_alpha_tab[,3] <- e1_alpha_tab[,1]/e1_alpha_tab[,2]
      e1_alpha_tab[,4] <- 1 - pchisq(e1_alpha_tab[,3]**2,1)
      e1_names_tab <- c(e1_names_tab, alpha_01.name)
    }
  }
  else{
    if(x$control$nb.alpha[1] >=2){
      e1_alpha_tab <- matrix(nrow = length(alpha_01)-1, ncol = 4)
      e1_alpha_tab[,1] <- alpha_01[-1]
      e1_alpha_tab[,2] <- alpha_01.se[-1]
      e1_alpha_tab[,3] <- e1_alpha_tab[,1]/e1_alpha_tab[,2]
      e1_alpha_tab[,4] <- 1 - pchisq(e1_alpha_tab[,3]**2,1)
      e1_names_tab <- c(e1_names_tab, alpha_01.name[-1])
    }
  }


  e1_bas_tab <- NULL
  if(x$control$hazard_baseline_01 == "Exponential"){
    e1_bas_tab <- matrix(nrow = 1, ncol = 4)
    e1_bas_tab[1,1] <- alpha_01[1]
    e1_bas_tab[1,2] <- alpha_01.se[1]
    e1_bas_tab[1,3] <- e1_bas_tab[,1]/e1_bas_tab[,2]
    e1_bas_tab[1,4] <- 1 - pchisq(e1_bas_tab[,3]**2,1)
    #e1_names_tab <- c(e1_names_tab, alpha_01.name[-1])
    rownames(e1_bas_tab) <- c("intercept")
    colnames(e1_bas_tab) <- c("Coeff", "SE", "Wald", "Pvalue")
  }
  if(x$control$hazard_baseline_01 == "Weibull"){
    e1_bas_tab <- matrix(nrow = 2, ncol = 4)
    e1_bas_tab[1,1] <- alpha_01[1]
    e1_bas_tab[1,2] <- alpha_01.se[1]
    #e1_bas_tab[1,] <- e1_alpha_tab[1,]
    #e1_alpha_tab <- e1_alpha_tab[-1,]
    #e1_names_tab <- c(e1_names_tab, alpha.name[-1])
    e1_bas_tab[2,1] <- shape_01
    e1_bas_tab[2,2] <- shape_01.se
    e1_bas_tab[,3] <- e1_bas_tab[,1]/e1_bas_tab[,2]
    e1_bas_tab[,4] <- 1 - pchisq(e1_bas_tab[,3]**2,1)
    rownames(e1_bas_tab) <- c("intercept",shape_01.name)
    colnames(e1_bas_tab) <- c("Coeff", "SE", "Wald", "Pvalue")
  }
  if(x$control$hazard_baseline_01 == "Splines"){
    e1_bas_tab <- matrix(nrow = length(gamma_01), ncol = 4)
    e1_bas_tab[,1] <- gamma_01
    e1_bas_tab[,2] <- gamma_01.se
    e1_bas_tab[,3] <- e1_bas_tab[,1]/e1_bas_tab[,2]
    e1_bas_tab[,4] <- 1 - pchisq(e1_bas_tab[,3]**2,1)
    rownames(e1_bas_tab) <- gamma_01.name
    colnames(e1_bas_tab) <- c("Coeff", "SE", "Wald", "Pvalue")
  }
  e1_surv_tab <- rbind(e1_var_tab, e1_share_random_tab, e1_share_current_tab, e1_share_slope_tab, e1_alpha_tab)
  rownames(e1_surv_tab) <- e1_names_tab
  colnames(e1_surv_tab) <- c("Coeff", "SE", "Wald", "Pvalue")

  e1_surv_tab <- round(e1_surv_tab, 4)
  e1_surv_tab[,4] <- ifelse(e1_surv_tab[,4] < 0.001, "<0.001", round(e1_surv_tab[,4],3))
  e1_bas_tab <- round(e1_bas_tab, 4)
  e1_bas_tab[,4] <- ifelse(e1_bas_tab[,4] < 0.001, "<0.001", round(e1_bas_tab[,4],3))


  if(nrow(e1_bas_tab)!=0){
    cat("\n")
    cat("       Regression:")
    cat("\n")
    print(e1_surv_tab)
  }
  cat("\n")
  cat(paste("     Baseline: ",x$control$hazard_baseline_01), "\n")
  cat("\n")
  print(e1_bas_tab)



  cat("\n")


  if(x$result_step1$istop != 1){
    cat("--------------------------------------------------------------------------------------------------------- \n")
    cat("WARNING", "\n")
    cat("The first step estimation didn't reach convergence because :")
    if(x$result_step1$istop==2) cat("Maximum number of iteration reached without convergence.")
    if(x$result_step1$istop==4) {cat("The program stopped abnormally. No results can be displayed. \n")}
    cat("\n")
    cat("We recommend to not interpret the previous results and to try to reach convergence for the first step. \n")

    cat("--------------------------------------------------------------------------------------------------------- \n")
    cat("\n")
  }
}
