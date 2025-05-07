logR_llh_lsjm_classicSingle <- function(param,hazard_baseline_01, sharedtype_01, func_sharedtype_01 = NULL,
                                   ord.splines, nb.beta, Zq, nb_pointsGK,
                                   nb.e.a, S, wk, rep_wk, sk_GK, nb.alpha,
                                   Matrices,
                                   left_trunc,
                                   knots.hazard_baseline.splines_01,
                                   index_beta_slope, index_b_slope, Ind,
                                   control){


  # Initialiser certains paramètres
  Gompertz.1_01 <- 0; Gompertz.2_01 <- 0;
  shape_01 <- 0;
  alpha.current_01 <- 0;
  alpha.slope_01 <- 0;
  alpha_01 <- c(0);
  gamma_01 <- c(0,3);
  beta_slope <- c(0); b_y_slope <- as.matrix(1);
  Xslope_T_i <- c(0); Uslope_T_i <- c(0); Xslope_GK_T_i <- as.matrix(1); Uslope_GK_T_i <- as.matrix(1);
  X_GK_T0_i<- as.matrix(1); U_GK_T0_i<- as.matrix(1); Xslope_GK_T0_i<- as.matrix(1); Uslope_GK_T0_i<- as.matrix(1);
  st_T0_i <- c(0); B_T_i_01 <- c(0);
  Bs_T_i_01 <- as.matrix(1);
  Bs_T0_i_01 <- as.matrix(1); Time_T0_i <- 0;
  st_T_i <- c(0); st_T0_i <- c(0); alpha_b_01 <- c(0);
  #Manage parameter
  curseur <- 1
  ## Risque 01
  ### Hazard baseline
  if(hazard_baseline_01 == "Weibull"){
    shape_01 <- param[curseur]**2
    curseur <- curseur + 1
  }
  if(hazard_baseline_01 == "Gompertz"){
    Gompertz.1_01 <- param[curseur]**2
    Gompertz.2_01 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(hazard_baseline_01 == "Splines"){
    gamma_01 <- param[(curseur):(curseur+ord.splines[1]+1)]
    curseur <- curseur + ord.splines[1] + 2
  }
  ### Covariables :
  nb.alpha_01 <- nb.alpha[1]
  if(nb.alpha_01 >=1){
    alpha_01 <-  param[(curseur):(curseur+nb.alpha_01-1)]
    curseur <- curseur+nb.alpha_01
  }
  ### Association
  if("random effects" %in% sharedtype_01){
    alpha_b_01 <- param[curseur:(curseur+nb.e.a-1)]
    curseur <- curseur + nb.e.a
  }
  if("value" %in% sharedtype_01){
    alpha.current_01 <-  param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% sharedtype_01){
    alpha.slope_01 <- param[curseur]
    curseur <- curseur + 1
  }
  ## Marker
  ### Fiexd effects
  beta <- param[curseur:(curseur+nb.beta-1)]
  curseur <- curseur+nb.beta
  sigma_epsilon <- param[curseur]
  curseur <- curseur +1
  C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
  C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
  MatCov <- C1
  MatCov <- as.matrix(MatCov)
  random.effects <- Zq%*%t(MatCov)
  b_y <- random.effects[,1:nb.e.a]
  b_y <- matrix(b_y, ncol = nb.e.a)
  if("slope" %in% sharedtype_01){
    b_y_slope <- as.matrix(b_y[,index_b_slope], ncol = length(index_b_slope))
    beta_slope <- beta[index_beta_slope]
  }
  ll_glob <- rep(NA, Ind)

  # Creations entrees rcpp
  sharedtype <- c("value" %in% sharedtype_01, "slope" %in% sharedtype_01,"random effects" %in% sharedtype_01
  )
  HB <- list(hazard_baseline_01)
  Weibull <- c(shape_01)
  Gompertz <- c(Gompertz.1_01, Gompertz.2_01)
  alpha_y_slope <- c(alpha.current_01, alpha.slope_01)
  alpha_z <- list(alpha_01)
  nb_points_integral <- c(S, nb_pointsGK)
  gamma_z0 <- list(gamma_01)

  func_assoc_01_type_int <- 0 # 0 pour linéaire (défaut)
  if(!is.null(func_sharedtype_01) && func_sharedtype_01 == "quad"){
    func_assoc_01_type_int <- 1 # 1 pour quadratique
  }


  delta1 <- Matrices[["delta1"]];  Z_01 <- Matrices[["Z_01"]]
  Time_T <- Matrices[["Time_T"]];
  st_T <- Matrices[["st_T"]];
  X_GK_T <- Matrices[["X_GK_T"]];U_GK_T <- Matrices[["U_GK_T"]]
  Xslope_GK_T <- Matrices[["Xslope_GK_T"]];Uslope_GK_T <- Matrices[["Uslope_GK_T"]]
  X_T <- Matrices[["X_T"]]; U_T <- Matrices[["U_T"]]; Xslope_T <- Matrices[["Xslope_T"]]; Uslope_T <- Matrices[["Uslope_T"]]
  X_base <- Matrices[["X_base"]]; U_base <- Matrices[["U_base"]]; y.new <- Matrices[["y.new"]]
  offset <- Matrices[["offset"]];
  B_T_01 <- Matrices[["B_T_01"]]; Bs_T_01 <- Matrices[["Bs_T_01"]];

    Time_T0 <- Matrices[["Time_T0"]]; st_T0 <- Matrices[["st_T0"]]; Bs_T0_01 <- Matrices[["Bs_T0_01"]]
    X_GK_T0 <- Matrices[["X_GK_T0"]];U_GK_T0 <- Matrices[["U_GK_T0"]]
    Xslope_GK_T0 <- Matrices[["Xslope_GK_T0"]];Uslope_GK_T0 <- Matrices[["Uslope_GK_T0"]]
    Bs_T0_01 <- Matrices[["Bs_T0_01"]]

  ll_glob <- log_llh_lsjm_classicSingle(sharedtype, HB, Gompertz, Weibull,
                                        nb_points_integral,
                                        func_assoc_01_type_int,
                                        alpha_y_slope,t(alpha_b_01), alpha_z,
                                        gamma_z0,beta, beta_slope, b_y, b_y_slope,
                                        wk, sigma_epsilon, delta1, Z_01,  X_T,  U_T,
                                        Xslope_T,  Uslope_T,  X_GK_T,  U_GK_T,  Xslope_GK_T,
                                        Uslope_GK_T,
                                        X_GK_T0,  U_GK_T0,  Xslope_GK_T0,  Uslope_GK_T0,
                                        Time_T,  Time_T0, st_T,  st_T0,
                                        Bs_T0_01,   left_trunc,
                                        X_base,  U_base,   y.new,  Ind,
                                        offset, B_T_01 , Bs_T_01

  )


  ll_glob2 <- sum(ll_glob)
  if(is.na(ll_glob2) || ll_glob2>0){
    #print(param)
    #print(ll_glob2)
    ll_glob2 <- -1E09
  }
  ll_glob2

}
