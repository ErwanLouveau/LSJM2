#' lsjm : Estimation of a location-scale joint model for longitudinal data with a flexible subject-specific variability and a complex survival design.
#'
#' This function fits joint models in which
#' we suppose that the variance of the residual error is subject-specific and the survival submodel contains either a single event, or competing risks or an illness-death model.
#' Nine differents models can be estimated (see details below).
#' Parameters are estimated through a maximum likelihood method, using a Marquardt-Levenberg algorithm.
#'
#' @param Objectlsmm The result of the lsmm function
#' @param survival_type a character defining which type of survival scheme to use either 'Single', or 'CR' or 'IDM'
#' @param formSurv_01 one-sided formula providing on the right-side the regression variable of the risk function for transition 0-1
#' @param formSurv_02 one-sided formula providing on the right-side the regression variable of the risk function for transition 0-2
#' @param formSurv_12 one-sided formula providing on the right-side the regression variable of the risk function for transition 1-2
#' @param sharedtype_01 a vector indicating the form(s) of the dependence structure. If the longitudinal model estimated is a standard mixed model, it should be included in c("value", "slope"). If it is a location-scale mixed model with a covariate and time-dependent variability one can add "current variability" and if it a model distinguishing within from between visits variabilities one can add c("variability inter", "variability intra") if there are subject-specifics.
#' @param sharedtype_02 a vector indicating the form(s) of the dependence structure. If the longitudinal model estimated is a standard mixed model, it should be included in c("value", "slope"). If it is a location-scale mixed model with a covariate and time-dependent variability one can add "current variability" and if it a model distinguishing within from between visits variabilities one can add c("variability inter", "variability intra") if there are subject-specifics.
#' @param sharedtype_12 a vector indicating the form(s) of the dependence structure. If the longitudinal model estimated is a standard mixed model, it should be included in c("value", "slope"). If it is a location-scale mixed model with a covariate and time-dependent variability one can add "current variability" and if it a model distinguishing within from between visits variabilities one can add c("variability inter", "variability intra") if there are subject-specifics.
#' @param hazardBase_01 a character providing the baseline hazard function, which is in c("Exponential", "Weibull", "Gompertz", "Splines")} of the risk function for transition 0-1
#' @param hazardBase_02 a character providing the baseline hazard function, which is in c("Exponential", "Weibull", "Gompertz", "Splines")} of the risk function for transition 0-2
#' @param hazardBase_12 a character providing the baseline hazard function, which is in c("Exponential", "Weibull", "Gompertz", "Splines")} of the risk function for transition 1-2
#' @param delta1 a one-sided formula with on the right-side the variable corresponding to the indicator of interest event (1 for subject experimenting the event and 0 for others).
#' @param delta2 a one-sided formula with on the right-side the variable corresponding to the indicator $\delta_{1i}$ of second event (1 for subject experimenting the event and 0 for others).
#' @param Time_T a one-sided formula with on the right-side the variable corresponding to the time of event
#' @param Time_L In case of an IDM model: a one-side formula with a variable giving the time of last time in state (0)
#' @param Time_R In case of an IDM model: a one-side formula with a variable providing the first time time known in state (1)
#' @param Time_T0 In case of delayed entry: is a a one-sided formula with on the right-side the variable corresponding to the time of entry
#' @param formSlopeFixed a one-sided formula corresponding to the derivative with respect to time of the function of fixed effects (if 'slope' is in sharedtype)
#' @param formSlopeRandom a one-sided formula corresponding to the derivative with respect to time of the function of random effects (if 'slope' is in sharedtype)
#' @param index_beta_slope a vector providing the index of beta parameters from the mixed model used in the slope
#' @param index_b_slope a vector providing the index of b parameters from the mixed model used in the slope
#' @param nb.knots.splines a vector providing the number of internal knots for Spline base hazard functions
#' @param nb_pointsGK the number of points for Gauss-Kronrod approximation : choice between 7 and 15. 15 by default.
#' @param S1 An integer : the number of QMC draws for the first step
#' @param S2 An integer : the number of QMC draws for the second step
#' @param binit optional initials parameters.
#' @param nproc An integer : the number of processors for parallel computing
#' @param clustertype one of the supported types from \code{makeCluster} function
#' @param maxiter
#' @param print.info
#' @param file
#' @param epsa
#' @param epsb
#' @param epsd
#'
#' @import survival
#' @import flexsurv
#' @import marqLevAlg
#' @import splines
#' @importFrom survival Surv
#' @export
#'
#'





lsjm <- function(Objectlsmm,  survival_type = c('Single', 'CR', 'IDM'),
                 formSurv_01, formSurv_02 = NULL, formSurv_12 = NULL,
                 sharedtype_01, sharedtype_02 = NULL, sharedtype_12 = NULL,
                 func_sharedtype_01 = NULL, func_sharedtype_02 = NULL, func_sharedtype_12 = NULL,
                 hazardBase_01, hazardBase_02 = NULL, hazardBase_12 = NULL,
                 delta1, delta2 = NULL, Time_T, Time_L = NULL, Time_R = NULL, Time_T0 = NULL,
                 formSlopeFixed = NULL, formSlopeRandom = NULL, index_beta_slope = NULL, index_b_slope = NULL, nb.knots.splines = c(1,1,1),
                 nb_pointsGK = 15, S1 = 1000, S2 = 5000, binit = NULL,
                 nproc = 1, clustertype = "SOCK", maxiter = 500,
                 print.info = FALSE, file = NULL, epsa = 1e-04, epsb = 1e-04, epsd = 1e-04
){
  if(missing(Objectlsmm)) stop("The argument Objectlsmm must be specified")
  if(!inherits(Objectlsmm, "lsmm_classic") && !inherits(Objectlsmm, "lsmm_covDep") && !inherits(Objectlsmm, "lsmm_interintra")) stop("use only \"lsmm\" objects")


  timeVar <- Objectlsmm$control$timeVar
  formVar <- Objectlsmm$control$formVar

  if(missing(survival_type)) stop("The argument survival_type must be specified")
  if(!inherits((survival_type),"character")) stop("The argument survival_type must be a character")
  if(length(survival_type) != 1) stop("The argument survival_type must be of length 1")
  if(!(survival_type %in% c('Single', 'CR', 'IDM'))) stop("The argument survival_type must be equal to 'Single' when there is only one event or 'CR' for competing risk model or 'IDM' for an illness-death model with possibly interval censoring")
  if(missing(formSurv_01)) stop("The argument formSurv_01 must be specified")
  if(!inherits((formSurv_01),"formula")) stop("The argument formSurv_01 must be a formula")
  if(survival_type %in% c('CR', 'IDM') && missing(formSurv_02)) stop("The argument formSurv_02 must be specified")
  if(survival_type %in% c('CR', 'IDM') && !inherits((formSurv_02),"formula")) stop("The argument formSurv_01 must be a formula")
  if(survival_type %in% c('IDM') && missing(formSurv_12)) stop("The argument formSurv_12 must be specified")
  if(survival_type %in% c('IDM') && !inherits((formSurv_12),"formula")) stop("The argument formSurv_12 must be a formula")

  if(!inherits((delta1),"formula")) stop("The argument delta1 must be a formula")
  if(survival_type %in% c('CR','IDM') && missing(delta2)) stop("The argument delta2 must be specified")
  if(survival_type %in% c('CR','IDM') && !inherits((delta2),"formula")) stop("The argument delta2 must be a formula")

  deltas <- list('delta1' = delta1,
                 'delta2' = delta2)

  if(missing(sharedtype_01)) stop("The argument sharedtype_01 must be specified")
  if(!inherits((sharedtype_01),"character")) stop("The argument sharedtype_01 must be a character")
  if(survival_type %in% c('CR', 'IDM') && missing(sharedtype_02)) stop("The argument sharedtype_02 must be specified")
  if(survival_type %in% c('CR', 'IDM') && !inherits((sharedtype_02),"character")) stop("The argument sharedtype_02 must be a character")
  if(survival_type %in% c('IDM') && missing(sharedtype_12)) stop("The argument sharedtype_12 must be specified")
  if(survival_type %in% c('IDM') && !inherits((sharedtype_12),"character")) stop("The argument sharedtype_12 must be a character")

  if(missing(hazardBase_01)) stop("The argument hazardBase_01 must be specified")
  if(!inherits((hazardBase_01),"character")) stop("The argument hazardBase_01 must be a character")
  if(survival_type %in% c('CR', 'IDM') && missing(hazardBase_02)) stop("The argument hazardBase_02 must be specified")
  if(survival_type %in% c('CR', 'IDM') && !inherits((hazardBase_02),"character")) stop("The argument hazardBase_02 must be a character")
  if(survival_type %in% c('IDM') && missing(hazardBase_12)) stop("The argument hazardBase_12 must be specified")
  if(survival_type %in% c('IDM') && !inherits((hazardBase_12),"character")) stop("The argument hazardBase_12 must be a character")
  if('slope' %in% c(sharedtype_01, sharedtype_02, sharedtype_12) && missing(formSlopeFixed)) stop("The argument formSlopeFixed must be specified")
  if('slope' %in% c(sharedtype_01, sharedtype_02, sharedtype_12) && !inherits((formSlopeFixed),"formula")) stop("The argument formSlopeFixed must be specified")
  if('slope' %in% c(sharedtype_01, sharedtype_02, sharedtype_12) && missing(formSlopeRandom)) stop("The argument formSlopeRandom must be specified")
  if('slope' %in% c(sharedtype_01, sharedtype_02, sharedtype_12) && !inherits((formSlopeRandom),"formula")) stop("The argument formSlopeRandom must be specified")
  if('slope' %in% c(sharedtype_01, sharedtype_02, sharedtype_12) && missing(index_beta_slope)) stop("The argument index_beta_slope must be specified")
  if('slope' %in% c(sharedtype_01, sharedtype_02, sharedtype_12) && missing(index_b_slope)) stop("The argument index_b_slope must be specified")

  Time <- list(Time_L = Time_L,
               Time_R = Time_R,
               Time_T = Time_T,
               Time_T0 = Time_T0)

  # Vérification du nouvel argument (simple vérification pour l'instant)
  if(!is.null(func_sharedtype_01) && func_sharedtype_01 != "quad"){
    warning("Currently only func_sharedtype_01 = 'quad' is supported. Using linear form.")
    func_sharedtype_01 <- NULL
  }


  if(formVar == 'classic'){
    if(!all(sharedtype_01 %in% c("value", "slope", "random effects"))) stop("The argument of sharedtype_01 must be in c('value', 'slope', 'random effects')")
    if(survival_type %in% c('CR', 'IDM') && !all(sharedtype_02 %in% c("value", "slope", "random effects"))) stop("The argument of sharedtype_02 must be in c('value', 'slope', 'random effects')")
    if(survival_type %in% c('IDM') && !all(sharedtype_12 %in% c("value", "slope", "random effects"))) stop("The argument of sharedtype_12 must be in c('value', 'slope', 'random effects')")
    if(survival_type == 'Single'){
      result <- lsjm_classicSingle(Objectlsmm, Time, deltas, hazardBase_01,  nb.knots.splines,
                                   formSurv_01,   nb_pointsGK, sharedtype_01, formSlopeFixed, formSlopeRandom,
                                   func_sharedtype_01 = func_sharedtype_01,
                                   index_beta_slope , index_b_slope, timeVar,
                                   S1, S2, binit, nproc , clustertype, maxiter,
                                   print.info, file , epsa , epsb , epsd )
    }
    else{
      if(survival_type == 'CR'){
        result <- lsjm_classicCR(Objectlsmm, Time, deltas, hazardBase_01, hazardBase_02, nb.knots.splines,
                                 formSurv_01, formSurv_02,  nb_pointsGK, sharedtype_01, sharedtype_02,
                                 formSlopeFixed, formSlopeRandom,
                                 index_beta_slope , index_b_slope, timeVar,
                                 S1, S2, binit, nproc , clustertype, maxiter,
                                 print.info, file , epsa , epsb , epsd)
      }
      else{
        result <- lsjm_classicIDM(Objectlsmm, Time, deltas, hazardBase_01, hazardBase_02, hazardBase_12, nb.knots.splines,
                                  formSurv_01, formSurv_02, formSurv_12, nb_pointsGK, sharedtype_01, sharedtype_02, sharedtype_12,
                                  formSlopeFixed, formSlopeRandom,
                                  index_beta_slope , index_b_slope, timeVar,
                                  S1, S2, binit, nproc , clustertype, maxiter,
                                  print.info, file , epsa , epsb , epsd)
      }
    }
  }
  else{
    if(formVar == 'cov-dependent'){
      if(!all(sharedtype_01 %in% c("value", "slope", "variability", "random effects"))) stop("The argument of sharedtype_01 must be in c('value', 'slope', 'variability', 'random effects')")
      if(survival_type %in% c('CR', 'IDM') && !all(sharedtype_02 %in% c("value", "slope", "variability", "random effects"))) stop("The argument of sharedtype_02 must be in c('value', 'slope', 'variability', 'random effects')")
      if(survival_type %in% c('IDM') && !all(sharedtype_12 %in% c("value", "slope", "variability", "random effects"))) stop("The argument of sharedtype_12 must be in c('value', 'slope', 'variability', 'random effects')")
      if(survival_type == 'Single'){
        result <- lsjm_covDepSingle(Objectlsmm, Time, deltas, hazardBase_01,  nb.knots.splines,
                                    formSurv_01,  nb_pointsGK, sharedtype_01,
                                    formSlopeFixed, formSlopeRandom,
                                    index_beta_slope , index_b_slope, timeVar,
                                    S1, S2, binit, nproc , clustertype, maxiter,
                                    print.info, file , epsa , epsb , epsd)
      }
      else{
        if(survival_type == 'CR'){
          result <- lsjm_covDepCR(Objectlsmm, Time, deltas, hazardBase_01, hazardBase_02, nb.knots.splines,
                                  formSurv_01, formSurv_02,  nb_pointsGK, sharedtype_01, sharedtype_02,
                                  formSlopeFixed, formSlopeRandom,
                                  index_beta_slope , index_b_slope,timeVar,
                                  S1, S2, binit, nproc , clustertype, maxiter,
                                  print.info, file , epsa , epsb , epsd)
        }
        else{
          result <- lsjm_covDepIDM(Objectlsmm, Time, deltas, hazardBase_01, hazardBase_02, hazardBase_12, nb.knots.splines,
                                   formSurv_01, formSurv_02, formSurv_12, nb_pointsGK, sharedtype_01, sharedtype_02, sharedtype_12,
                                   formSlopeFixed, formSlopeRandom,
                                   index_beta_slope , index_b_slope,timeVar,
                                   S1, S2, binit, nproc , clustertype, maxiter,
                                   print.info, file , epsa , epsb , epsd)
        }
      }
    }
    else{
      if(Objectlsmm$control$var_inter && Objectlsmm$control$var_intra){
        if(!all(sharedtype_01 %in% c("value", "slope", "variability inter", "variability intra", "random effects"))) stop("The argument of sharedtype_01 must be in c('value', 'slope', 'variability inter', 'variability intra', 'random effects')")
        if(survival_type %in% c('CR', 'IDM') && !all(sharedtype_02 %in% c("value", "slope", "variability inter", "variability intra", "random effects"))) stop("The argument of sharedtype_02 must be in c('value', 'slope', 'variability inter', 'variability intra', 'random effects')")
        if(survival_type %in% c('IDM') && !all(sharedtype_12 %in% c("value", "slope", "variability inter", "variability intra", "random effects"))) stop("The argument of sharedtype_12 must be in c('value', 'slope', 'variability inter', 'variability intra', 'random effects')")
      }
      else{
        if(Objectlsmm$control$var_inter){
          if(!all(sharedtype_01 %in% c("value", "slope", "variability inter", "random effects"))) stop("The argument of sharedtype_01 must be in c('value', 'slope', 'variability inter', 'random effects')")
          if(survival_type %in% c('CR', 'IDM') && !all(sharedtype_02 %in% c("value", "slope", "variability inter", "random effects"))) stop("The argument of sharedtype_02 must be in c('value', 'slope', 'variability inter', 'random effects')")
          if(survival_type %in% c('IDM') && !all(sharedtype_12 %in% c("value", "slope", "variability inter", "random effects"))) stop("The argument of sharedtype_12 must be in c('value', 'slope', 'variability inter', 'random effects')")
        }
        else{
          if(Objectlsmm$control$var_intra){
            if(!all(sharedtype_01 %in% c("value", "slope", "variability intra", "random effects"))) stop("The argument of sharedtype_01 must be in c('value', 'slope', 'variability intra', 'random effects')")
            if(survival_type %in% c('CR', 'IDM') && !all(sharedtype_02 %in% c("value", "slope", "variability intra", "random effects"))) stop("The argument of sharedtype_02 must be in c('value', 'slope', 'variability intra', 'random effects')")
            if(survival_type %in% c('IDM') && !all(sharedtype_12 %in% c("value", "slope", "variability intra", "random effects"))) stop("The argument of sharedtype_12 must be in c('value', 'slope', 'variability intra', 'random effects')")
          }
          else{
            if(Objectlsmm$control$var_intra){
              if(!all(sharedtype_01 %in% c("value", "slope", "random effects"))) stop("The argument of sharedtype_01 must be in c('value', 'slope', 'random effects')")
              if(survival_type %in% c('CR', 'IDM') && !all(sharedtype_02 %in% c("value", "slope", "random effects"))) stop("The argument of sharedtype_02 must be in c('value', 'slope', 'random effects')")
              if(survival_type %in% c('IDM') && !all(sharedtype_12 %in% c("value", "slope", "random effects"))) stop("The argument of sharedtype_12 must be in c('value', 'slope', 'random effects')")
            }
          }
        }
      }
      if(survival_type == 'Single'){
        result <- lsjm_interintraSingle(Objectlsmm, Time, deltas, hazardBase_01, nb.knots.splines,
                                        formSurv_01, nb_pointsGK, sharedtype_01,
                                        formSlopeFixed, formSlopeRandom,
                                        index_beta_slope , index_b_slope,timeVar,
                                        S1, S2, binit,nproc , clustertype, maxiter,
                                        print.info, file , epsa , epsb , epsd)
      }
      else{
        if(survival_type == 'CR'){
          result <- lsjm_interintraCR(Objectlsmm, Time, deltas, hazardBase_01, hazardBase_02, nb.knots.splines,
                                      formSurv_01, formSurv_02,  nb_pointsGK, sharedtype_01, sharedtype_02,
                                      formSlopeFixed, formSlopeRandom,
                                      index_beta_slope , index_b_slope,timeVar,
                                      S1, S2, binit,nproc , clustertype, maxiter,
                                      print.info, file , epsa , epsb , epsd)
        }
        else{
          result <- lsjm_interintraIDM(Objectlsmm, Time, deltas, hazardBase_01, hazardBase_02, hazardBase_12, nb.knots.splines,
                                       formSurv_01, formSurv_02, formSurv_12, nb_pointsGK, sharedtype_01, sharedtype_02, sharedtype_12,
                                       formSlopeFixed, formSlopeRandom,
                                       index_beta_slope , index_b_slope,timeVar,
                                       S1, S2, binit,nproc , clustertype, maxiter,
                                       print.info, file , epsa , epsb , epsd)
        }
      }

    }

  }
  result
}
