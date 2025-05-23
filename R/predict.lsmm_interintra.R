#' @rdname predict
#' @export
#'

predict.lsmm_interintra <- function(object, which = "RE", Objectranef = NULL, data.long = NULL){

  Objectlsmm <- object
  if(missing(Objectlsmm)) stop("The argument Objectlsmm must be specified")
  if(!inherits((Objectlsmm),"lsmm_interintra")) stop("use only \"lsmm_interintra\" objects")
  #if(missing(data.long)) stop("The argument data.long must be specified")
  #if(!inherits((data.long),"data.frame")) stop("use only \"data.frame\" objects")
  if(missing(which)) stop("The argument which must be specified")
  if(!inherits((which),"character")) stop("The argument which must be a character object")

  x <- Objectlsmm
  if(x$result_step1$istop != 1|| (!is.null(x$result_step2) && x$result_step2$istop !=1)){
    stop("The model didn't reach convergence.")
  }
  if(is.null(x$result_step2)){
    param <- x$result_step1$b
  }
  else{
    param <- x$result_step2$b
  }

  #x$control$nproc <- 1

  mu.inter <- 0; sigma.epsilon.inter <-0; mu.intra <- 0;sigma.epsilon.intra <- 0
  curseur <- 1
  ### Fixed effects
  beta <- param[curseur:(curseur+x$control$nb.beta-1)]
  curseur <- curseur+x$control$nb.beta
  if(x$control$var_inter){
    mu.inter <- param[curseur]
    curseur <- curseur + 1
  }
  else{
    sigma.epsilon.inter <- param[curseur]
    curseur <- curseur +1
  }
  if(x$control$var_intra){
    mu.intra <- param[curseur]
    curseur <- curseur + 1
  }
  else{
    sigma.epsilon.intra <- param[curseur]
    curseur <- curseur +1
  }


  ### Cholesky matrix for random effects => A revoir
  if(x$control$var_inter && x$control$var_intra){
    if(x$control$correlated_re){

      C1 <- matrix(rep(0,(x$control$nb.e.a+2)**2),nrow=x$control$nb.e.a+2,ncol=x$control$nb.e.a+2)
      C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
      Cholesky <- C1
    }
    else{
      borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
      C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      C2 <-matrix(c(param[(borne1+1)], 0,param[borne1+2], param[borne1+3]),nrow=2,ncol=2, byrow = TRUE)
      C3 <- matrix(rep(0,2*x$control$nb.e.a), ncol = x$control$nb.e.a)
      C4 <- matrix(rep(0,2*x$control$nb.e.a), nrow = x$control$nb.e.a)
      Cholesky <- rbind(cbind(C1,C4),cbind(C3,C2))
      Cholesky <- as.matrix(Cholesky)
    }
  }
  else{
    if(x$control$var_inter){
      if(x$control$correlated_re){
        C1 <- matrix(rep(0,(x$control$nb.e.a+1)**2),nrow=x$control$nb.e.a+1,ncol=x$control$nb.e.a+1)
        C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
        Cholesky <- C1
      }
      else{
        borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
        C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
        C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
        C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
        C3 <- matrix(rep(0,x$control$nb.e.a), ncol = 1)
        C4 <- matrix(rep(0,x$control$nb.e.a), nrow = 1)
        Cholesky <- rbind(cbind(C1,C3),cbind(C4,C2))
        Cholesky <- as.matrix(Cholesky)
      }
    }
    else{
      if(x$control$var_intra){
        if(x$control$correlated_re){
          C1 <- matrix(rep(0,(x$control$nb.e.a+1)**2),nrow=x$control$nb.e.a+1,ncol=x$control$nb.e.a+1)
          C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
          Cholesky <- C1
        }
        else{
          borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
          C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
          C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
          C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
          C3 <- matrix(rep(0,x$control$nb.e.a), ncol = 1)
          C4 <- matrix(rep(0,x$control$nb.e.a), nrow = 1)
          Cholesky <- rbind(cbind(C1,C3),cbind(C4,C2))
          Cholesky <- as.matrix(Cholesky)
        }
      }
    }
    if(!x$control$var_inter && !x$control$var_intra){
      C1 <- matrix(rep(0,(length(param)-curseur)**2),nrow=length(param)-curseur,ncol=length(param)-curseur)
      C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
      Cholesky <- C1
    }

  }

  MatCov <- Cholesky%*%t(Cholesky)
  if(is.null(data.long)){
    data.long <- x$control$data.long
  }
  if(x$control$var_inter && x$control$var_intra){
    random.effects.Predictions <- matrix(NA, nrow = length(unique(data.long$id)), ncol = x$control$nb.e.a+2+1)
    binit <- matrix(0, nrow = 1, ncol = x$control$nb.e.a+2)
  }
  else{
    if(x$control$var_inter || x$control$var_intra){
      random.effects.Predictions <- matrix(NA, nrow = length(unique(data.long$id)), ncol = x$control$nb.e.a+1+1)
      binit <- matrix(0, nrow = 1, ncol = x$control$nb.e.a+1)
    }
    else{
      random.effects.Predictions <- matrix(NA, nrow = length(unique(data.long$id)), ncol = x$control$nb.e.a+1)
      binit <- matrix(0, nrow = 1, ncol = x$control$nb.e.a)
    }
  }

  data.id <- data.long[!duplicated(data.long$id),]
  list.long <- data.manag.long(x$control$formGroup,x$control$formFixed, x$control$formRandom,data.long)
  X_base <- list.long$X; U_base <- list.long$U; y.new <- list.long$y.new
  time.measures <- data.long[,x$control$timeVar]
  ID.visit <- data.long[all.vars(x$control$formGroupVisit)][,1]; offset <- list.long$offset

  offset_ID <- c()
  len_visit <- c(0)
  for(oo in 1:length(unique(data.long$id))){
    ID.visit_i <- ID.visit[offset[oo]:(offset[oo+1]-1)]
    offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
    len_visit <- c(len_visit,length(unique(ID.visit_i)))
    frame_offset_ID_i <- cbind(offset_ID_i, rep(oo, length(offset_ID_i)))
    offset_ID <- rbind(offset_ID, frame_offset_ID_i)
  }
  offset_position <- as.vector(c(1, 1 + cumsum(tapply(offset_ID[,2], offset_ID[,2], length))))

  cv.Pred <- c()

  if(is.null(Objectranef) || ('RE' %in% which)){

    n_cores <- min(x$control$nproc,detectCores() - 1)   # Utiliser tous les cœurs sauf 1
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    results <- foreach(id.boucle = 1:length(unique(data.long$id)),
                       .combine = function(...) {
                         lst <- list(...)
                         list(random.effects.Predictions = do.call(rbind, lapply(lst, `[[`, 1)),
                              cv.Pred = do.call(rbind, lapply(lst, `[[`, 2)))
                       },
                       .multicombine = TRUE,
                       .packages = c("mvtnorm", "marqLevAlg")) %dopar% {

                         # Extraire les données spécifiques à l'individu
                         X_base_i <- X_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
                         X_base_i <- matrix(X_base_i, nrow = offset[id.boucle+1]-offset[id.boucle])
                         X_base_i <- unique(X_base_i)
                         U_i <- U_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
                         U_i <- matrix(U_i, nrow = offset[id.boucle+1]-offset[id.boucle])
                         U_base_i <- unique(U_i)
                         y_i <- y.new[offset[id.boucle]:(offset[id.boucle+1]-1)]
                         ID.visit_i <- ID.visit[offset[id.boucle]:(offset[id.boucle+1]-1)]
                         offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
                         len_visit_i <- length(unique(ID.visit_i))

                         # Estimation des effets aléatoires
                         random.effects_i <- marqLevAlg(binit, fn = re_lsmm_interintra, minimize = FALSE, nb.e.a = x$control$nb.e.a, variability_inter_visit = x$control$var_inter,
                                                        variability_intra_visit = x$control$var_intra, Sigma.re = MatCov,
                                                        beta = beta,
                                                        mu.inter = mu.inter , sigma.epsilon.inter = sigma.epsilon.inter, mu.intra = mu.intra,sigma.epsilon.intra = sigma.epsilon.intra,
                                                        len_visit_i = len_visit_i, X_base_i = X_base_i, U_base_i = U_base_i, y_i = y_i, offset_ID_i = offset_ID_i,
                                                        nproc = 1, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                                        file = "", blinding = FALSE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4)

                         while(random.effects_i$istop !=1){
                           binit <- mvtnorm::rmvnorm(1, mean = rep(0, ncol(MatCov)), MatCov)
                           random.effects_i <- marqLevAlg(binit, fn = re_lsmm_interintra, minimize = FALSE, nb.e.a = x$control$nb.e.a, variability_inter_visit = x$control$var_inter,
                                                          variability_intra_visit = x$control$var_intra, Sigma.re = MatCov,
                                                          beta = beta,
                                                          mu.inter = mu.inter , sigma.epsilon.inter = sigma.epsilon.inter, mu.intra = mu.intra,sigma.epsilon.intra = sigma.epsilon.intra,
                                                          len_visit_i = len_visit_i, X_base_i = X_base_i, U_base_i = U_base_i, y_i = y_i, offset_ID_i = offset_ID_i,
                                                          nproc = 1, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                                          file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
                           if(x$control$var_inter && x$control$var_intra){
                             binit <-matrix(0, nrow = 1, ncol = x$control$nb.e.a+2)
                           }
                           else{
                             if(x$control$var_inter || x$control$var_intra){
                               binit <-matrix(0, nrow = 1, ncol = x$control$nb.e.a+1)
                             }
                             else{
                               binit <-matrix(0, nrow = 1, ncol = x$control$nb.e.a)
                             }
                           }
                         }


                         random_effects_pred <- c(data.id$id[id.boucle], random.effects_i$b)
                         if ('Y' %in% which) {
                           CV <- X_base_i %*% beta + U_base_i %*% random.effects_i$b[1:(x$control$nb.e.a)]
                           time.measures_i <- unique(time.measures[offset[id.boucle]:(offset[id.boucle+1]-1)])
                           if(x$control$var_inter && x$control$var_intra){
                             Varia.inter <- exp(mu.inter + random.effects_i$b[x$control$nb.e.a+1])
                             Varia.intra <- exp(mu.intra + random.effects_i$b[x$control$nb.e.a+2])
                           }
                           else{
                             if(x$control$var_inter){
                               Varia.inter <- exp(mu.inter + random.effects_i$b[x$control$nb.e.a+1])
                               Varia.intra <- sigma.epsilon.intra
                             }
                             else{
                               if(x$control$var_intra){
                                 Varia.intra <- exp(mu.intra + random.effects_i$b[x$control$nb.e.a+1])
                                 Varia.inter <- sigma.epsilon.inter
                               }
                               else{
                                 Varia.intra <- sigma.epsilon.intra
                                 Varia.inter <- sigma.epsilon.inter
                               }
                             }
                           }

                           cv_pred <- cbind(rep(data.id$id[id.boucle], length(CV)),
                                                           time.measures_i, CV, rep(Varia.inter,length(time.measures_i)), rep(Varia.intra,length(time.measures_i)))
                         } else {
                           cv_pred <- NULL
                         }

                         # Retourner les deux résultats
                         list(random_effects_pred, cv_pred)

                       }

    # Extraire les deux tables finales
    random.effects.Predictions <- results$random.effects.Predictions
    cv.Pred <- results$cv.Pred

    # Fermer le cluster
    stopCluster(cl)

  }
  else{
    if('Y' %in% which){
      for(id.boucle in 1:length(unique(data.long$id))){
        X_base_i <- X_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
        X_base_i <- matrix(X_base_i, nrow = offset[id.boucle+1]-offset[id.boucle])
        X_base_i <- unique(X_base_i)
        U_i <- U_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
        U_i <- matrix(U_i, nrow = offset[id.boucle+1]-offset[id.boucle])
        U_base_i <- unique(U_i)
        y_i <- y.new[offset[id.boucle]:(offset[id.boucle+1]-1)]
        ID.visit_i <- ID.visit[offset[id.boucle]:(offset[id.boucle+1]-1)]
        offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
        len_visit_i <- length(unique(ID.visit_i))


        random.effects_i <- as.matrix(Objectranef[id.boucle,-1], nrow = 1)
        CV <- X_base_i %*% beta + U_base_i %*% random.effects_i[1,1:(x$control$nb.e.a)]
        time.measures_i <- unique(time.measures[offset[id.boucle]:(offset[id.boucle+1]-1)])

        if(x$control$var_inter && x$control$var_intra){
          Varia.inter <- exp(mu.inter + random.effects_i[x$control$nb.e.a+1])
          Varia.intra <- exp(mu.intra + random.effects_i[x$control$nb.e.a+2])
        }
        else{
          if(x$control$var_inter){
            Varia.inter <- exp(mu.inter + random.effects_i[x$control$nb.e.a+1])
            Varia.intra <- sigma.epsilon.intra
          }
          else{
            if(x$control$var_intra){
              Varia.intra <- exp(mu.intra + random.effects_i[x$control$nb.e.a+1])
              Varia.inter <- sigma.epsilon.inter
            }
            else{
              Varia.intra <- sigma.epsilon.intra
              Varia.inter <- sigma.epsilon.inter
            }
          }
        }
        cv_pred <- cbind(rep(data.id$id[id.boucle], length(CV)),
                                        time.measures_i, CV, rep(Varia.inter,length(time.measures_i)), rep(Varia.intra, length(time.measures_i)))



        cv.Pred <- rbind(cv.Pred,cv_pred)

      }


    }
  }

  if('RE' %in% which){
    random.effects.Predictions <- as.data.frame(random.effects.Predictions)
    name_b <- grep("*cov*", rownames(x$table.res), value = TRUE)
    colnames(random.effects.Predictions) <- c("id",unique(unique(gsub("\\*.*", "", gsub("__", "_", name_b)))))
  }
  else{
    random.effects.Predictions <- Objectranef
  }
  if("Y" %in% which){
    cv.Pred <- as.data.frame(cv.Pred)
    colnames(cv.Pred) <- c("id", "time", "predY", "predSD_inter", "predSD_intra")
  }
  else{
    cv.Pred <- NULL
  }


  resultat <- list(predictRE = random.effects.Predictions, predictY = cv.Pred)
  class(resultat) <- c("predict.lsmm_interintra")
  resultat









}
