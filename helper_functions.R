simulate_data <- function(age, beta, disp = 1, family) {
  X <- model.matrix(~ age)
  n <- length(age)
  mu <- exp(X %*% beta)
  y <- switch(family,
              "poisson" = rpois(n, mu),
              "nb1" = rNBII(n, mu, disp),
              "nb2" = rNBI(n, mu, disp),
              "dp" = rDPO(n, mu, disp))
  sim_data <- data.frame(y = y, age = age, obs = factor(1:n))
  return(sim_data)
}

fit_model_suite <- function(data) {
  fit1 <- glm(y ~ age, family = poisson, data = data) 
  fit2 <- glm(y ~ age, family = quasipoisson, data = data)
  fit3 <- glmmTMB(y ~ age, family = nbinom1, data = data)
  fit4 <- glmmTMB(y ~ age, family = nbinom2, data = data)
  fit5 <- glmmTMB(y ~ age, family = compois, data = data)
  fit6 <- glmmTMB(y ~ age, family = genpois, data = data)
  fit7 <- glmmTMB(y ~ age + (1 | obs), family = poisson, data = data)
  return(list("poisson" = fit1,
              "qp" = fit2,
              "nb1" = fit3,
              "nb2" = fit4,
              "comp" = fit5,
              "gp" = fit6,
              "pn" = fit7))
}

n_out_hnp <- function(obj, family) {
  hnp_obj <- adapted_hnp(obj, family = family)
  perc_out <- hnp_obj$out/hnp_obj$total * 100
  return(perc_out)
}

adapted_hnp <- function(obj, family) {
  if(family %in% c("poisson","quasipoisson"))  {
    return(hnp(obj, resid.type = "pearson", 
               how.many.out = TRUE, plot.sim = FALSE))
  } else {
    obj_data <- obj$frame
    dfun <- function(obj) resid(obj, type = "pearson")
    sfun <- function(n, obj) simulate(obj)$sim_1
    if(family != "poisson_normal") {
      ffun <- function(response) {
        fit <- try(glmmTMB(response ~ age, family = family, data = obj_data), silent = TRUE)
        while(class(fit) == "try-error") {
          response2 <- sfun(1, obj)
          fit <- try(glmmTMB(response2 ~ age, family = family, data = obj_data), silent = TRUE)
        }
        return(fit)
      }
    } else {
      ffun <- function(response) {
        fit <- try(glmmTMB(response ~ age + (1 | obs), family = poisson, data = obj_data), silent = TRUE)
        while(class(fit) == "try-error") {
          response2 <- sfun(1, obj)
          fit <- try(glmmTMB(response2 ~ age + (1 | obs), family = poisson, data = obj_data), silent = TRUE)
        }
        return(fit)
      }
    }
    return(hnp(obj, newclass = TRUE,
               diagfun = dfun, simfun = sfun, fitfun = ffun,
               how.many.out = TRUE, plot.sim = FALSE))
  }
}

get_family <- function(obj) {
  if(class(obj)[1] == "glm") {
    return(paste(obj$family$family))
  } else {
    fam <- paste(obj$call$family)
    if(fam == "poisson") fam <- "poisson_normal"
    return(fam)
  }
}

get_coef <- function(obj) {
  if(class(obj)[1] == "glm") {
    return(c(coef(obj), summary(obj)$dispersion))
  } else {
    return(as.numeric(obj$fit$par))
  }
}

AICc <- function(obj) {
  aic <- AIC(obj)
  if(class(obj)[1] == "glm") {
    np <- length(obj$coef)
    n <- nrow(obj$data)
    if(obj$family$family == "quasipoisson") {
      np <- np + 1
    }
  } else {
    np <- length(obj$fit$par)
    n <- nrow(obj$frame)
  }
  aicc <- aic + (2 * np^2 + 2 * np)/(n - np - 1)
  return(aicc)
}

run_study <- function(age, beta, disp, family, n_sim) {
  all_data <- lapply(1:n_sim,
                     function(x) simulate_data(age = age, beta = beta, disp = disp, family = family))
  all_models <- lapply(all_data,
                       function(x) fit_model_suite(x))
  all_coef <- lapply(all_models,
                     function(x) unlist(sapply(x,
                                               function(y) get_coef(y))))
  all_aic <- lapply(all_models,
                    function(x) unlist(lapply(x,
                                              function(y) AIC(y))))
  all_hnp <- lapply(all_models,
                    function(x) unlist(lapply(x,
                                              function(y) n_out_hnp(y, family = get_family(y)))))
  return(list("data" = all_data,
              "models" = all_models,
              "coef" = all_coef,
              "aic" = all_aic,
              "hnp" = all_hnp))
}