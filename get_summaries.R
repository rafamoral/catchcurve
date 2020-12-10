load("results_simulations.RData")

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

aic_summary1 <- do.call(data.frame, results1$aic)
colnames(aic_summary1) <- 1:n_sim
aic_summary1[aic_summary1 < 0] <- NA
aic1 <- apply(aic_summary1, 1, mean, na.rm = TRUE)

aic_summary2 <- do.call(data.frame, results2$aic)
colnames(aic_summary2) <- 1:n_sim
aic_summary2[aic_summary2 < 0] <- NA
aic2 <- apply(aic_summary2, 1, mean, na.rm = TRUE)

aic_summary3 <- do.call(data.frame, results3$aic)
colnames(aic_summary3) <- 1:n_sim
aic_summary3[aic_summary3 < 0] <- NA
aic3 <- apply(aic_summary3, 1, mean, na.rm = TRUE)

aic_summary4 <- do.call(data.frame, results4$aic)
colnames(aic_summary4) <- 1:n_sim
aic_summary4[aic_summary4 < 0] <- NA
aic4 <- apply(aic_summary4, 1, mean, na.rm = TRUE)

results1$aicc <- lapply(results1$models,
                        function(x) unlist(lapply(x,
                                                  function(y) AICc(y))))
results2$aicc <- lapply(results2$models,
                        function(x) unlist(lapply(x,
                                                  function(y) AICc(y))))
results3$aicc <- lapply(results3$models,
                        function(x) unlist(lapply(x,
                                                  function(y) AICc(y))))
results4$aicc <- lapply(results4$models,
                        function(x) unlist(lapply(x,
                                                  function(y) AICc(y))))

aicc_summary1 <- do.call(data.frame, results1$aicc)
colnames(aicc_summary1) <- 1:n_sim
aicc_summary1[aicc_summary1 < 0] <- NA
aicc1 <- apply(aicc_summary1, 1, mean, na.rm = TRUE)

aicc_summary2 <- do.call(data.frame, results2$aicc)
colnames(aicc_summary2) <- 1:n_sim
aicc_summary2[aicc_summary2 < 0] <- NA
aicc2 <- apply(aicc_summary2, 1, mean, na.rm = TRUE)

aicc_summary3 <- do.call(data.frame, results3$aicc)
colnames(aicc_summary3) <- 1:n_sim
aicc_summary3[aicc_summary3 < 0] <- NA
aicc3 <- apply(aicc_summary3, 1, mean, na.rm = TRUE)

aicc_summary4 <- do.call(data.frame, results4$aicc)
colnames(aicc_summary4) <- 1:n_sim
aicc_summary4[aicc_summary4 < 0] <- NA
aicc4 <- apply(aicc_summary4, 1, mean, na.rm = TRUE)

hnp_summary1 <- do.call(data.frame, results1$hnp)
colnames(hnp_summary1) <- 1:n_sim
na_index <- is.na(aic_summary1)
na_index[1:2,] <- FALSE
hnp_summary1[na_index] <- NA
hnp1 <- apply(hnp_summary1, 1, mean, na.rm = TRUE)

hnp_summary2 <- do.call(data.frame, results2$hnp)
colnames(hnp_summary2) <- 1:n_sim
na_index <- is.na(aic_summary2)
na_index[1:2,] <- FALSE
hnp_summary2[na_index] <- NA
hnp2 <- apply(hnp_summary2, 1, mean, na.rm = TRUE)

hnp_summary3 <- do.call(data.frame, results3$hnp)
colnames(hnp_summary3) <- 1:n_sim
na_index <- is.na(aic_summary3)
na_index[1:2,] <- FALSE
hnp_summary3[na_index] <- NA
hnp3 <- apply(hnp_summary3, 1, mean, na.rm = TRUE)

hnp_summary4 <- do.call(data.frame, results4$hnp)
colnames(hnp_summary4) <- 1:n_sim
na_index <- is.na(aic_summary4)
na_index[1:2,] <- FALSE
hnp_summary4[na_index] <- NA
hnp4 <- apply(hnp_summary4, 1, mean, na.rm = TRUE)

coef_all1 <- do.call(data.frame, lapply(results1$coef, function(x) t(x[-3,])))
coef_all2 <- do.call(data.frame, lapply(results2$coef, function(x) t(x[-3,])))
coef_all3 <- do.call(data.frame, lapply(results3$coef, function(x) t(x[-3,])))
coef_all4 <- do.call(data.frame, lapply(results4$coef, function(x) t(x[-3,])))

intercept1 <- coef_all1[,seq(1,39,2)]
intercept2 <- coef_all2[,seq(1,39,2)]
intercept3 <- coef_all3[,seq(1,39,2)]
intercept4 <- coef_all4[,seq(1,39,2)]
slope1 <- coef_all1[,seq(2,40,2)]
slope2 <- coef_all2[,seq(2,40,2)]
slope3 <- coef_all3[,seq(2,40,2)]
slope4 <- coef_all4[,seq(2,40,2)]

est_intercept1 <- apply(intercept1, 1, mean)
est_intercept2 <- apply(intercept2, 1, mean)
est_intercept3 <- apply(intercept3, 1, mean)
est_intercept4 <- apply(intercept4, 1, mean)
est_slope1 <- apply(slope1, 1, mean)
est_slope2 <- apply(slope2, 1, mean)
est_slope3 <- apply(slope3, 1, mean)
est_slope4 <- apply(slope4, 1, mean)

sd_intercept1 <- apply(intercept1, 1, sd)
sd_intercept2 <- apply(intercept2, 1, sd)
sd_intercept3 <- apply(intercept3, 1, sd)
sd_intercept4 <- apply(intercept4, 1, sd)
sd_slope1 <- apply(slope1, 1, sd)
sd_slope2 <- apply(slope2, 1, sd)
sd_slope3 <- apply(slope3, 1, sd)
sd_slope4 <- apply(slope4, 1, sd)

bias_intercept1 <- apply(intercept1 - true_beta[1], 1, mean)
bias_intercept2 <- apply(intercept2 - true_beta[1], 1, mean)
bias_intercept3 <- apply(intercept3 - true_beta[1], 1, mean)
bias_intercept4 <- apply(intercept4 - true_beta[1], 1, mean)
bias_slope1 <- apply(slope1 - true_beta[2], 1, mean)
bias_slope2 <- apply(slope2 - true_beta[2], 1, mean)
bias_slope3 <- apply(slope3 - true_beta[2], 1, mean)
bias_slope4 <- apply(slope4 - true_beta[2], 1, mean)

mse_intercept1 <- apply((intercept1 - true_beta[1])^2, 1, sum)
mse_intercept2 <- apply((intercept2 - true_beta[1])^2, 1, sum)
mse_intercept3 <- apply((intercept3 - true_beta[1])^2, 1, sum)
mse_intercept4 <- apply((intercept4 - true_beta[1])^2, 1, sum)
rmse_intercept1 <- sqrt(mse_intercept1)
rmse_intercept2 <- sqrt(mse_intercept2)
rmse_intercept3 <- sqrt(mse_intercept3)
rmse_intercept4 <- sqrt(mse_intercept4)
mse_slope1 <- apply((slope1 - true_beta[2])^2, 1, sum)
mse_slope2 <- apply((slope2 - true_beta[2])^2, 1, sum)
mse_slope3 <- apply((slope3 - true_beta[2])^2, 1, sum)
mse_slope4 <- apply((slope4 - true_beta[2])^2, 1, sum)
rmse_slope1 <- sqrt(mse_slope1)
rmse_slope2 <- sqrt(mse_slope2)
rmse_slope3 <- sqrt(mse_slope3)
rmse_slope4 <- sqrt(mse_slope4)

get_se <- function(obj) {
  if(class(obj)[1] == "glm") {
    return(summary(obj)$coefficients[,2])
  } else {
    return(summary(obj)$coefficients$cond[,2])
  }
}

se1 <- lapply(results1$models, function(x) unlist(lapply(x, get_se)))
se2 <- lapply(results2$models, function(x) unlist(lapply(x, get_se)))
se3 <- lapply(results3$models, function(x) unlist(lapply(x, get_se)))
se4 <- lapply(results4$models, function(x) unlist(lapply(x, get_se)))
se1_all <- do.call(cbind, se1)
se2_all <- do.call(cbind, se2)
se3_all <- do.call(cbind, se3)
se4_all <- do.call(cbind, se4)
se1_int <- se1_all[seq(1,13,2),]
se2_int <- se2_all[seq(1,13,2),]
se3_int <- se3_all[seq(1,13,2),]
se4_int <- se4_all[seq(1,13,2),]
se1_slo <- se1_all[seq(2,14,2),]
se2_slo <- se2_all[seq(2,14,2),]
se3_slo <- se3_all[seq(2,14,2),]
se4_slo <- se4_all[seq(2,14,2),]
se_intercept1 <- apply(se1_int, 1, mean)
se_intercept2 <- apply(se2_int, 1, mean)
se_intercept3 <- apply(se3_int, 1, mean)
se_intercept4 <- apply(se4_int, 1, mean)
se_slope1 <- apply(se1_slo, 1, mean)
se_slope2 <- apply(se2_slo, 1, mean)
se_slope3 <- apply(se3_slo, 1, mean)
se_slope4 <- apply(se4_slo, 1, mean)

aic_all <- cbind(aic1, aic2, aic3, aic4)
aicc_all <- cbind(aicc1, aicc2, aicc3, aicc4)
hnp_all <- cbind(hnp1, hnp2, hnp3, hnp4)
intercept_all <- cbind(est_intercept1, est_intercept2, est_intercept3, est_intercept4)
slope_all <- cbind(est_slope1, est_slope2, est_slope3, est_slope4)
sd_intercept_all <- cbind(sd_intercept1, sd_intercept2, sd_intercept3, sd_intercept4)
sd_slope_all <- cbind(sd_slope1, sd_slope2, sd_slope3, sd_slope4)
bias_intercept_all <- cbind(bias_intercept1, bias_intercept2, bias_intercept3, bias_intercept4)
bias_slope_all <- cbind(bias_slope1, bias_slope2, bias_slope3, bias_slope4)
mse_intercept_all <- cbind(mse_intercept1, mse_intercept2, mse_intercept3, mse_intercept4)
mse_slope_all <- cbind(mse_slope1, mse_slope2, mse_slope3, mse_slope4)
rmse_intercept_all <- cbind(rmse_intercept1, rmse_intercept2, rmse_intercept3, rmse_intercept4)
rmse_slope_all <- cbind(rmse_slope1, rmse_slope2, rmse_slope3, rmse_slope4)
se_intercept_all <- cbind(se_intercept1, se_intercept2, se_intercept3, se_intercept4)
se_slope_all <- cbind(se_slope1, se_slope2, se_slope3, se_slope4)

colnames(aic_all) <- colnames(aicc_all) <- colnames(hnp_all) <- 
  colnames(intercept_all) <- colnames(slope_all) <- colnames(bias_intercept_all) <-
  colnames(bias_slope_all) <- colnames(mse_intercept_all) <- colnames(mse_slope_all) <-
  colnames(rmse_intercept_all) <- colnames(rmse_slope_all) <- 
  colnames(se_intercept_all) <- colnames(se_slope_all) <- c("equidisp","linear_overdisp","quad_overdisp","underdisp")

rownames(se_intercept_all) <- rownames(se_slope_all) <- rownames(rmse_intercept_all)