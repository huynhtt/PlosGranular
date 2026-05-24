

# =========================
# 0. Packages
# =========================
library("readxl")
library("tidyr")
library("dplyr")
library("lubridate")
library("zoo")
library("forecast")
library("fracdiff")
library("tseries")
library("urca")
library("FinTS")
library("moments")
library("lmtest")


# =========================
# 1. Load data
# =========================
file_path <- "C:/Users/trong/OneDrive/HUYNHTT/Thesis 2025/Paper/Plos One/Granular Fractional Model/Code Data/data_FX.xlsx"

fx  <- read_excel(file_path, sheet = "fx")
vix <- read_excel(file_path, sheet = "vix")
wti <- read_excel(file_path, sheet = "wti")

# =========================
# 2. Convert yyyymmdd integer to Date
# =========================
convert_yyyymmdd <- function(x) {
  as.Date(as.character(as.integer(x)), format = "%Y%m%d")
}

fx <- fx %>%
  mutate(date = convert_yyyymmdd(date)) %>%
  arrange(date)

vix <- vix %>%
  mutate(date = convert_yyyymmdd(date)) %>%
  arrange(date)

wti <- wti %>%
  mutate(date = convert_yyyymmdd(date)) %>%
  arrange(date)

# =========================
# 3. Rename variables
# =========================
fx <- fx %>%
  rename(
    close_fx = close,
    high_fx  = high,
    low_fx   = low
  )

vix <- vix %>%
  rename(close_vix = close)

wti <- wti %>%
  rename(close_wti = close)

# =========================
# 4. Merge daily data
# =========================
daily_data <- fx %>%
  left_join(vix, by = "date") %>%
  left_join(wti, by = "date") %>%
  arrange(date)

# Fill missing VIX / WTI caused by different trading calendars
daily_data <- daily_data %>%
  mutate(
    close_vix = na.locf(close_vix, na.rm = FALSE),
    close_wti = na.locf(close_wti, na.rm = FALSE)
  ) %>%
  drop_na()

# =========================
# 5. Daily returns and ranges
# =========================
daily_data <- daily_data %>%
  mutate(
    r_fx_daily     = log(close_fx) - lag(log(close_fx)),
    range_fx_daily = log(high_fx) - log(low_fx),
    r_wti_daily    = log(close_wti) - lag(log(close_wti)),
    d_vix_daily    = close_vix - lag(close_vix)
  ) %>%
  drop_na()

# =========================
# 6. Aggregate to monthly data
# =========================
monthly_data <- daily_data %>%
  mutate(month = as.yearmon(date)) %>%
  group_by(month) %>%
  summarise(
    date = max(date),
    close_fx = last(close_fx),
    high_fx = max(high_fx, na.rm = TRUE),
    low_fx = min(low_fx, na.rm = TRUE),
    close_wti = last(close_wti),
    close_vix = last(close_vix),
    range_fx = log(max(high_fx, na.rm = TRUE)) - log(min(low_fx, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(date) %>%
  mutate(
    r_fx  = log(close_fx) - lag(log(close_fx)),
    r_wti = log(close_wti) - lag(log(close_wti)),
    d_vix = close_vix - lag(close_vix)
  ) %>%
  drop_na()

# =========================
# 7. Descriptive statistics
# =========================
desc_stats <- monthly_data %>%
  summarise(
    across(
      c(r_fx, range_fx, r_wti, d_vix),
      list(
        n = ~sum(!is.na(.)),
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        min = ~min(., na.rm = TRUE),
        max = ~max(., na.rm = TRUE),
        skew = ~skewness(., na.rm = TRUE),
        kurtosis = ~kurtosis(., na.rm = TRUE) - 3
      ),
      .names = "{.col}_{.fn}"
    )
  )

print(desc_stats)

# =========================
# 8. Correlation matrix with p-values
# =========================
vars <- monthly_data %>%
  select(r_fx, range_fx, r_wti, d_vix)

cor_matrix <- cor(vars, use = "complete.obs")

p_matrix <- matrix(NA, ncol = ncol(vars), nrow = ncol(vars))
colnames(p_matrix) <- colnames(vars)
rownames(p_matrix) <- colnames(vars)

for (i in 1:ncol(vars)) {
  for (j in 1:ncol(vars)) {
    p_matrix[i, j] <- cor.test(vars[[i]], vars[[j]])$p.value
  }
}

print(cor_matrix)
print(p_matrix)

# =========================
# 9. Unit root tests
# ADF, PP, KPSS
# =========================
unit_root_tests <- function(x, name) {
  x <- na.omit(as.numeric(x))
  
  adf_res  <- adf.test(x)
  pp_res   <- pp.test(x)
  kpss_res <- kpss.test(x, null = "Level")
  
  data.frame(
    variable = name,
    ADF_stat = as.numeric(adf_res$statistic),
    ADF_p = adf_res$p.value,
    PP_stat = as.numeric(pp_res$statistic),
    PP_p = pp_res$p.value,
    KPSS_stat = as.numeric(kpss_res$statistic),
    KPSS_p = kpss_res$p.value
  )
}

unit_root_table <- bind_rows(
  unit_root_tests(monthly_data$r_fx, "r_fx"),
  unit_root_tests(monthly_data$range_fx, "range_fx"),
  unit_root_tests(monthly_data$r_wti, "r_wti"),
  unit_root_tests(monthly_data$d_vix, "d_vix")
)

print(unit_root_table)

# =========================
# 10. BDS nonlinearity test
# =========================
bds_table <- list(
  r_fx = bds.test(na.omit(monthly_data$r_fx)),
  range_fx = bds.test(na.omit(monthly_data$range_fx)),
  r_wti = bds.test(na.omit(monthly_data$r_wti)),
  d_vix = bds.test(na.omit(monthly_data$d_vix))
)

print(bds_table)

# =========================
# 11. Fuzzification function
# =========================
make_fuzzy <- function(data, a = 1, xi = 0.5) {
  data %>%
    mutate(
      v_lower  = r_fx - (1 - xi) * a * range_fx,
      v_center = r_fx,
      v_upper  = r_fx + (1 - xi) * a * range_fx
    )
}

# =========================
# 12. Estimate full-sample fuzzy AR(2)
# with WTI and VIX controls
# =========================
estimate_far2 <- function(data, branch) {
  
  y <- data[[branch]]
  
  df <- data.frame(
    y = y,
    lag1 = lag(y, 1),
    lag2 = lag(y, 2),
    r_wti = data$r_wti,
    d_vix = data$d_vix
  ) %>%
    drop_na()
  
  model <- lm(y ~ lag1 + lag2 + r_wti + d_vix, data = df)
  return(model)
}

fuzzy_data <- make_fuzzy(monthly_data, a = 1, xi = 0.5)

model_lower  <- estimate_far2(fuzzy_data, "v_lower")
model_center <- estimate_far2(fuzzy_data, "v_center")
model_upper  <- estimate_far2(fuzzy_data, "v_upper")

summary(model_lower)
summary(model_center)
summary(model_upper)

far2_param_table <- data.frame(
  Branch = c("Lower", "Center", "Upper"),
  Intercept = c(coef(model_lower)[1], coef(model_center)[1], coef(model_upper)[1]),
  Phi1 = c(coef(model_lower)[2], coef(model_center)[2], coef(model_upper)[2]),
  Phi2 = c(coef(model_lower)[3], coef(model_center)[3], coef(model_upper)[3]),
  gamma_wti = c(coef(model_lower)[4], coef(model_center)[4], coef(model_upper)[4]),
  gamma_vix = c(coef(model_lower)[5], coef(model_center)[5], coef(model_upper)[5]),
  R2 = c(
    summary(model_lower)$r.squared,
    summary(model_center)$r.squared,
    summary(model_upper)$r.squared
  )
)

print(far2_param_table)

# =========================
# 13. Rolling FAR(2) forecast
# =========================
rolling_far2_forecast <- function(data, a = 1, xi = 0.5, window = 60) {
  
  df <- make_fuzzy(data, a = a, xi = xi) %>%
    mutate(
      v_lower_lag1  = lag(v_lower, 1),
      v_lower_lag2  = lag(v_lower, 2),
      v_center_lag1 = lag(v_center, 1),
      v_center_lag2 = lag(v_center, 2),
      v_upper_lag1  = lag(v_upper, 1),
      v_upper_lag2  = lag(v_upper, 2)
    ) %>%
    drop_na()
  
  n <- nrow(df)
  
  out <- data.frame(
    date = df$date,
    actual = df$r_fx,
    pred_lower = NA,
    pred_center = NA,
    pred_upper = NA,
    pred_rw = NA,
    pred_arima = NA,
    pred_arfima = NA
  )
  
  for (i in (window + 1):(n - 1)) {
    
    train <- df[(i - window + 1):i, ]
    test  <- df[i + 1, ]
    
    fit_l <- lm(
      v_lower ~ v_lower_lag1 + v_lower_lag2 + r_wti + d_vix,
      data = train
    )
    
    fit_c <- lm(
      v_center ~ v_center_lag1 + v_center_lag2 + r_wti + d_vix,
      data = train
    )
    
    fit_u <- lm(
      v_upper ~ v_upper_lag1 + v_upper_lag2 + r_wti + d_vix,
      data = train
    )
    
    out$pred_lower[i + 1] <- predict(
      fit_l,
      newdata = data.frame(
        v_lower_lag1 = test$v_lower_lag1,
        v_lower_lag2 = test$v_lower_lag2,
        r_wti = test$r_wti,
        d_vix = test$d_vix
      )
    )
    
    out$pred_center[i + 1] <- predict(
      fit_c,
      newdata = data.frame(
        v_center_lag1 = test$v_center_lag1,
        v_center_lag2 = test$v_center_lag2,
        r_wti = test$r_wti,
        d_vix = test$d_vix
      )
    )
    
    out$pred_upper[i + 1] <- predict(
      fit_u,
      newdata = data.frame(
        v_upper_lag1 = test$v_upper_lag1,
        v_upper_lag2 = test$v_upper_lag2,
        r_wti = test$r_wti,
        d_vix = test$d_vix
      )
    )
    
    # Lagged-return random walk
    out$pred_rw[i + 1] <- df$r_fx[i]
    
    # ARIMA benchmark
    fit_arima <- auto.arima(
      train$v_center,
      seasonal = FALSE,
      stepwise = FALSE,
      approximation = FALSE
    )
    
    out$pred_arima[i + 1] <- as.numeric(
      forecast(fit_arima, h = 1)$mean
    )
    
    # ARFIMA benchmark
    out$pred_arfima[i + 1] <- tryCatch({
      fit_arfima <- fracdiff(train$v_center, nar = 2, nma = 0)
      ar_coef <- fit_arfima$ar
      
      if (length(ar_coef) >= 2) {
        ar_coef[1] * df$v_center[i] + ar_coef[2] * df$v_center[i - 1]
      } else if (length(ar_coef) == 1) {
        ar_coef[1] * df$v_center[i]
      } else {
        mean(train$v_center, na.rm = TRUE)
      }
    }, error = function(e) {
      mean(train$v_center, na.rm = TRUE)
    })
  }
  
  out <- out %>% drop_na()
  return(out)
}

forecast_result <- rolling_far2_forecast(
  monthly_data,
  a = 1,
  xi = 0.5,
  window = 60
)
# =========================
# 14. Evaluation metrics
# =========================
rmse <- function(actual, pred) {
  sqrt(mean((actual - pred)^2, na.rm = TRUE))
}

mae <- function(actual, pred) {
  mean(abs(actual - pred), na.rm = TRUE)
}

evaluate_models <- function(fc) {
  
  data.frame(
    Model = c("Fuzzy AR(2)", "Random Walk", "ARIMA", "ARFIMA"),
    RMSE = c(
      rmse(fc$actual, fc$pred_center),
      rmse(fc$actual, fc$pred_rw),
      rmse(fc$actual, fc$pred_arima),
      rmse(fc$actual, fc$pred_arfima)
    ),
    MAE = c(
      mae(fc$actual, fc$pred_center),
      mae(fc$actual, fc$pred_rw),
      mae(fc$actual, fc$pred_arima),
      mae(fc$actual, fc$pred_arfima)
    )
  )
}

performance_table <- evaluate_models(forecast_result)
print(performance_table)

# =========================
# 15. Diebold-Mariano tests
# =========================
dm_far2_rw <- dm.test(
  e1 = forecast_result$actual - forecast_result$pred_center,
  e2 = forecast_result$actual - forecast_result$pred_rw,
  alternative = "less",
  h = 1,
  power = 2
)

dm_far2_arima <- dm.test(
  e1 = forecast_result$actual - forecast_result$pred_center,
  e2 = forecast_result$actual - forecast_result$pred_arima,
  alternative = "two.sided",
  h = 1,
  power = 2
)

dm_far2_arfima <- dm.test(
  e1 = forecast_result$actual - forecast_result$pred_center,
  e2 = forecast_result$actual - forecast_result$pred_arfima,
  alternative = "two.sided",
  h = 1,
  power = 2
)

dm_table <- data.frame(
  Comparison = c("FAR2 vs RW", "FAR2 vs ARIMA", "FAR2 vs ARFIMA"),
  DM_stat = c(
    as.numeric(dm_far2_rw$statistic),
    as.numeric(dm_far2_arima$statistic),
    as.numeric(dm_far2_arfima$statistic)
  ),
  p_value = c(
    dm_far2_rw$p.value,
    dm_far2_arima$p.value,
    dm_far2_arfima$p.value
  )
)

print(dm_table)

# =========================
# 16. Fuzzy interval coverage
# =========================
coverage <- mean(
  forecast_result$actual >= forecast_result$pred_lower &
    forecast_result$actual <= forecast_result$pred_upper,
  na.rm = TRUE
)

avg_width <- mean(
  forecast_result$pred_upper - forecast_result$pred_lower,
  na.rm = TRUE
)

interval_result <- data.frame(
  Coverage = coverage,
  Average_width = avg_width
)

print(interval_result)

# =========================
# 17. Sensitivity analysis for fuzziness parameter a
# =========================
a_grid <- c(0.25, 0.50, 0.75, 1.00, 1.50, 2.00)

sensitivity_table <- lapply(a_grid, function(a_value) {
  
  fc <- rolling_far2_forecast(
    monthly_data,
    a = a_value,
    xi = 0.5,
    window = 60
  )
  
  dm_res <- dm.test(
    e1 = fc$actual - fc$pred_center,
    e2 = fc$actual - fc$pred_rw,
    alternative = "less",
    h = 1,
    power = 2
  )
  
  data.frame(
    a = a_value,
    RMSE = rmse(fc$actual, fc$pred_center),
    MAE = mae(fc$actual, fc$pred_center),
    Coverage = mean(
      fc$actual >= fc$pred_lower &
        fc$actual <= fc$pred_upper,
      na.rm = TRUE
    ),
    Band_width = mean(fc$pred_upper - fc$pred_lower, na.rm = TRUE),
    DM_vs_RW = as.numeric(dm_res$statistic),
    p_value = dm_res$p.value
  )
  
}) %>%
  bind_rows()

print(sensitivity_table)

# =========================
# 18. Export tables
# =========================
#write.csv(desc_stats, "summary_statistics.csv", row.names = FALSE)
#write.csv(cor_matrix, "correlation_matrix.csv")
#write.csv(p_matrix, "correlation_pvalues.csv")
#write.csv(unit_root_table, "unit_root_tests.csv", row.names = FALSE)
#write.csv(far2_param_table, "far2_parameters.csv", row.names = FALSE)
#write.csv(performance_table, "forecast_performance.csv", row.names = FALSE)
#write.csv(dm_table, "dm_tests.csv", row.names = FALSE)
#write.csv(sensitivity_table, "sensitivity_fuzziness_a.csv", row.names = FALSE)
#write.csv(forecast_result, "rolling_forecast_results.csv", row.names = FALSE)








############################################################
# 19. HELPER FUNCTION FOR CLEAN PRINTING
############################################################

round_df <- function(x, digits = 4) {
  x %>%
    mutate(
      across(where(is.numeric), ~round(., digits))
    )
}

############################################################
# 20. SUMMARY STATISTICS TABLE
############################################################

summary_table <- monthly_data %>%
  summarise(
    across(
      c(r_fx, range_fx, r_wti, d_vix),
      list(
        n = ~sum(!is.na(.)),
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        min = ~min(., na.rm = TRUE),
        max = ~max(., na.rm = TRUE),
        skew = ~moments::skewness(., na.rm = TRUE),
        kurtosis = ~moments::kurtosis(., na.rm = TRUE) - 3
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "(.+)_(n|mean|sd|min|max|skew|kurtosis)"
  )

############################################################
# 21. CORRELATION MATRIX WITH P-VALUES
############################################################

cor_with_p <- matrix(
  "",
  nrow = nrow(cor_matrix),
  ncol = ncol(cor_matrix)
)

rownames(cor_with_p) <- rownames(cor_matrix)
colnames(cor_with_p) <- colnames(cor_matrix)

for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    
    cor_with_p[i, j] <- paste0(
      round(cor_matrix[i, j], 4),
      " (",
      round(p_matrix[i, j], 4),
      ")"
    )
  }
}

############################################################
# 22. COMPACT BDS TABLE
############################################################

extract_bds <- function(bds_obj, varname) {
  
  stat <- as.data.frame(as.table(bds_obj$statistic))
  pval <- as.data.frame(as.table(bds_obj$p.value))
  
  names(stat) <- c("Dimension", "Epsilon", "BDS_stat")
  names(pval) <- c("Dimension", "Epsilon", "p_value")
  
  left_join(stat, pval,
            by = c("Dimension", "Epsilon")) %>%
    mutate(
      Variable = varname,
      BDS_stat = round(as.numeric(BDS_stat), 4),
      p_value = round(as.numeric(p_value), 4)
    ) %>%
    select(
      Variable,
      Dimension,
      Epsilon,
      BDS_stat,
      p_value
    )
}

bds_compact_table <- bind_rows(
  extract_bds(bds_table$r_fx, "r_fx"),
  extract_bds(bds_table$range_fx, "range_fx"),
  extract_bds(bds_table$r_wti, "r_wti"),
  extract_bds(bds_table$d_vix, "d_vix")
)

############################################################
# 23. FORECAST SAMPLE TABLE
############################################################

forecast_preview <- forecast_result %>%
  select(
    date,
    actual,
    pred_lower,
    pred_center,
    pred_upper,
    pred_rw,
    pred_arima,
    pred_arfima
  )

############################################################
# 24. MODEL RANKING TABLES
############################################################

rmse_ranking <- performance_table %>%
  arrange(RMSE) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, everything())

mae_ranking <- performance_table %>%
  arrange(MAE) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, everything())

############################################################
# 25. SHOW RESULTS
############################################################

cat("\n==================================================\n")
cat("DATA OVERVIEW\n")
cat("==================================================\n")

cat("\nShow number of observations:\n")
print(nrow(monthly_data))

cat("\nShow sample period:\n")
print(range(monthly_data$date))

cat("\nShow variables:\n")
print(names(monthly_data))


cat("\n==================================================\n")
cat("SUMMARY STATISTICS\n")
cat("==================================================\n")

cat("\nShow descriptive statistics:\n")
print(round_df(summary_table, 4))


cat("\n==================================================\n")
cat("CORRELATION MATRIX\n")
cat("==================================================\n")

cat("\nShow correlation matrix:\n")
print(round(cor_matrix, 4))

cat("\nShow p-values:\n")
print(round(p_matrix, 4))

cat("\nShow formatted correlation matrix:\n")
print(cor_with_p, quote = FALSE)


cat("\n==================================================\n")
cat("UNIT ROOT TESTS\n")
cat("==================================================\n")

cat("\nShow ADF / PP / KPSS results:\n")
print(round_df(unit_root_table, 4))

cat("\nNote:\n")
cat("- Small ADF/PP p-values => reject unit root.\n")
cat("- Large KPSS p-values => stationarity supported.\n")


cat("\n==================================================\n")
cat("BDS NONLINEARITY TESTS\n")
cat("==================================================\n")

cat("\nShow compact BDS table:\n")
print(bds_compact_table)

cat("\nNote:\n")
cat("- Small BDS p-values => nonlinear dependence exists.\n")
cat("- Supports fuzzy / nonlinear modeling.\n")


cat("\n==================================================\n")
cat("FULL-SAMPLE FUZZY AR(2)\n")
cat("==================================================\n")

cat("\nShow estimated parameters:\n")
print(round_df(far2_param_table, 6))

cat("\nShow detailed regression summaries:\n")

cat("\nLOWER BRANCH\n")
print(summary(model_lower))

cat("\nCENTER BRANCH\n")
print(summary(model_center))

cat("\nUPPER BRANCH\n")
print(summary(model_upper))


cat("\n==================================================\n")
cat("ROLLING FORECAST PERFORMANCE\n")
cat("==================================================\n")

cat("\nShow RMSE and MAE:\n")
print(round_df(performance_table, 4))

cat("\nShow Diebold-Mariano tests:\n")
print(round_df(dm_table, 4))

cat("\nNote:\n")
cat("- Lower RMSE/MAE => better forecast accuracy.\n")
cat("- Negative DM statistic => FAR(2) outperforms benchmark.\n")


cat("\n==================================================\n")
cat("FUZZY INTERVAL RESULTS\n")
cat("==================================================\n")

cat("\nShow interval coverage and width:\n")
print(round_df(interval_result, 4))

cat("\nNote:\n")
cat("- Coverage measures proportion inside fuzzy interval.\n")
cat("- Width measures uncertainty band size.\n")


cat("\n==================================================\n")
cat("SENSITIVITY ANALYSIS\n")
cat("==================================================\n")

cat("\nShow sensitivity results for parameter a:\n")
print(round_df(sensitivity_table, 4))

cat("\nNote:\n")
cat("- Stable RMSE/MAE across a => point forecasts robust.\n")
cat("- Increasing coverage/width => uncertainty representation controlled by a.\n")


cat("\n==================================================\n")
cat("FORECAST SAMPLE CHECK\n")
cat("==================================================\n")

cat("\nShow number of forecasts:\n")
print(nrow(forecast_result))

cat("\nShow forecast period:\n")
print(range(forecast_result$date))

cat("\nShow first 10 forecasts:\n")
print(
  forecast_preview %>%
    head(10) %>%
    round_df(5)
)

cat("\nShow last 10 forecasts:\n")
print(
  forecast_preview %>%
    tail(10) %>%
    round_df(5)
)


cat("\n==================================================\n")
cat("MODEL RANKING\n")
cat("==================================================\n")

cat("\nShow ranking by RMSE:\n")
print(round_df(rmse_ranking, 4))

cat("\nShow ranking by MAE:\n")
print(round_df(mae_ranking, 4))


cat("\n==================================================\n")
cat("FINAL NOTES\n")
cat("==================================================\n")

cat("\n1. Stationary returns + significant BDS tests support nonlinear modeling.\n")
cat("2. FAR(2) should outperform RW if conditional predictability exists.\n")
cat("3. Similar FAR(2), ARIMA, ARFIMA performance implies competitive rather than dominant forecasting power.\n")
cat("4. Fuzzy intervals provide explicit uncertainty representation.\n")
cat("5. Sensitivity analysis evaluates robustness of fuzziness parameter a.\n")

############################################################
# END
############################################################
############################################################
# END
############################################################