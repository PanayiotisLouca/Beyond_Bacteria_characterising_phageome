## Author: 
# Panayiotis Louca 

## Purpose of script: 
#  

## Date Created: 
# 11 April 2025 

## Notes: 
#  

## Clear environment 
rm(list = ls()) 

## Set seed 
set.seed(1234)

## Set functions: 
#  

## load up packages: 

### core 
library(tidyverse)

# mixed effect modelling 
library(lmerTest)
library(glmmTMB)
library(broom.mixed)

# parallelisation 
library(parallel)
library(future.apply)

# -------------------------------------------------------------------------- # 

# ************************* # 
#   IMPORT & PREP DATA   ---- 
# ************************* # 

# -------------------------------------------------------------------------- #  

##   Dataset ---- 

path <- file.path("/scratch/users/k2480753/phages/phages_SCFA/phages_microbes_metabs_full_DATASET_SCFA_v2.csv")
df <- read.csv(path) %>% as.data.frame(.)

# vectorise columns 
phage_columns = grep("phage_ivig_|phage_uvig_", names(df), value = TRUE)

bile_acid_names = 
  c( # serum secondary bile acid id's 
    'serum_metab_M12261', 
    'serum_metab_M18477',
    "serum_metab_M31912",
    'serum_metab_M32620',
    'serum_metab_M36850',
    "serum_metab_M32599",
    "serum_metab_M32807",
    "serum_metab_M57577",
    "serum_metab_M39379",
    "serum_metab_M42574",
    "serum_metab_M52975",
    "serum_metab_M63610",
    "serum_metab_M63731",
    "serum_metab_M62526",
    "serum_metab_M62784",
    "serum_metab_M63688",
    # stool secondary bile acid id's 
    "stool_metab_M1114",
    "stool_metab_M1483",
    "stool_metab_M12261",
    "stool_metab_M18477",
    "stool_metab_M43835",
    "stool_metab_M31904",
    "stool_metab_M31912",
    "stool_metab_M31889",
    "stool_metab_M31891",
    "stool_metab_M46336",
    "stool_metab_M43830",
    "stool_metab_M44644",
    "stool_metab_M34093",
    "stool_metab_M32620",
    "stool_metab_M36850",
    "stool_metab_M32599",
    "stool_metab_M32807",
    "stool_metab_M57577",
    "stool_metab_M39379",
    "stool_metab_M39378",
    "stool_metab_M57573",
    "stool_metab_M43501",
    "stool_metab_M52969",
    "stool_metab_M52975",
    "stool_metab_M52970",
    "stool_metab_M57774",
    "stool_metab_M62526",
    "stool_metab_M62527",
    "stool_metab_M63118",
    "stool_metab_M63119",
    "stool_metab_M62890",
    "stool_metab_M63688"
  )

ba_columns = bile_acid_names[which(bile_acid_names %in% names(df))]

table(df$gm_batch)

covar_columns = c('gm_age','gm_sex','gm_BMI','gm_batch')

# -------------------------------------------------------------------------- #  

# clean variable types 

# phages 
df %>% select(all_of(phage_columns)) %>%
  sapply(., class) %>%
  table


# bile acids 
df %>% select(all_of(ba_columns)) %>%
  sapply(., class) %>%
  table


# covars 
df %>% select(all_of(covar_columns)) %>%
  glimpse

df <- df %>%
  mutate(gm_sex = factor(gm_sex),
         fid = factor(fid),
         gm_batch = factor(gm_batch))

# -------------------------------------------------------------------------- #  

df_analyse <- df %>%
  select(
    iid,
    all_of(covar_columns),
    fid,
    # phages 
    all_of(phage_columns),
    # Bile acid columns  
    all_of(ba_columns)
  )

# -------------------------------------------------------------------------- #  

# Scale/normalise 
#- qunatile normalisation function 
rankit = function(vect) {
  # Nan position 
  vect.position = is.na(vect)
  # Rank non Nans 
  qnormvect = qnorm((rank(vect[!vect.position]) - 1/2) / length(vect[!vect.position]))
  # Replace non-Nans with rank transformed values 
  vect[!vect.position] = qnormvect
  return(vect)
}

# z score age and BMI 
df_analyse$gm_age = (df_analyse$gm_age - mean(df_analyse$gm_age))/sd(df_analyse$gm_age)
df_analyse$gm_BMI = (df_analyse$gm_BMI - mean(df_analyse$gm_BMI))/sd(df_analyse$gm_BMI)

# CLR normalisation function for phages 
CLRnorm <- function(features) {
  if (!is.data.frame(features) & !is.matrix(features)) {
    stop("Input must be a data frame or matrix.")
  }
  features_norm <- as.matrix(features)
  features_norm[features_norm == 0] <- 1e-6
  features_CLR <- chemometrics::clr(features_norm)
  as.data.frame(features_CLR, col.names = colnames(features))
}

# -------------------------------------------------------------------------- #  

###   setup data ---- 

# gut phage data only 
dat_phage <- df_analyse %>%
  select(iid,
         all_of(phage_columns)
  )

# convert phage data to decimal proportions 
dat_phage[ ,-1] <- sapply(dat_phage[ ,-1], function(x) (x/100))

dat_phage <- # move iid to rownames 
  dat_phage %>%
  tibble::column_to_rownames('iid')

# Apply CLR transformation to phage columns  
dat_phage <- CLRnorm(dat_phage)

# Apply inverse normalisation to phage columns 
dat_phage <- dat_phage %>%
  mutate(across(everything(), rankit))

# Bile acid data only 
dat_metabs <- df_analyse %>%
  select(iid,
         all_of(ba_columns))

dat_metabs <- dat_metabs %>%
  tibble::column_to_rownames('iid')

# z score metabolites 
dat_metabs <- lapply(dat_metabs, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))

# -------------------------------------------------------------------------- #  

# Function for error and warning reporting 
# REF: https://stackoverflow.com/a/4952908/10360530 
catch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(v=value, w=warn, e=err)
}

# -------------------------------------------------------------------------- #  

set.seed(1234)

######   CREATE LMER FUNCTION  ---- 

# Create formulas for unadjusted and adjusted models 
formula_unadjusted <- "y ~ x"
formula_adjusted <- "y ~ x + gm_age + gm_sex + gm_BMI + (1 | fid) + (1 | gm_batch)"

# REGRESSION FUNCTION 
run_regression_parallel <- function(dependent_vars,
                                    independ_vars,
                                    covars = NULL,
                                    cores = parallel::detectCores() - 1) {
  start_time <- Sys.time()
  
  dependent_vars <- as.data.frame(dependent_vars)
  independ_vars <- as.data.frame(independ_vars)
  
  # Setup parallel backend 
  future::plan("multisession", workers = cores)
  progressr::handlers(global = TRUE)
  
  # All (i,j) combinations 
  combos <- expand.grid(
    depend_idx = seq_along(dependent_vars),
    independ_idx = seq_along(independ_vars)
  )
  
  results_list <- progressr::with_progress({
    p <- progressr::progressor(along = 1:nrow(combos))
    future_lapply(1:nrow(combos), function(k) {
      i <- combos$depend_idx[k]
      j <- combos$independ_idx[k]
      p()  # update progress
      
      depend_var_name <- colnames(dependent_vars)[i]
      independ_var_name <- colnames(independ_vars)[j]
      y <- dependent_vars[[i]]
      x <- independ_vars[[j]]
      
      df_model <- covars
      df_model$y <- y
      df_model$x <- x
      
      n_total <- nrow(df_model)
      n_not_0_depend_var <- sum(y != 0, na.rm = TRUE)
      n_not_0_independ_var <- sum(x != 0, na.rm = TRUE)
      
      # Skip if too few observations 
      if (n_total < 10) {
        return(data.frame(
          depend_var = depend_var_name,
          independ_var = independ_var_name,
          estimate_unadjusted = NA,
          std_err_unadjusted = NA,
          ci_low_unadjusted = NA,
          ci_upper_unadjusted = NA,
          pval_unadjusted = NA,
          estimate_adjusted = NA,
          std_err_adjusted = NA,
          ci_low_adjusted = NA,
          ci_upper_adjusted = NA,
          pval_adjusted = NA,
          formula_unadjusted = formula_unadjusted,
          formula_adjusted = formula_adjusted,
          n_total = n_total,
          n_not_0_depend_var = n_not_0_depend_var,
          n_not_0_independ_var = n_not_0_independ_var,
          error_message_unadjusted = "Too few observations",
          error_message_adjusted = "Too few observations",
          warning_message_unadjusted = NA,
          warning_message_adjusted = NA,
          stringsAsFactors = FALSE
        ))
      }
      
      # UNADJUSTED 
      unadj_result <- tryCatch({
        model_unadj <- glmmTMB(as.formula(formula_unadjusted),
                               data = df_model, REML = FALSE)
        model_output_unadj <- broom.mixed::tidy(model_unadj, conf.int = TRUE) %>%
          dplyr::filter(term == "x") %>%
          dplyr::select(estimate, std.error, conf.low, conf.high, p.value) %>%
          dplyr::rename(
            estimate_unadjusted = estimate,
            std_err_unadjusted = std.error,
            ci_low_unadjusted = conf.low,
            ci_upper_unadjusted = conf.high,
            pval_unadjusted = p.value
          )
        if (nrow(model_output_unadj) == 0)
          stop("No unadjusted model output for variable")
        list(output = model_output_unadj, error = NA, warning = NA)
      }, warning = function(w) {
        list(output = data.frame(
          estimate_unadjusted = NA,
          std_err_unadjusted = NA,
          ci_low_unadjusted = NA,
          ci_upper_unadjusted = NA,
          pval_unadjusted = NA
        ), error = NA, warning = w$message)
      }, error = function(e) {
        list(output = data.frame(
          estimate_unadjusted = NA,
          std_err_unadjusted = NA,
          ci_low_unadjusted = NA,
          ci_upper_unadjusted = NA,
          pval_unadjusted = NA
        ), error = e$message, warning = NA)
      })
      
      # ADJUSTED 
      adj_result <- tryCatch({
        model_adj <- glmmTMB(as.formula(formula_adjusted),
                             data = df_model, REML = FALSE)
        model_output_adj <- broom.mixed::tidy(model_adj, conf.int = TRUE) %>%
          dplyr::filter(term == "x") %>%
          dplyr::select(estimate, std.error, conf.low, conf.high, p.value) %>%
          dplyr::rename(
            estimate_adjusted = estimate,
            std_err_adjusted = std.error,
            ci_low_adjusted = conf.low,
            ci_upper_adjusted = conf.high,
            pval_adjusted = p.value
          )
        if (nrow(model_output_adj) == 0)
          stop("No adjusted model output for variable")
        list(output = model_output_adj, error = NA, warning = NA)
      }, warning = function(w) {
        list(output = data.frame(
          estimate_adjusted = NA,
          std_err_adjusted = NA,
          ci_low_adjusted = NA,
          ci_upper_adjusted = NA,
          pval_adjusted = NA
        ), error = NA, warning = w$message)
      }, error = function(e) {
        list(output = data.frame(
          estimate_adjusted = NA,
          std_err_adjusted = NA,
          ci_low_adjusted = NA,
          ci_upper_adjusted = NA,
          pval_adjusted = NA
        ), error = e$message, warning = NA)
      })
      
      # Combine results
      combined_results <- cbind(
        depend_var = depend_var_name,
        independ_var = independ_var_name,
        unadj_result$output,
        adj_result$output,
        formula_unadjusted = formula_unadjusted,
        formula_adjusted = formula_adjusted,
        n_total = n_total,
        n_not_0_depend_var = n_not_0_depend_var,
        n_not_0_independ_var = n_not_0_independ_var,
        error_message_unadjusted = unadj_result$error,
        error_message_adjusted = adj_result$error,
        warning_message_unadjusted = unadj_result$warning,
        warning_message_adjusted = adj_result$warning
      )
      return(combined_results)
    })
  })
  
  # Combine all into single dataframe 
  return(bind_rows(results_list))
}

# -------------------------------------------------------------------------- #  

######   SPLIT INTO TASKS  ---- 

options(future.globals.maxSize = 2000 * 1024^2)  # Increase to 2 GiB 

cores <- as.integer(Sys.getenv("SLURM_NTASKS")) # get n cores 

args <- commandArgs(trailingOnly = TRUE)
index <- as.integer(args[1])

# define ranges of columns to test for each job 

# split n phages to test into 25 batches 
n_jobs <- 25 # Number of array jobs 

n_cols <- length(phage_columns) # Number of phage columns 

cols_per_job_large <- 670  

# Create a list to store the ranges 
ranges <- list()

# First 27 jobs: 1,000 columns each 
for (i in 1:25) {
  start_idx <- (i - 1) * cols_per_job_large + 1
  end_idx <- i * cols_per_job_large
  ranges[[as.character(i)]] <- start_idx:end_idx
}

# Verify the split 
sapply(ranges, length)  # Should show 25 batches of 670 

# -------------------------------------------------------------------------- #  

selected_columns <- ranges[[index]]
if (is.null(selected_columns)) stop("Invalid index: no column range defined.")


# Run the regression 
regression_output <- run_regression_parallel(
  independ_vars = dat_phage[ ,selected_columns],
  dependent_vars = dat_metabs,
  covars = df_analyse %>% select(gm_age, gm_sex, gm_BMI, fid, gm_batch),
  cores = cores
)


# -------------------------------------------------------------------------- #  

# Save results for this job 


path <- file.path("/scratch/users/k2480753/phages/phages_bile_acids", paste0("phages_bile_acid_RESULTS_", index, ".csv"))

write.csv(regression_output,
          path,
          row.names = FALSE)



