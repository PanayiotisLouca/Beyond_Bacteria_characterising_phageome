## Author: 
#  Panayiotis Louca 

## Purpose of script: 
#  

## Date Created: 
#  28 February 2025 

## Notes: 
#  

## Clear environment 

## Set seed 
set.seed(1234)

## Set functions: 


## load up packages: 

### core 
library(tidyverse)

# Maaslin2 
library(Maaslin2)

# -------------------------------------------------------------------------- # 

# ************************* # 
#   IMPORT & PREP DATA   ---- 
# ************************* # 

# -------------------------------------------------------------------------- #  

##   Dataset ---- 

path <- file.path("/scratch/users/k2480753/phages/data/phages_microbes_metabs_full_DATASET.csv")
df <- read.csv(path) %>% as.data.frame(.)

# vectorise columns 
phage_columns = grep("phage_ivig|phage_uvig", names(df), value = TRUE)
covar_columns = c('gm_age','gm_sex','gm_BMI','gm_batch')

# -------------------------------------------------------------------------- #  

# clean variable types 

# phages 
df %>% select(all_of(phage_columns)) %>%
  sapply(., class) %>%
  table

# covars 
df %>% select(all_of(covar_columns)) %>%
  glimpse

df <- df %>%
  mutate(gm_sex = factor(gm_sex),
         gm_batch = factor(gm_batch))

# Calculate Triglyceride-Glucose (TyG) index 
# logarithm (triglycerides (mg/dl) × glucose level (mg/dl)/2) 
# Triglycerides (mmol/L) × 88.57 = mg/dL 
# Glucose (mmol/L) × 18.0182 = mg/dL 
sapply(df %>% select(serum_TG, serum_Gluc), class)

df$TyG_index <- log((df$serum_TG * 88.57 * df$serum_Gluc * 18.0182) / 2)

# -------------------------------------------------------------------------- #  

# *************************************************************************** # 
#   ANALYSE LINKS BETWEEN PHAGES AND METABOLIC PARAMETERS (parallelised)   ---- 
# *************************************************************************** # 

# -------------------------------------------------------------------------- #  

df_analyse <- df %>%
  select(
    iid,
    all_of(covar_columns),
    fid,
    gm_BMI,
    serum_TG,
    TyG_index,
    serum_Gluc,
    # phage viral contigs 
    all_of(phage_columns)
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

df_analyse$serum_TG_transf = rankit(df_analyse$serum_TG)
df_analyse$serum_Gluc_transf = rankit(df_analyse$serum_Gluc)
df_analyse$gm_BMI_transf = rankit(df_analyse$gm_BMI)
df_analyse$TyG_index_tranf = rankit(df_analyse$TyG_index)

# z score age 
df_analyse$gm_age = (df_analyse$gm_age - mean(df_analyse$gm_age))/sd(df_analyse$gm_age)

# CLR normalise phages 

# Create a function to CLR normalise 
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

##   SETUP MAASLIN2 DATASETS ---- 
set.seed(1234)


# -------------------------------------------------------------------------- #  

###   split data ---- 

# gut microbiome data only 
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

# filter to metadata and outcomes 
dat_meta <- df_analyse %>%
  select(iid,
         fid,
         gm_BMI_transf,
         serum_TG_transf,
         serum_Gluc_transf,
         TyG_index_tranf,
         all_of(covar_columns)) 

dat_meta <-  # move iid to rownames 
  dat_meta %>%
  tibble::column_to_rownames('iid')

# ------------------------------------------------------------------------------------------------------------------- # 

# Detect number of cores from SLURM 
cores <- as.integer(Sys.getenv("SLURM_NTASKS"))

###   Run Maaslin2 ---- 

# setup loop 

# Define variables and folders 
variables <- c("gm_BMI_transf", "serum_TG_transf", "serum_Gluc_transf",'TyG_index_tranf')
folders <- c("BMI", "TG", "Glucose", 'TyG_index')
base_dir <- "/scratch/users/k2480753/phages/phage_metab_health"

# Loop through each variable 
for (i in seq_along(variables)) {
  
  var <- variables[i]
  out_dir <- file.path(base_dir, folders[i])
  
  # Run Maaslin2 
  maaslin_results <- Maaslin2(
    input_data = dat_phage,
    input_metadata = dat_meta,
    output = out_dir,
    fixed_effects = c(var, "gm_age", "gm_sex"),
    random_effects = c("gm_batch", "fid"),
    plot_scatter = FALSE,
    plot_heatmap = FALSE,
    normalization = "NONE",
    transform = "NONE",
    min_abundance = 0,
    max_significance = 0.1,
    min_prevalence = 0,
    cores = cores
  )
}    

# -------------------------------------------------------------------------- #  




