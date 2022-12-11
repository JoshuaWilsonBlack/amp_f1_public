# This script generates a distribution plot for the first four PCs given
# intervals of the size defined below and permuted data.
# rewritten from original to perform PCA on rank transformed data.

library(here)
library(tidyverse)
library(glue)

# Set parameters for test.

INTERVAL_SIZE <- 60 # size of desired interval in seconds.
CUT_OFF <- 7 # Number of vowel types required to accept interval. (used by imputation function)
ITERATIONS <- 1000 # Number of times to repeat test.
INCLUDE_AMPLITUDE <- TRUE # Change according to whether amplitude is included as a variable.

set.seed(2000) # Set seed for reproducibility

# Copy required functions from `interval_representation.Rmd`
data_load <- function(filename) {
  df <- read_rds(filename) %>%
    select(
      Speaker,
      Vowel,
      Target.segments.start,
      F1_50,
      F2_50,
      word,
      MeanPitch,
      utterance.articulation.rate,
      participant_gender,
      participant_nz_ethnic,
      participant_age_category,
      following_segment_category,
      intensity_max
    ) %>%
    rename(
      time_unscaled = Target.segments.start,
      articulation_rate = utterance.articulation.rate,
      ethnicity = participant_nz_ethnic,
      gender = participant_gender,
      age_category = participant_age_category
    ) %>%
    ungroup()
}


permute_data <- function(df) {
  df %>%
    group_by(Speaker) %>%
    mutate(
      num_obs = n(),
      time_unscaled = sample(time_unscaled, num_obs, replace=FALSE)
    ) %>%
    arrange(
      time_unscaled, .by_group = TRUE
    ) %>%
    ungroup() %>%
    select(-num_obs)
}


# Interval creation is modified to just fit the desired interval length
# (`INTERVAL_SIZE`).
create_intervals <- function(df) {
  df %>%
    group_by(Speaker) %>%
    mutate(
      interval = as.numeric(
        as.factor(
          cut(
            time_unscaled, 
            breaks = seq(0, max(time_unscaled) + INTERVAL_SIZE, INTERVAL_SIZE))
        )
      )*INTERVAL_SIZE,
    ) %>%
    # Trim terminal intervals with insufficient data. Replaces bad interval values with NA.
    # note: still grouped by Speaker
    mutate(
      speaker_length = max(time_unscaled),
      remaining_from_start = speaker_length - (interval - INTERVAL_SIZE),
      interval = if_else(
        remaining_from_start >= 3/4 * INTERVAL_SIZE, 
        interval, 
        NA_real_
      ),
    ) %>%
    ungroup() %>% 
    # Scale amplitude and articulation rate data at the speaker level.
    group_by(Speaker) %>%
    mutate(
      intensity_max = scale(intensity_max),
      articulation_rate = scale(articulation_rate)
    ) %>%
    # Scale F1 and F2 values by speaker and vowel
    group_by(Speaker, Vowel) %>%
    mutate(
      F1_50 = scale(F1_50),
      F2_50 = scale(F2_50)
    ) %>%
    # Reshape for convenience - now two rows for each vowel token.
    pivot_longer(
      F1_50:F2_50,
      names_to = "formant_type",
      values_to = "formant_value"
    ) %>%
    ungroup() %>%
    # Take summary value for formants
    group_by(Speaker, Vowel, interval, formant_type) %>%
    mutate(
      n_int = n(),
      mean_int = mean(formant_value),
    ) %>%
    ungroup() %>%
    # Take summary value for intensity and articulation rate
    group_by(Speaker, interval) %>%
    mutate(
      art_int = mean(articulation_rate, na.rm=TRUE),
      amp_int = mean(intensity_max, na.rm=TRUE)
    ) %>%
    ungroup()
}

# Converting step in interval_representation.Rmd to a function.
generate_interval_df <-  function(in_df) {
  in_df %>%
    filter(!is.na(interval)) %>% # filter out bad intervals
    group_by(Speaker, Vowel, formant_type, interval) %>%
    summarise(
      n_int = first(n_int),
      mean_int = first(mean_int),
      art_int = first(art_int),
      amp_int = first(amp_int)
    ) %>%
    ungroup() %>%
    # Some vowels may be missing amplitude data, so we add this
    # group_by + mutate to ensure that we get at least one.
    group_by(Speaker, interval) %>%
    mutate(
      amp_int = first(amp_int)
    ) %>%
    ungroup() %>%
    filter( # Some intervals have articulation rate (~70) but no freq. We remove these.
      !is.na(mean_int)
    )
}

# Remove bad speakers (those with no data for a whole vowel type). From 
# `interval_representation.Rmd` but turned into a function.
remove_bad_speakers <- function(in_df) {
  bad_speakers <- in_df %>%
    group_by(Speaker) %>%
    summarise(
      n_vowels = n_distinct(Vowel)
    ) %>%
    filter(
      n_vowels < 10
    ) %>%
    pull(Speaker)
  
  in_df <- in_df %>%
    filter(
      !Speaker %in% bad_speakers
    )
}


impute_formants <- function(df, cut_off) {
  
  # Inputs: 
  #   1) df: dataframe with speaker, interval, vowel, formant_type, mean_int,
  #   art_int, amp_int, n_int columns.
  #   2) cut_off: integer value representing minimum number of vowel types 
  #   required to keep interval.
  
  
  df <- df %>%
    
    # First calculate how many vowel types in each interval.
    group_by(
      Speaker, interval
    ) %>%
    mutate(
      n_vowels = n_distinct(Vowel)
    ) %>%
    ungroup() %>%
    
    # Remove all intervals with less vowel types than the cut off value.
    filter(
      n_vowels >= cut_off
    ) %>%
    
    mutate(
      formant_type = as.factor(formant_type)
    ) %>%
    
    # Group and use summary to create a dataframe with a row for each
    # Speaker, interval, vowel, and formant_type, including those 
    # vowels for which there is no value (note .drop=FALSE ensures
    # that there is a row for all levels of the Vowel factor).
    group_by(Speaker, interval, Vowel, formant_type, .drop=FALSE) %>%
    summarise(
      mean_int = first(mean_int),
      art_int = first(art_int),
      amp_int = first(amp_int),
      n_vowels = first(n_vowels),
      n_int = first(n_int)
    ) %>%
    ungroup() %>%
    
    # Ensure that interval level data is maintained for all levels of Vowel.
    # (It will be set to NA for any of the absent levels in the previous step)
    group_by(Speaker, interval) %>% 
    mutate(
      art_int = max(art_int, na.rm=TRUE),
      amp_int = max(amp_int, na.rm=TRUE),
    ) %>%
    ungroup() %>%
    
    mutate(
      # Calculate imputed value for intervals with 1 or 2 tokens.
      imputed_2 = mean_int * 2/3, # Imputed value if there are 2 real tokens.
      imputed_3 = mean_int * 1/3, # Imputed value if there is 1 real token.
      
      # Any NA values of n_tokens will be those intervals with 0 tokens. We 
      # make this clear now.
      n_int = if_else(is.na(n_int), 0L, n_int),
      
      # Set value to 0 for intervals with no vowel tokens.
      mean_int = if_else(n_int == 0, 0, mean_int), 
      
      # Set value to relevant value for intervals with 1 or 2 tokens of given 
      # vowel.
      mean_int = if_else(n_int == 1, imputed_3, mean_int), 
      mean_int = if_else(n_int == 2, imputed_2, mean_int)
    ) %>%
    
    # Remove intermediate variables.
    select(-c(imputed_2, imputed_3))
}

### Now functions from `corpus_pca.Rmd`
# reshape_for_pca <- function(in_df) {
#   out_df <- in_df %>%
#     ungroup() %>%
#     mutate(
#       vowel_formant = str_c(Vowel, "_", formant_type),
#       speaker_interval = str_c(Speaker, "_", interval)
#     ) %>%
#     select(speaker_interval, vowel_formant, mean_int) %>%
#     pivot_wider(
#       names_from = vowel_formant,
#       values_from = mean_int
#     ) %>%
#     unnest(cols = everything())
# }
reshape_for_pca <- function(in_df, include_amplitude = FALSE) {
  
  if (include_amplitude == FALSE) {
    desired_variables <- c(
      "speaker_interval", "vowel_formant", "mean_int"
    )
  } else {
    desired_variables <- c(
      "speaker_interval", "vowel_formant", "amp_int", "mean_int"
    )
  }
  
  out_df <- in_df %>%
    ungroup() %>%
    mutate(
      vowel_formant = str_c(Vowel, "_", formant_type),
      speaker_interval = str_c(Speaker, "_", interval)
    ) %>%
    select(all_of(desired_variables)) %>%
    pivot_wider(
      names_from = vowel_formant,
      values_from = mean_int
    ) %>%
    unnest(cols = everything())
}

rank_transform <- function(in_df, ties_method = "average") {
  out_df <- in_df %>%
    transmute(
      across(.cols = everything(), .fns = ~ (rank(.x, ties.method=ties_method))) 
    )
}

whole_corpus_pca <- function(in_df) {
  out_pca <- in_df %>%
    rank_transform() %>%
    prcomp(scale=TRUE)
}

# Data load

qb_vowels <- data_load(here('processed_data', 'Quakebox_filtered.rds'))

# Remove foot
qb_vowels <- qb_vowels %>%
  filter(
    !Vowel == "FOOT"
  ) %>%
  mutate(
    Vowel = fct_drop(Vowel)
  )

# Taken from https://stackoverflow.com/questions/13112238/a-matrix-version-of-cor-test
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y, method="spearman", exact=FALSE)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}


### Set up loop to iterate through 1000 times, collecting the amount of variance
### explained by the first 4 PCs and the number of significant, at the 0.05 level
### correlations.

variances_explained <- matrix(ncol=5, nrow=ITERATIONS)
sig_correlations <- rep(0, ITERATIONS)
for (i in 1:ITERATIONS) {
  reshaped <- qb_vowels %>%
    permute_data() %>%
    create_intervals() %>%
    generate_interval_df() %>%
    remove_bad_speakers() %>%
    impute_formants(CUT_OFF) %>%
    reshape_for_pca(include_amplitude = INCLUDE_AMPLITUDE) %>%
    select(-speaker_interval)
    
  # Calculate number of significant pairwise correlations in data.
  cor_test_matrix <- cor.test.p(as.matrix(reshaped))
  sig_correlations[i] <- length(cor_test_matrix[cor_test_matrix < 0.05])
  
  pcaed <- reshaped %>%
    whole_corpus_pca()
  
  variances_explained[i, ] <- (pcaed$sdev^2 / sum(pcaed$sdev^2))[1:5]
}

# 30 times with 60 second intervals.
# 1: Values are not uniquely identified; output will contain list-cols.
# * Use `values_fn = list` to suppress this warning.
# * Use `values_fn = length` to identify where the duplicates arise
# * Use `values_fn = {summary_fun}` to summarise duplicates

variances_explained <- as_tibble(variances_explained) %>%
  rename( # There's a much nicer way to do this
    PC1 = V1,
    PC2 = V2,
    PC3 = V3,
    PC4 = V4,
    PC5 = V5
  )

if (INCLUDE_AMPLITUDE == TRUE) {
  amp_tag <- "_with_amplitude"
} else {
  amp_tag <- ""
}

write_rds(
  variances_explained, 
  here('processed_data', glue('variances_explained_{INTERVAL_SIZE}{amp_tag}.rds'))
)

write_rds(
  sig_correlations, 
  here('processed_data', glue('sig_cors_{INTERVAL_SIZE}{amp_tag}.rds'))
)
