library(tidyverse)
library(here)
library(scales)
library(lme4)

qb_vowels <- read_rds(
  here('processed_data', 'Quakebox_filtered.rds')
)

# copy scaling code

qb_vowels <- qb_vowels %>%
  group_by(Speaker) %>%
  rename(
    time = Target.segments.start,
    art_rate = utterance.articulation.rate,
    pitch = MeanPitch
  ) %>%
  # Within speaker scaling.
  mutate(
    speaker_scaled_time = rescale(time, to = c(0, 1)),
    speaker_length = max(time),
    speaker_scaled_amp_max = scale(intensity_max),
    speaker_scaled_art_rate = scale(art_rate),
    speaker_scaled_pitch = scale(pitch) 
  ) %>%
  ungroup() %>%
  # Across speaker scaling.
  mutate(
    scaled_art_rate = scale(art_rate),
    scaled_pitch = scale(pitch),
    scaled_length = scale(speaker_length)
    # We don't scale amplitude across speakers as we can't control for recording
    # variation.
  )

# Scale formant data
qb_vowels <- qb_vowels %>%
  group_by(Speaker, Vowel) %>%
  mutate(
    speaker_scaled_F1 = scale(F1_50),
    speaker_scaled_F2 = scale(F2_50)
  ) %>%
  # Remove rows with missing F1
  filter(
    !is.na(speaker_scaled_F1)
  ) %>%
  ungroup() %>%
  # We may use across-speaker scaled F1 as a response.
  mutate(
    scaled_F1 = scale(F1_50),
    scaled_F1 = scale(F2_50)
  )

# COPY FILTERING CODE
qb_vowels <- qb_vowels %>%
  filter(
    abs(speaker_scaled_amp_max) <= 2.5
  )

qb_vowels <- qb_vowels %>%
  filter(
    abs(speaker_scaled_art_rate) <= 2.5,
    abs(speaker_scaled_pitch) <= 2.5
  )

qb_vowels <- qb_vowels %>%
  group_by(Speaker, Vowel) %>%
  mutate(
    n_obs = n()
  ) %>%
  ungroup() %>%
  group_by(Speaker) %>%
  mutate(
    min_obs = min(n_obs)
  ) %>%
  filter(
    min_obs >= 5
  )

# COPY DATA REP CODE

qb_vowels <- qb_vowels %>%
  mutate(
    Vowel = as.factor(Vowel),
    Speaker = as.factor(Speaker),
    participant_gender = as.factor(participant_gender)
  )

# COPY TOPIC WRANGELING CODE
qb_vowels <- qb_vowels %>%
  mutate(
    topic = str_to_lower(str_trim(type)), # Remove white space and switch to lower case.,
    topic = str_replace(topic, "earthquakes", "earthquake"),
    topic = str_replace_all(topic, " experience| to the earthquake| of the earthquake| earthquake", "")
  )

qb_vowels %>% pull(topic) %>% unique()

qb_vowels <- qb_vowels %>%
  mutate(
    topic = str_replace_all(topic, "\\{other\\}", "")
  )

qb_vowels %>% pull(topic) %>% unique()

qb_vowels <- qb_vowels %>%
  mutate(
    topic_previous = lag(topic),
    change = if_else(topic_previous != topic, 1, 0)
  ) %>%
  # We now remove untagged sections. Doing this here means we have no risk of
  # missing a gaps between two sections on the same topic.
  filter(
    topic != "" 
  ) %>%
  group_by(Speaker) %>%
  mutate(
    topic_no = cumsum(change)
  ) %>%
  ungroup()

qb_vowels <- qb_vowels %>%
  group_by(Speaker, topic_no) %>%
  mutate(
    topic_length = max(time) - min(time),
    n_in_topic = n()
  ) %>%
  ungroup()

qb_vowels <- qb_vowels %>%
  group_by(Speaker, topic_no) %>%
  mutate(
    topic_time_scaled = rescale(time, to=c(0,1))
  )

# Begin FAKE TOPIC GENERATION CODE
topic_lengths <- qb_vowels %>%
  group_by(Speaker, topic_no) %>%
  summarise(
    topic_length = first(topic_length)
  ) %>%
  filter(
    !is.na(topic_no)
  )

topic_lengths


collect_chunk <- function(speaker, chunk_start, chunk_end, vowel_data) {
  out_df <- vowel_data %>%
    filter(
      Speaker == speaker,
      time > chunk_start,
      time < chunk_end
    ) %>%
    ungroup() %>%
    select(
      c(
        time, art_rate, Vowel, pitch, speaker_scaled_time, 
        speaker_scaled_amp_max, topic_no, intensity_max, 
        speaker_scaled_art_rate, speaker_scaled_pitch
      )
    )
}

generate_chunks <- function(vowel_data) {
  chunks <- vowel_data %>%
    group_by(Speaker, topic_no) %>%
    summarise(
      topic_length = first(topic_length),
      speaker_length = first(speaker_length)
    ) %>%
    filter(
      !is.na(topic_no)
    ) %>%
    mutate(
      chunk_start = map2_int(
        topic_length, 
        speaker_length, 
        ~ sample(seq(0, round(.y - .x)), 1)
      ),
      chunk_end = chunk_start + topic_length
    )
}

### START OF LOOP

mod <- read_rds(here('models', 'glmm_pitch_fit_chunks.rds'))
a <- summary(mod)
coeff_names <- names(a$coefficients[,1])

iterations <- 1000
coefs = matrix(nrow=iterations, ncol=5, dimnames=list(list(), coeff_names))
ses = matrix(nrow=iterations, ncol=5, dimnames=list(list(), coeff_names))
t_vals = matrix(nrow=iterations, ncol=5, dimnames=list(list(), coeff_names))
for (i in 1:iterations) {
  chunks <- generate_chunks(qb_vowels)
  
  chunks <- chunks %>%
    mutate(
      chunk_df = pmap(
        list(Speaker, chunk_start, chunk_end),
        ~ collect_chunk(..1, ..2, ..3, qb_vowels)  
      )
    ) %>%
    rename(
      chunk = topic_no # We may need to compare overlap of real and fake topics.
    ) %>%
    unnest(chunk_df)
  
  chunks <- chunks %>%
    group_by(Speaker, chunk) %>%
    mutate(
      chunk_time_scaled = rescale(time, to=c(0,1)),
      chunk_part = cut(
        chunk_time_scaled, 
        breaks = c(-0.1, 0.33, 0.66, 1.1), 
        labels = c("start", "middle", "end")
      )
    ) %>%
    group_by(Speaker, chunk, chunk_part) %>%
    mutate(
      chunk_part_n = n()
    ) %>%
    ungroup()
  
  chunks_to_filter <- chunks %>%
    filter(
      chunk_part_n < 5
    ) %>%
    select(Speaker, chunk) %>%
    unique()
  
  speaker_chunks_to_filter <- chunks_to_filter %>%
    mutate(
      speaker_chunk = str_c(Speaker, "_", chunk)
    ) %>%
    pull(speaker_chunk)
  
  # Also insisting on having more than 2 chunks
  chunks_filtered <- chunks %>%
    group_by(Speaker) %>%
    mutate(
      speaker_chunk_n = n_distinct(chunk)  
    ) %>%
    ungroup() %>%
    mutate(
      speaker_chunk = str_c(Speaker, "_", chunk)
    ) %>%
    filter(
      !speaker_chunk %in% speaker_chunks_to_filter,
      speaker_chunk_n > 2
    )
  
  glmm_fit <- lmer(
    speaker_scaled_amp_max ~ 
      chunk_part + 
      speaker_scaled_time + 
      speaker_scaled_pitch +
      (0+speaker_scaled_time|Speaker) +
      (0+chunk_part|speaker_chunk), 
    data=chunks_filtered %>%
      mutate(
        speaker_scaled_time = speaker_scaled_time - 0.5
      ))
  
  summ_glmm_fit <- summary(glmm_fit)
  
  coefs[i,] <- summ_glmm_fit$coefficients[,1]
  ses[i,] <- summ_glmm_fit$coefficients[,2]
  t_vals[i,] <- summ_glmm_fit$coefficients[,3]
}

# Think before running me.
write_rds(coefs, here('processed_data', 'glmm_pitch_coeffs_perm.rds'))
write_rds(ses, here('processed_data', 'glmm_pitch_ses_perm.rds'))
write_rds(t_vals, here('processed_data', 'glmm_pitch_tvals_perm.rds'))

# Or me.
coefs <- read_rds(here('processed_data', 'glmm_pitch_coeffs_perm.rds'))
ses <- read_rds(here('processed_data', 'glmm_pitch_ses_perm.rds'))
t_vals <- read_rds(here('processed_data', 'glmm_pitch_tvals_perm.rds'))

# Hard coding coefficients from a model fit on real topics with the 
# same structure as the model fit above.
real_model <- tibble(
  'name' = c('(Intercept)', 'chunk_partmiddle', 'chunk_partend', 'speaker_scaled_time'),
  'value' = c(0.05013, -0.01897, -0.10491, -0.24328)
)

real_model_t <- tibble(
  'name' = c('(Intercept)', 'chunk_partmiddle', 'chunk_partend', 'speaker_scaled_time'),
  'value' = c(2.142, -1.154, -5.283, -4.185)
)

real_model_se <- tibble(
  'name' = c('(Intercept)', 'chunk_partmiddle', 'chunk_partend', 'speaker_scaled_time'),
  'value' = c(0.01589, 0.01910, 0.02192, 0.06591)
)


# Let's plot
coefs %>%
  as_tibble() %>%
  pivot_longer(
    cols = everything(),
    values_to = "value"
  ) %>%
  mutate(
    name = factor(
      name, 
      levels = c(
        '(Intercept)', 'chunk_partmiddle', 'chunk_partend', 'speaker_scaled_time'
      )
    )
  ) %>%
  ggplot(
    aes(
      x = name,
      y = value
    )
  ) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(colour = "red", data = real_model)

violin_three_way <- function(plot_data, comparison_data = list()) {
  out_plot <- plot_data %>%
    as_tibble() %>%
    pivot_longer(
      cols = everything(),
      values_to = "value"
    ) %>%
    mutate(
      name = factor(
        name, 
        levels = c(
          '(Intercept)', 'chunk_partmiddle', 'chunk_partend', 'speaker_scaled_time'
        )
      )
    ) %>%
    ggplot(
      aes(
        x = name,
        y = value
      )
    ) +
    geom_violin(draw_quantiles = c(0.05, 0.5, 0.95))
  
  if (length(comparison_data) > 1) {
    out_plot <- out_plot +
      geom_point(colour = "red", data = comparison_data)
  }
  
  out_plot
}

violin_three_way(coefs, real_model)
violin_three_way(t_vals, comparison_data = real_model_t) +
  labs(
    title = "Distribution of t-values for Random and Topical Segments",
    caption = "Red points indicate values for topical segments",
    y = "t-value",
    x = "Coefficient"
  )
violin_three_way(ses, real_model_se)

coefs %>%
  as_tibble() %>%
  pivot_longer(
    cols = `(Intercept)`:chunk_partend,
    values_to = "coefficient"
  ) %>%
  ggplot(
    aes(
      x = speaker_scaled_time,
      y = coefficient,
      colour = name
    )
  ) +
  geom_point()