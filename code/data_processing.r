###################################################
##### explore_data.r
##### Last Modified by S.Zhang on 2024/10/04
###################################################

library(dplyr)
library(data.table)
library(purrr)

p_data <- "../data"
p_data_mimic_iv <- "../data/mimic_iv"


################################################################################
##### Process mimic IV data
################################################################################

### Load Data
# Vital information
df_vital <- read.csv(file.path(p_data_mimic_iv,"vitalsign.csv"), stringsAsFactors = F)%>%mutate(unique_id = paste0(subject_id, "_", stay_id))
# Diagnosis
df_diag <- read.csv(file.path(p_data_mimic_iv,"diagnosis.csv"), stringsAsFactors = F)%>%mutate(unique_id = paste0(subject_id, "_", stay_id))
# Demographics
df_edstay <- read.csv(file.path(p_data_mimic_iv,"edstays.csv"), stringsAsFactors = F)%>%mutate(unique_id = paste0(subject_id, "_", stay_id))


# Get patients demographic 
table(df_edstay$race)
df_patients_demo <- df_edstay%>%
  mutate(race_r = case_when(grepl("ASIAN", race) ~ "ASIAN_API",
                            grepl("BLACK", race) ~ "BLACK",
                            grepl("HISPANIC", race) ~ "HISPANIC",
                            grepl("WHITE", race) ~ "WHITE",
                            grepl("ISLANDER", race) ~ "ASIAN_API",
                            grepl("ALASKA", race) ~ "ASIAN_API",
                            grepl("UNKNOWN", race) ~ "UNKNOWN_OTHER",
                            grepl("DECLINED", race) ~ "UNKNOWN_OTHER",
                            grepl("UNABLE", race) ~ "UNKNOWN_OTHER",
                            TRUE ~ "UNKNOWN_OTHER"))%>%
  select(unique_id, gender, race_r)%>%unique()


# Get list of patients who diagnosis sepsis
df_patients_sepsis <- df_diag%>%
  mutate(flag_septis = ifelse(grepl("sepsis", icd_title) |grepl("Sepsis", icd_title), 1, 0))%>%
  group_by(unique_id)%>%
  summarise(sepsis = max(flag_septis))
table(df_patients_sepsis$sepsis)
sum(is.na(df_patients_sepsis$sepsis) )

# Get patient that has vital information
df_patients_vital <- df_vital%>%
  select(unique_id)%>%
  unique()


# Get stay length and admission time
df_staylength <- df_edstay%>%
  mutate(time_admission = as.POSIXct(intime, format = "%Y-%m-%d %H:%M:%S"),
         time_discharge = as.POSIXct(outtime, format = "%Y-%m-%d %H:%M:%S"),
         stay_length_hrs = as.numeric(difftime(time_discharge, time_admission, units = "hours")))%>%
  inner_join(df_patients_sepsis)%>%
  inner_join(df_patients_vital)%>%
  inner_join(df_patients_demo)%>%
  select(unique_id, gender, race_r,
         sepsis, time_admission, time_discharge, stay_length_hrs )
table(df_staylength$sepsis)
# 0      1 
# 405796   1464 

### Subsample step1: Excluding Short Stays
# Patients with > 8 hours ED visits
df_stay_limit_hrs <- df_staylength%>%
  filter(stay_length_hrs>=8)
table(df_stay_limit_hrs$sepsis)
sum(is.na(df_stay_limit_hrs$sepsis))
# 0      1 
# 113146    422 

# PS matching on length of stay, race, gender
# Fit a propensity score model (set ratio 1:10 for now)
match_out <- matchit(sepsis ~ stay_length_hrs + race_r + gender, 
                     data = df_stay_limit_hrs, 
                     method = "nearest",  # Use nearest-neighbor matching
                     distance = "logit",
                     ratio = 10,             
                     replace = TRUE)
# View matching summary
summary(match_out)

# Create a matched dataset
matched_data <- match.data(match_out)
table(matched_data$race_r)
glimpse(matched_data)



### Data correction steps, remove outliers (using winsorize way)
df_vital_need <- df_vital%>%
  filter(unique_id%in%matched_data$unique_id)
  
# Summary stats
fn_stats <- function(datain, varneed){
  dataneed <- datain%>%
    mutate(groupind = "allpatients")%>%
    group_by(groupind)%>%
    summarise(varname = varneed,
              mean = mean(!!sym(varneed), na.rm = T),
              min = min(!!sym(varneed), na.rm = T),
              max = max(!!sym(varneed), na.rm = T),
              q001 = quantile(!!sym(varneed), 0.001, na.rm = T),
              q01 = quantile(!!sym(varneed), 0.01, na.rm = T),
              q90 = quantile(!!sym(varneed), 0.90, na.rm = T),
              q95 = quantile(!!sym(varneed), 0.95, na.rm = T),
              q99 = quantile(!!sym(varneed), 0.99, na.rm = T),
              q999 = quantile(!!sym(varneed), 0.999, na.rm = T))
  
  return(dataneed)
}

list_vital <- c("temperature", "heartrate", "resprate", "o2sat", "sbp", "dbp")
data_stats <- data.frame()
for (varneed in list_vital){
  data_stats_h <- fn_stats(df_vital_need, varneed)
  data_stats <- rbind(data_stats, data_stats_h)
}


### Winsorize dataset at 99.9th percentile
fn_winsorize <- function(datain,varneed){
  upperbound <- as.numeric(data_stats[data_stats$varname == varneed, "q999"])
  lowerbound <- as.numeric(data_stats[data_stats$varname == varneed, "q01"])
  dataout <- datain%>%
    mutate(newvar = !!sym(varneed),
           newvar = ifelse(newvar > upperbound, upperbound, newvar),
           newvar = ifelse(newvar < lowerbound, lowerbound, newvar))
  
  names(dataout) <- gsub("newvar", paste0(varneed, "_r"), names(dataout))
  return(dataout)
}

df_vital_winsorize <- df_vital_need
for (varneed in list_vital){
  df_vital_winsorize <- fn_winsorize(df_vital_winsorize, varneed)
}
df_vital_winsorize <- df_vital_winsorize%>%
  mutate(charttime_r = as.POSIXct(charttime, format = "%Y-%m-%d %H:%M:%S"))%>%
  select(unique_id, charttime_r, 
         temperature_r, heartrate_r, resprate_r, o2sat_r, sbp_r, dbp_r)

# Merge all data together
df_master_raw <- df_vital_winsorize%>%
  inner_join(df_staylength)


################################################################################
##### Time encoding
################################################################################
df_master_with_time <- df_master_raw%>%
  arrange(unique_id, charttime_r)%>%
  group_by(unique_id)%>%
  mutate(timepoints = row_number())%>%
  ungroup()%>%
  mutate(time_diff = as.numeric(difftime(charttime_r, time_admission, units = "hours")))%>%
  mutate(time_encoding_1 = sin(time_diff/(8*1)),
         time_encoding_2 = cos(time_diff/(8*2)),
         time_encoding_3 = sin(time_diff/(8*3)),
         time_encoding_4 = cos(time_diff/(8*4)),
         time_encoding_5 = sin(time_diff/(8*5)),
         time_encoding_6 = cos(time_diff/(8*6)),
         time_encoding_7 = sin(time_diff/(8*7)),
         time_encoding_8 = cos(time_diff/(8*8)))


################################################################################
##### Event embedding
################################################################################

# Function to create event embeddings
create_event_embeddings <- function(values, varname, embedding_dim = 8) {
  # Sort the values and replace them with their ranks
  ordered_values <- rank(values, ties.method = "average") / length(values)
  
  # Divide ranks into 10 quantile groups
  quantile_groups <- cut(ordered_values, breaks = seq(0, 1, by = 0.1), labels = FALSE, include.lowest = TRUE)
  
  # Initialize embeddings matrix
  embeddings <- matrix(0, nrow = length(values), ncol = embedding_dim)
  
  # Create embeddings for each event based on its quantile group
  for (i in 1:length(quantile_groups)) {
    group_id <- quantile_groups[i]
    embeddings[i, ] <- rnorm(embedding_dim, mean = group_id, sd = 0.1)  # Example embedding using normal distribution
  }
  
  # Convert matrix to data frame
  df_event_embeddings <- data.frame(embeddings)
  names(df_event_embeddings) <- paste0(varname, 1:8)
  
  return(df_event_embeddings)
}



# Temperature
temperature_embeddings <- create_event_embeddings(df_master_with_time$temperature_r, "temperature_")
# heartrate_
heartrate_embeddings <- create_event_embeddings(df_master_with_time$heartrate_r, "heartrate_")
# resprate_
resprate_embeddings <- create_event_embeddings(df_master_with_time$resprate_r, "resprate_")
# o2sat_
o2sat_embeddings <- create_event_embeddings(df_master_with_time$o2sat_r, "o2sat_")
# sbp_
sbp_embeddings <- create_event_embeddings(df_master_with_time$sbp_r, "sbp_")
# dbp_
dbp_embeddings <- create_event_embeddings(df_master_with_time$dbp_r, "dbp_")


df_analytical_ready <- cbind(df_master_with_time,
                             temperature_embeddings, heartrate_embeddings, resprate_embeddings,
                             o2sat_embeddings, sbp_embeddings, dbp_embeddings)
                             

write.csv(df_analytical_ready, file.path(p_data, "analytical_ready_1115.csv"), row.names = F)



################################################################################
##### Create Table 1
################################################################################
df_demo_with_flag <- df_analytical_ready%>%
  select(unique_id, sepsis, 
         gender, race_r)%>%
  unique()

df_demo_table_gender <- df_demo_with_flag%>%
  mutate(description = "gender",
         group = gender)%>%
  group_by(description, sepsis, group)%>%
  summarise(n_bene = n())

df_demo_table_race <- df_demo_with_flag%>%
  mutate(description = "race",
         group = race_r)%>%
  group_by(description, sepsis, group)%>%
  summarise(n_bene = n())

n_sepsis_yes <- sum(df_demo_with_flag$sepsis)
n_sepsis_no <- nrow(df_demo_with_flag) - n_sepsis_yes
df_demo_table1_h0 <- rbind(df_demo_table_gender, df_demo_table_race)%>%filter(sepsis == 0)%>%
  mutate(pct_bene_no = round(n_bene/n_sepsis_no,3),
         n_bene_no = n_bene)%>%ungroup()%>%
  select(description, group,  n_bene_no, pct_bene_no)
df_demo_table1_h1 <- rbind(df_demo_table_gender, df_demo_table_race)%>%filter(sepsis == 1)%>%
  mutate(pct_bene_yes = round(n_bene/n_sepsis_yes,3),
         n_bene_yes = n_bene)%>%ungroup()%>%
  select(description, group,  n_bene_yes, pct_bene_yes)

df_demo_table <- df_demo_table1_h0%>%
  left_join(df_demo_table1_h1)


##### Create Time segment #####
df_sepsis_qsofa <- df_analytical_ready%>%
  group_by(unique_id)%>%
  filter(max(flag_onset_qsofa) >= 1)%>%
  ungroup()

df_sepsis_qsofa_admission <- df_sepsis_qsofa%>%
  group_by(unique_id)%>%
  arrange(charttime_r)%>%
  slice(1)%>%ungroup()%>%
  mutate(charttime_admission = charttime_r)%>%
  select(unique_id, charttime_admission)


df_sepsis_qsofa_onset <- df_sepsis_qsofa%>%
  filter(flag_onset_qsofa == 1)%>%
  mutate(charttime_onset = charttime_r)%>%
  select(unique_id, charttime_onset)


df_sepsis_qsofa_segment <- df_sepsis_qsofa_admission%>%
  left_join(df_sepsis_qsofa_onset)%>%
  mutate(segment_length = as.numeric(charttime_onset - charttime_admission))


df_sepsis_check_timepoint <- df_sepsis_qsofa%>%
  group_by(unique_id)%>%
  summarise(n_timepoint = n())





