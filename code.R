
# INIT ----
#+++++++++++++++++++++++++++++

DATADIR='D:\\yehuda\\data\\'
RESDIR='D:\\yehuda\\results\\'

STUDY_MAX_PCR_TO_HOSP = 21
STUDY_MIN_PCR_TO_HOSP = -2
STUDY_MAX_PCR_TO_DEATH = 60
STUDY_MIN_TO_DEATH = 0
STUDY_MIN_AGE = 12
MIN_N = 5

AGE_GROUPS_BREAKS_NARROW = c( 12, 35, 45, 55, 65, 75, 85, Inf)
AGE_GROUPS_LABELS_NARROW = c('12-34', '35-44', '45-54', '55-64', '65-74', '75-84', '85+')

AGE_GROUPS_BREAKS = c(12, 55, 65, Inf)
AGE_GROUPS_LABELS = c('12-54', '55-64', '65+')

REASONOFDEATHNAME_HEB <- c('לא ידוע', 'נפטר שלילי לקורונה', 'נפטר מסיבה אחרת', 'נפטר כחשוד', 'נפטר מקורונה')
REASONOFDEATHNAME <- c('unknown', 'negative', 'other', 'suspect', 'positive')

MEDICAL_SITUATION_HEB <- c('לא ידוע', 'החלים', 'קל', 'בינוני', 'קשה', 'קריטי', 'נפטר')
MEDICAL_SITUATION <- c('unknown', 'recovered', 'mild', 'moderate', 'severe', 'critical', 'dead')

DIVISIONTYPENAME_HEB = c('רגילה', 'ייעודית')
DIVISIONTYPENAME = c('stndard', 'dedicated')

STUDY_START = '2020-12-20'
STUDY_END = '2021-10-07'
CUTTOF_3T4WAVE = '2021-05-01'

# Put output in file
ofile <- file(paste0(RESDIR, 'output.txt'), open = "wt")
sink(ofile, type = "output")
sink(ofile, type = "message")

library(plyr)
library(tidyr)
library(MatchIt)
library(ggplot2)
library(survival)
library(survminer)
library(MASS)
library(dplyr)
library(cobalt)
library(stringr)

set.seed(12345)

# READ TABLES ----
#+++++++++++++++++++++++++++++

source.master <- read.table(paste0(DATADIR, 'galon_master.csv'), header = TRUE, sep='|', encoding = 'UTF-8', na.strings = 'NULL')
source.hospitalizations <- read.table(paste0(DATADIR, 'galon_hospitalizations.csv'), header = TRUE, sep='|', encoding = 'UTF-8', na.strings = 'NULL')
source.monitoring <- read.table(paste0(DATADIR, 'galon_monitoring.csv'), header = TRUE, sep='|', encoding = 'UTF-8', na.strings = 'NULL')
source.quarantine <- read.table(paste0(DATADIR, 'galon_quarantine.csv'), header = TRUE, sep='|', encoding = 'UTF-8', na.strings = 'NULL')

unique.master <- source.master %>% distinct()
unique.hospitalizations <- source.hospitalizations %>% distinct()
unique.monitoring <- source.monitoring %>% distinct()
unique.quarantine <- source.quarantine %>% distinct()

data.master <- unique.master
data.hospitalizations <- unique.hospitalizations
data.monitoring <- unique.monitoring
data.quarantine <- unique.quarantine

# CONVERT VALUES ----
#+++++++++++++++++++++++++++++

## master ----
#.............................

data.master['first_positive_result_test_date'][data.master['first_positive_result_test_date']== '2000-01-01'] <- NA
data.master$gender <- mapvalues(
  data.master$gendercodename, from = c('1', '2', '999'), to = c('male', 'female', NA)
)
data.master$ind_death <- mapvalues(
  data.master$ind_death, from = c(NA, '0', '1'), to = c('0', '0', '1')
)
data.master$positive_pcr_within_5 <- mapvalues(
  data.master$positive_pcr_within_5, from = c(NA, 'Y'), to = c('0', '1')
)
data.master$positive_pcr_before_vaccine <- mapvalues(
  data.master$positive_pcr_before_vaccine, from = c(NA, 'Y'), to = c('0', '1')
)
data.master$positive_pcr_no_vaccine <- mapvalues(
  data.master$positive_pcr_no_vaccine, from = c(NA, 'Y'), to = c('0', '1')
)
data.master$new_reasonofdeathname <- mapvalues(
  data.master$new_reasonofdeathname, from = REASONOFDEATHNAME_HEB, to = REASONOFDEATHNAME
)

## monitoring ----
#.............................

data.monitoring$first_update_medical_situation <- mapvalues(
  data.monitoring$first_update_medical_situation, from = MEDICAL_SITUATION_HEB, to = MEDICAL_SITUATION
)

data.monitoring$last_update_medical_situation <- mapvalues(
  data.monitoring$last_update_medical_situation, from = MEDICAL_SITUATION_HEB, to = MEDICAL_SITUATION
)

## hospitalizations ----
#.............................

data.hospitalizations$new_divisiontypename <- mapvalues(
  data.hospitalizations$new_divisiontypename, from = DIVISIONTYPENAME_HEB, to = DIVISIONTYPENAME
)

# SET TYPES ----
#+++++++++++++++++++++++++++++

## master ----
#.............................

data.master$first_positive_result_test_date <- as.Date(data.master$first_positive_result_test_date)
data.master$positive_pcr_within_5 <- as.factor(data.master$positive_pcr_within_5)
data.master$positive_pcr_before_vaccine <- as.factor(data.master$positive_pcr_before_vaccine)
data.master$positive_pcr_no_vaccine <- as.factor(data.master$positive_pcr_no_vaccine)
data.master$new_birthdate <- as.numeric(data.master$new_birthdate)
data.master$gender <- as.factor(data.master$gender)
data.master$ind_death <- as.factor(data.master$ind_death)
data.master$new_reasonofdeathname <- as.factor(data.master$new_reasonofdeathname)
data.master$first_vaccincationdatetime <- as.Date(data.master$first_vaccincationdatetime)
data.master$second_vaccinationdatetime <- as.Date(data.master$second_vaccinationdatetime)

## hospitalizations ----
#.............................

data.hospitalizations$new_deathdate <- as.Date(data.hospitalizations$new_deathdate)
data.hospitalizations$new_arrivaldate <- as.Date(data.hospitalizations$new_arrivaldate)
data.hospitalizations$new_releasedate <- as.Date(data.hospitalizations$new_releasedate)
data.hospitalizations$duration_hospitalization <- as.numeric(data.hospitalizations$duration_hospitalization)
data.hospitalizations$new_divisiontypename <- as.factor(data.hospitalizations$new_divisiontypename)

## monitoring ----
#.............................

data.monitoring$first_datemonitoring <- as.Date(data.monitoring$first_datemonitoring)
data.monitoring$first_update_medical_situation <- as.factor(data.monitoring$first_update_medical_situation)
data.monitoring$last_datemonitoring <- as.Date(data.monitoring$last_datemonitoring)
data.monitoring$last_update_medical_situation <- as.factor(data.monitoring$last_update_medical_situation)

## quarantine ----
#.............................

data.quarantine$new_fever <- as.factor(data.quarantine$new_fever)
data.quarantine$isolationstartday <- as.Date(as.character(data.quarantine$isolationstartdaykey), format = '%Y%m%d')

# ADD ERRORS ----
#+++++++++++++++++++++++++++++

## Init err columns ----
#.............................

data.master$err <- ''
data.hospitalizations$err <- ''
data.monitoring$err <- ''
data.quarantine$err <- ''

## master ----
#.............................

data.master <- data.master %>%
# Set err if duplicate patient id in master table
  mutate(
    err = case_when(
      (patient_id %in% patient_id[duplicated(patient_id)]) ~ paste0(err, ' duplicated'),
      TRUE ~ err
    )
  ) %>%
  # Add error if no patient_id
  mutate(
    err = case_when(
      is.na(patient_id) ~ paste0(err, ' no_patient_id'),
      TRUE ~ err
    )
  ) %>%
  # Add error if no first_positive_result_test_date
  mutate(
    err = case_when(
      is.na(first_positive_result_test_date) ~ paste0(err, ' no_first_positive_result_test_date'),
      TRUE ~ err
    )
  ) %>%
  # Add error if no gendercodename
  mutate(
    err = case_when(
      is.na(gendercodename) ~ paste0(err, ' no_gendercodename'),
      TRUE ~ err
    )
  ) %>%
  # Add error if no new_birthdate
  mutate(
    err = case_when(
      is.na(new_birthdate) ~ paste0(err, ' no_new_birthdate'),
      TRUE ~ err
    )
  )

## hospitalizations ----
#.............................

data.hospitalizations <- data.hospitalizations %>%
  # Add error if different death dates for same patient
  mutate(
    err = case_when(
      (patient_id %in% filter(distinct(drop_na(.[,c('patient_id', 'new_deathdate')])), duplicated(patient_id))[,'patient_id']) ~
        paste0(err, ' diff_deathdate'),
      TRUE ~ err
    )
  ) %>%
  # Error if releasedate < arrivaldate
  mutate(
    err = case_when(
      !is.na(new_releasedate) & !is.na(new_arrivaldate) & new_releasedate < new_arrivaldate ~
        paste0(err, ' releasedate_arrivaldate'),
      TRUE ~ err
    )
  ) %>%
  # Error if arrivaldate - deathdate < STUDY_MIN_TO_DEATH days
  mutate(
    err = case_when(
      !is.na(new_deathdate) & !is.na(new_arrivaldate) & difftime(new_deathdate, new_arrivaldate, units = 'days') < STUDY_MIN_TO_DEATH ~
        paste0(err, ' deathdate_arrivaldate'),
      TRUE ~ err
    )
  ) %>%
  # Get first_positive_result_test_date from master table. !patient id is unique without errors!
  merge(
    data.master[data.master$err == '',c('patient_id', 'first_positive_result_test_date')],
    by = 'patient_id', all.x = TRUE, all.y = FALSE
  ) %>%
  # TODO split so it'll be more readable.
  # Error if first_positive_result_test_date - deathdate < STUDY_MIN_TO_DEATH days
  mutate(
    err = case_when(
      !is.na(new_deathdate) & !is.na(first_positive_result_test_date) & difftime(new_deathdate, first_positive_result_test_date, units = 'days') < STUDY_MIN_TO_DEATH ~
        paste0(err, ' deathdate_first_positive_result_test_date'),
      TRUE ~ err
    )
  ) %>%
  # Remove first_positive_result_test_date column
  select(-first_positive_result_test_date)

# MAKE DATA TABLE ----
#+++++++++++++++++++++++++++++

## master ----
#.............................

data.comb <- data.master %>% filter(err == '') %>% select(-err)

## death dates ----
#.............................

# Generate table with death dates for all patient who die from corona (in hospital)
tmp.hospitalizations.deathdate <- data.hospitalizations %>%
  # TODO change to only if include 'deathdate' instead of twice.
  # Remove errors
  filter(!grepl('.* diff_deathdate.*|.* deathdate_first_positive_result_test_date.*', err)) %>%
  # Select the relevant columns
  dplyr::select(c('patient_id', 'new_deathdate')) %>%
  # Remove incomplete cases
  drop_na() %>%
  # Remove duplicates (if patient hospitalized twice)
  distinct() %>%
  # Leave only patient exist in master table with ind_death=1
  merge(
    data.comb[data.comb$ind_death == '1',c('patient_id'), drop = FALSE],
    by = 'patient_id',
    all = FALSE
  )

# Add deathdate to data.comb
data.comb <- data.comb %>%
  merge(tmp.hospitalizations.deathdate, by = 'patient_id', all.x = TRUE, all.y = FALSE)

## hospitalizations ----
#.............................

# Generate table with relevant hospitalization for corona
tmp.hospitalizations <- data.hospitalizations %>%
  # Leave only hospitalization in dedicated division with new_arrivaldate
  filter(new_divisiontypename == 'dedicated' & !is.na(new_arrivaldate)) %>%
  # Select the relevant columns
  dplyr::select(c('patient_id', 'new_arrivaldate', 'new_releasedate', 'err')) %>%
  # Remove duplicates
  distinct() %>%
  # Get first_positive_result_test_date from master table.
  merge(data.comb[, c('patient_id', 'first_positive_result_test_date')], by = 'patient_id', all = FALSE) %>%
  # Leave only arrival date during cutoff from positive PCR and release date after positive PCR
  filter(
    between(as.numeric(difftime(new_arrivaldate, first_positive_result_test_date, units = 'days')), STUDY_MIN_PCR_TO_HOSP, STUDY_MAX_PCR_TO_HOSP) &
      (is.na(new_releasedate) | first_positive_result_test_date <= new_releasedate)
  )

## hospitalization indication ----
#.............................

# TODO maybe take max so if patient have at least one correct hospitalization he'll indicated as hospitalized.
# Generate table with indication for hospitalization for corona
tmp.hospitalizations.ind <- tmp.hospitalizations %>%
  # Create hospitalization indication: -1 on error and 1 if hospitalized
  mutate(
    ind_hospitalized = case_when(
      err != '' ~ '-1',
      TRUE ~ '1'
    )
  ) %>%
  # Select the relevant columns
  dplyr::select(c('patient_id', 'ind_hospitalized')) %>%
  # group by patient_id and take -1 if exists
  group_by(patient_id) %>%
  dplyr::summarise(ind_hospitalized = min(ind_hospitalized))

# Add ind_hospitalized to data.comb
data.comb <- data.comb %>%
  merge(tmp.hospitalizations.ind, by = 'patient_id', all.x = TRUE, all.y = FALSE) %>%
  # Update all NA (not hospitalized) to 0
  mutate(ind_hospitalized = as.factor(case_when(is.na(ind_hospitalized) ~ '0', TRUE ~ ind_hospitalized)))

## hospitalization duration ----
#.............................

# Generate table with hospitalization duration
tmp.hospitalizations.duration <- tmp.hospitalizations %>%
  # Remove errors
  filter(err == '') %>%
  # Select the relevant columns
  dplyr::select(c('patient_id', 'new_arrivaldate', 'new_releasedate', 'first_positive_result_test_date')) %>%
  # Remove incomplete cases
  drop_na() %>%
  # Merge overlapping
  arrange(patient_id, new_arrivaldate, new_releasedate) %>%
  group_by(patient_id) %>%
  dplyr::mutate(idx= c(0, cumsum(as.numeric(lead(new_arrivaldate)) >
                          cummax(as.numeric(new_releasedate)))[-n()])) %>%
  group_by(patient_id, idx) %>%
  dplyr::summarise(new_arrivaldate = min(new_arrivaldate), new_releasedate = max(new_releasedate), .groups = 'drop') %>%
  dplyr::select(-idx) %>%
  group_by(patient_id) %>%
  dplyr::summarise(duration_hospitalization = sum(difftime(new_releasedate, new_arrivaldate, units = 'days')))

# Tests
# /WKAesFL1iaYv1NkYXb6+w== NA -overlap [19]
# //k2nt0e3IvsKNJNv6P2QA== 1 
# /2DoR3ckPWhD0adlFPZ/8Q== 5

# Add duration_hospitalization to data.comb
data.comb <- data.comb %>%
  merge(tmp.hospitalizations.duration, by = 'patient_id', all.x = TRUE, all.y = FALSE)

## hospitalization first ----
#.............................

# Generate table with the first hospitalization in the cuttof range from positive PCR
tmp.hospitalizations.first <- tmp.hospitalizations %>%
  # Remove errors
  filter(err == '') %>%
  # Select the relevant columns
  dplyr::select(c('patient_id', 'new_arrivaldate')) %>%
  # Remove incomplete cases
  drop_na() %>%
  # Get the first hospitalization for each patient
  group_by(patient_id) %>%
  dplyr::summarise(new_arrivaldate_first = min(new_arrivaldate))

# Add new_arrivaldate_first to data.comb
data.comb <- data.comb %>%
  merge(tmp.hospitalizations.first, by = 'patient_id', all.x = TRUE, all.y = FALSE)

## Add columns ----
#.............................

data.comb <- data.comb %>%
  mutate(
    age = 2021 - new_birthdate,
    days_first_vaccincation_to_first_positive_result = as.numeric(difftime(
      first_positive_result_test_date,
      first_vaccincationdatetime,
      units = 'days'
    )),
    days_first_positive_result_to_death = as.numeric(difftime(
      new_deathdate,
      first_positive_result_test_date,
      units = 'days'
    )),
    days_first_positive_result_to_hospitalization = as.numeric(difftime(
      new_arrivaldate_first,
      first_positive_result_test_date,
      units = 'days'
    )),
    first_positive_result_test_week = strftime(first_positive_result_test_date, format = "%Y%V"),
    first_positive_result_test_month = strftime(first_positive_result_test_date, format = "%Y%m"),
    case = as.factor(case_when(
      positive_pcr_within_5 == '1' ~ '1',
      positive_pcr_no_vaccine =='1' ~ '0'
    )),
    age_group = cut(
      age,
      right = FALSE,
      breaks = AGE_GROUPS_BREAKS,
      labels = AGE_GROUPS_LABELS
    ),
    age_group_narrow = cut(
      age,
      right = FALSE,
      breaks = AGE_GROUPS_BREAKS_NARROW,
      labels = AGE_GROUPS_LABELS_NARROW
    )
  )

## Make matched study table ----
#.............................

# case = all patient who get positive pcr test in 5 days since the first vaccination during study period
# control = all patient who didn't get the vaccine and get positive pcr test during study period

# Get study table before match
study.pre.match <- data.comb %>%
  # Leave only complete cases or control
  drop_na(case, age_group, gender, first_positive_result_test_date) %>%
  # Leave only cases with positive pcr during study period
  filter(between(first_positive_result_test_date, as.Date(STUDY_START), as.Date(STUDY_END))) %>%
  # Leave only patients aged STUDY_MIN_AGE and above
  filter(STUDY_MIN_AGE <= age)

# Do match
study.matched <- matchit(
  case ~ age_group_narrow + gender + age + first_positive_result_test_date,
  data = study.pre.match, method = 'nearest', ratio = 1, exact = c("age_group_narrow", "gender")
)

# Make study table
study <- match.data(study.matched) %>%
  arrange(desc(case), age, gender)

# ANALYSIS ----
#+++++++++++++++++++++++++++++

## Match Covariate ----
#.............................

plt.covariate <- love.plot(
  bal.tab(study.matched),
  stats = 'mean.diffs',
  stars = "raw",
  var.names = data.frame(
    old=c(
      paste("age_group_narrow_",levels(study.pre.match$age_group_narrow),sep=""),
      'gender_male',
      'gender_female',
      'age',
      'first_positive_result_test_date'
    ),
    new=c(
      levels(study.pre.match$age_group_narrow),
      'Gender',
      'Gender',
      'Age',
      'PCR-positive\n result date'
    )
  ),
  thresholds = 0.1
)

## Case control months distribution ----
#.............................

# Here we test the distribution of positive pcr over time of our case and control

plt.first_positive_result_test_month <- data.comb %>%
  drop_na(case) %>%
  {ggplot(., aes(x = first_positive_result_test_month, fill = case)) +
    geom_bar(position = 'identity', alpha = 0.4,  stat="count") +
    theme_light() +
    scale_fill_discrete(name="",labels=c("Unvaccinated","Recently Injected")) +
    xlab("Month of first positive Covid19 test")+ ylab("Counts") +
    facet_grid(case ~ ., scales = 'free') +
    theme(strip.text.y =element_blank(), axis.text.x = element_text(angle = 45)) +
    scale_x_discrete(
      breaks = unique(.$first_positive_result_test_month),
      labels = gsub("(\\d{4})(\\d{2})", "\\2-\\1", unique(.$first_positive_result_test_month))
    )}

## General days to hospitalization distribution ----
#.............................

# Here we test the distribution of the days from positive pcr to first hospitalization of our case and control
plt.days_first_positive_result_to_hospitalization <- data.hospitalizations %>%
  # Leave only rows without errors and in corona dedicated division
  filter(err == '' & new_divisiontypename == 'dedicated') %>%
  # Add data from master table
  merge(
    data.comb[, c('patient_id', 'first_positive_result_test_date', 'positive_pcr_within_5', 'positive_pcr_no_vaccine')],
    by = 'patient_id', all = FALSE
  ) %>%
  # Add coluns
  mutate(
    # Add case control indication column
    case = as.factor(case_when(
      positive_pcr_within_5 == '1' ~ '1',
      positive_pcr_no_vaccine =='1' ~ '0')),
    # Add days from first positive PCR to hospitalization
    days_first_positive_result_to_hospitalization = as.numeric(difftime(
      new_arrivaldate,
      first_positive_result_test_date,
      units = 'days'))
  ) %>%
  # Remove incomplete cases
  drop_na(days_first_positive_result_to_hospitalization) %>%
  # filter only in 20 days range from positive PCR to hospitalization for the graph to be readable
  filter(between(days_first_positive_result_to_hospitalization, -20, 20)) %>%
  # Plot
  {ggplot(., aes(x = days_first_positive_result_to_hospitalization, fill = case)) +
      geom_bar(position = 'identity', alpha = 0.4,  stat="count") +
      theme_light() +
      scale_fill_discrete(name="",labels=c("Unvaccinated","Recently Injected")) +
      xlab("Number of days from positive Covid19 test to hospitalization")+
      ylab("Counts") +
      theme(strip.text.y =element_blank())}

## Study days to hospitalization distribution ----
#.............................

# Here we test the distribution of the days from positive pcr to first hospitalization of our case and control in our study data

# Generate table with hospitalization
plt.days_first_positive_result_to_hospitalization_study <- study %>%
  # Remove incomplete cases
  drop_na(days_first_positive_result_to_hospitalization) %>%
  # filter only in 20 days range from positive PCR to hospitalization for the graph to be readable
  filter(between(days_first_positive_result_to_hospitalization, -20, 20)) %>%
  # Plot
  {ggplot(., aes(x = days_first_positive_result_to_hospitalization, fill = case)) +
      geom_bar(position = 'identity', alpha = 0.4,  stat="count") +
      theme_light() +
      scale_fill_discrete(name="",labels=c("Unvaccinated","Recently Injected")) +
      xlab("Number of days from positive Covid19 test to hospitalization in study") +
      ylab("Counts") +
      facet_grid(case ~ .) +
      theme(strip.text.y =element_blank())}

## Death ----
#.............................

# Here we test she correlation between vaccination and the probability to die from corona

anl.death <- bind_rows(
  
  # Narrow age group
  study %>%
    mutate(description='Narrow age group', age_group=age_group_narrow) %>%
    group_by(description, age_group, gender) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
    .groups = 'drop'),
  
  # Wider age group
  study %>%
    mutate(description='Wider age group') %>%
    group_by(description, age_group, gender) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Wider age group no gender
  study %>%
    mutate(description='Wider age group') %>%
    group_by(description, age_group) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Aged 55 and above by gender
  study %>%
    # Leave only patient Aged 55 and above
    filter(age >= 55) %>%
    mutate(description='Patient aged 55 and above') %>%
    group_by(description, gender) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Aged 55 and above
  study %>%
    # Leave only patient Aged 55 and above
    filter(age >= 55) %>%
    mutate(description='Patient aged 55 and above') %>%
    group_by(description) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop'),
)

anl.death$description <- as.factor(anl.death$description)

## Death when < 2 days vaccination to positive pcr ----
#.............................

# Here we test across gender and age groups the correlation between vaccination and the probability to die from corona

# Make study table
study.2d <- study %>%
  # Leave only cases with less then 2 days from vaccination to positive pcr
  filter(case == '1' & between(days_first_vaccincation_to_first_positive_result, 0, 1)) %>%
  # Add the matched control
  bind_rows(study[study$case == '0' & study$subclass %in% .$subclass,]) %>%
  # Sort table
  arrange(desc(case), age, gender)

# Make result table
anl.death2d <- bind_rows(
  
  # Narrow age group
  study.2d %>%
    mutate(description='Narrow age group', age_group=age_group_narrow) %>%
    group_by(description, age_group, gender) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Wider age group
  study.2d %>%
    mutate(description='Wider age group') %>%
    group_by(description, age_group, gender) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Wider age group no gender
  study.2d %>%
    mutate(description='Wider age group') %>%
    group_by(description, age_group) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Aged 55 and above by gender
  study.2d %>%
    # Leave only patient Aged 55 and above
    filter(age >= 55) %>%
    mutate(description='Patient aged 55 and above') %>%
    group_by(description, gender) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Aged 55 and above
  study.2d %>%
    # Leave only patient Aged 55 and above
    filter(age >= 55) %>%
    mutate(description='Patient aged 55 and above') %>%
    group_by(description) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop')
)

anl.death2d$description <- as.factor(anl.death2d$description)

## Death when 2-5 days vaccination to positive pcr ----
#.............................

# Here we test across gender and age groups the correlation between vaccination and the probability to die from corona

# Make study table
study.2t5d <- study %>%
  # Leave only cases with 2-5 days from vaccination to positive pcr
  filter(case == '1' & between(days_first_vaccincation_to_first_positive_result, 2, 5)) %>%
  # Add the matched control
  bind_rows(study[study$case == '0' & study$subclass %in% .$subclass,]) %>%
  # Sort table
  arrange(desc(case), age, gender)

# Make result table
anl.death2t5d <- bind_rows(
  
  # Narrow age group
  study.2t5d %>%
    mutate(description='Narrow age group', age_group=age_group_narrow) %>%
    group_by(description, age_group, gender) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Wider age group
  study.2t5d %>%
    mutate(description='Wider age group') %>%
    group_by(description, age_group, gender) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Wider age group no gender
  study.2t5d %>%
    mutate(description='Wider age group') %>%
    group_by(description, age_group) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Aged 55 and above by gender
  study.2t5d %>%
    # Leave only patient Aged 55 and above
    filter(age >= 55) %>%
    mutate(description='Patient aged 55 and above') %>%
    group_by(description, gender) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Aged 55 and above
  study.2t5d %>%
    # Leave only patient Aged 55 and above
    filter(age >= 55) %>%
    mutate(description='Patient aged 55 and above') %>%
    group_by(description) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.death.n=sum(case=='1' & ind_death=='1'),
      control.death.n=sum(case=='0' & ind_death=='1'),
      case.death.p=sum(case=='1' & ind_death=='1')/sum(case=='1'),
      control.death.p=sum(case=='0' & ind_death=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, ind_death)$p.value,
      fisher.test.or=fisher.test(case, ind_death)$estimate,
      fisher.test.conf.int=paste(fisher.test(case, ind_death)$conf.int, collapse = '-'),
      .groups = 'drop')
)

anl.death2t5d$description <- as.factor(anl.death2t5d$description)

## Survival ----
#.............................

# Make data for survival analysis to cutoff time
# Make sure nrow(study %>% filter(ind_death == '1' & is.na(new_deathdate))) == 0

study.survival <- study %>%
  mutate(
    ind_death = as.numeric(case_when(
      days_first_positive_result_to_death > STUDY_MAX_PCR_TO_DEATH ~ '0',
      TRUE ~ as.character(ind_death)
    )),
    days_first_positive_result_to_death = case_when(
      is.na(days_first_positive_result_to_death) ~ STUDY_MAX_PCR_TO_DEATH,
      days_first_positive_result_to_death > STUDY_MAX_PCR_TO_DEATH ~ STUDY_MAX_PCR_TO_DEATH,
      TRUE ~ days_first_positive_result_to_death
    ),
    days_first_positive_result_to_death = days_first_positive_result_to_death - STUDY_MIN_TO_DEATH
  )

# General

plt.survival <- study.survival %>%
  {ggsurvplot(
    fit = survfit(
      Surv(days_first_positive_result_to_death, ind_death) ~ case,
      data= .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    risk.table = TRUE,
    risk.table.title = "No. at Risk",
    table.heigth=.45,
    tables.theme = theme_cleantable(),
    pval=TRUE,
    pval.coord = c(1,0.65),
    ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since positive PCR",
    ylab = "Survival probability",
    break.time.by = 10,
    xlim = c(0, STUDY_MAX_PCR_TO_DEATH - STUDY_MIN_TO_DEATH + 5)
  )}

# Aged 55 and above

plt.survival55p <- study.survival %>%
  # Leave only patient Aged 55 and above
  filter(age >= 55) %>%
  {ggsurvplot(
    fit=survfit(
      Surv(days_first_positive_result_to_death, ind_death) ~ case,
      data=  .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    risk.table = TRUE,
    risk.table.title = "No. at Risk",
    table.heigth=.45,
    tables.theme = theme_cleantable(),
    pval=TRUE,
    pval.coord = c(1,0.65),
    ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since positive PCR",
    ylab = "Survival probability",
    break.time.by = 10,
    xlim = c(0, STUDY_MAX_PCR_TO_DEATH - STUDY_MIN_TO_DEATH + 5)
  )}

# Aged 65 and above

plt.survival65p <- study.survival %>%
  # Leave only patient Aged 65 and above
  filter(age >= 65) %>%
  {ggsurvplot(
    fit=survfit(
      Surv(days_first_positive_result_to_death, ind_death) ~ case,
      data=  .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    risk.table = TRUE,
    risk.table.title = "No. at Risk",
    table.heigth=.45,
    tables.theme = theme_cleantable(),
    pval=TRUE,
    pval.coord = c(1,0.65),
    ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since positive PCR",
    ylab = "Survival probability",
    break.time.by = 10,
    xlim = c(0, STUDY_MAX_PCR_TO_DEATH - STUDY_MIN_TO_DEATH + 5)
  )}

# Aged 55 - 65

plt.survival55_65 <- study.survival %>%
  # Leave only patient Aged 55-65 and above
  filter(between(age, 55, 64)) %>%
  {ggsurvplot(
    fit=survfit(
      Surv(days_first_positive_result_to_death, ind_death) ~ case,
      data=  .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    risk.table = TRUE,
    risk.table.title = "No. at Risk",
    table.heigth=.45,
    tables.theme = theme_cleantable(),
    pval=TRUE,
    pval.coord = c(1,0.65),
    ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since positive PCR",
    ylab = "Survival probability",
    break.time.by = 10,
    xlim = c(0, STUDY_MAX_PCR_TO_DEATH - STUDY_MIN_TO_DEATH + 5)
  )}

# By groups + gender

plt.survival_ngrp_gen <- study.survival %>%
  mutate(age_group_narrow=factor(.$age_group_narrow, levels = levels(.$age_group_narrow), labels = sub("^(.*)", "Age group: \\1", levels(.$age_group_narrow)))) %>%
  {ggsurvplot_facet(
    short.panel.labs = TRUE,
    facet.by = c('age_group_narrow', 'gender'),
    fit = survfit(
      Surv(days_first_positive_result_to_death, ind_death) ~ case,
      data= .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    pval=TRUE,
    #pval.coord = c(1,0.65),
    #ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since positive PCR",
    ylab = "Survival probability",
    break.time.by = 10,
    labeller = c(c("gender: female"), c("a"))
  ) + theme(panel.spacing = unit(1, "lines"), plot.margin = margin(20, 20, 20, 20))}

plt.survival_wgrp_gen <- study.survival %>%
  mutate(age_group=factor(.$age_group, levels = levels(.$age_group), labels = sub("^(.*)", "Age group: \\1", levels(.$age_group)))) %>%
  {ggsurvplot_facet(
    short.panel.labs = TRUE,
    facet.by = c('age_group', 'gender'),
    fit = survfit(
      Surv(days_first_positive_result_to_death, ind_death) ~ case,
      data= .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    pval=TRUE,
    pval.coord = c(1,0.65),
    ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since positive PCR",
    ylab = "Survival probability",
    break.time.by = 10
  ) + theme(panel.spacing = unit(1, "lines"), plot.margin = margin(20, 20, 20, 20))}

# By groups

plt.survival_ngrp <- study.survival %>%
  mutate(age_group_narrow=factor(.$age_group_narrow, levels = levels(.$age_group_narrow), labels = sub("^(.*)", "Age group: \\1", levels(.$age_group_narrow)))) %>%
  {ggsurvplot_facet(
    short.panel.labs = TRUE,
    facet.by = c('age_group_narrow'),
    fit = survfit(
      Surv(days_first_positive_result_to_death, ind_death) ~ case,
      data= .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    pval=TRUE,
    # pval.coord = c(1,0.65),
    #ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since positive PCR",
    ylab = "Survival probability",
    break.time.by = 10
  ) + theme(panel.spacing = unit(1, "lines"), plot.margin = margin(20, 20, 20, 20))}

plt.survival_wgrp <- study.survival %>%
  mutate(age_group=factor(.$age_group, levels = levels(.$age_group), labels = sub("^(.*)", "Age group: \\1", levels(.$age_group)))) %>%
  {ggsurvplot_facet(
    # ncol = 1,
    short.panel.labs = TRUE,
    facet.by = c('age_group'),
    fit = survfit(
      Surv(days_first_positive_result_to_death, ind_death) ~ case,
      data= .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    pval=TRUE,
    pval.coord = c(1,0.65),
    ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since positive PCR",
    ylab = "Survival probability",
    break.time.by = 10
  ) + theme(panel.spacing = unit(1, "lines"), plot.margin = margin(20, 20, 20, 20))}

## Hospitalization kaplan meier ----
#.............................

# Make data for kaplan meier analysis to cutoff time
# Make sure nrow(study %>% filter(ind_hospitalized == '1' & is.na(new_arrivaldate_first))) == 0

study.hosp_km <- study %>%
  mutate(
    ind_hospitalized = as.numeric(case_when(
      days_first_positive_result_to_hospitalization > STUDY_MAX_PCR_TO_HOSP ~ '0',
      TRUE ~ as.character(ind_hospitalized)
    )),
    days_first_positive_result_to_hospitalization = case_when(
      is.na(days_first_positive_result_to_hospitalization) ~ STUDY_MAX_PCR_TO_HOSP,
      days_first_positive_result_to_hospitalization > STUDY_MAX_PCR_TO_HOSP ~ STUDY_MAX_PCR_TO_HOSP,
      TRUE ~ days_first_positive_result_to_hospitalization
    ),
    days_first_positive_result_to_hospitalization = days_first_positive_result_to_hospitalization - STUDY_MIN_PCR_TO_HOSP
  )

# General

plt.hosp_km <- study.hosp_km %>%
  {ggsurvplot(
    fit = survfit(
      Surv(days_first_positive_result_to_hospitalization, ind_hospitalized) ~ case,
      data= .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    risk.table = TRUE,
    risk.table.title = "No. at Risk",
    table.heigth=.45,
    tables.theme = theme_cleantable(),
    pval=TRUE,
    pval.coord = c(1,0.65),
    ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since 2-days before positive PCR",
    ylab = "Hospitalization probability",
    #break.time.by = 10,
    xlim = c(0, STUDY_MAX_PCR_TO_HOSP - STUDY_MIN_PCR_TO_HOSP + 1)
  )}

# Aged 55 and above

plt.hosp_km55p <- study.hosp_km %>%
  # Leave only patient Aged 55 and above
  filter(age >= 55) %>%
  {ggsurvplot(
    fit=survfit(
      Surv(days_first_positive_result_to_hospitalization, ind_hospitalized) ~ case,
      data=  .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    risk.table = TRUE,
    risk.table.title = "No. at Risk",
    table.heigth=.45,
    tables.theme = theme_cleantable(),
    pval=TRUE,
    pval.coord = c(1,0.65),
    ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since 2-days before positive PCR",
    ylab = "Hospitalization probability",
    #break.time.by = 10,
    xlim = c(0, STUDY_MAX_PCR_TO_HOSP - STUDY_MIN_PCR_TO_HOSP + 1)
  )}

# Aged 65 and above

plt.hosp_km65p <- study.hosp_km %>%
  # Leave only patient Aged 65 and above
  filter(age >= 65) %>%
  {ggsurvplot(
    fit=survfit(
      Surv(days_first_positive_result_to_hospitalization, ind_hospitalized) ~ case,
      data=  .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    risk.table = TRUE,
    risk.table.title = "No. at Risk",
    table.heigth=.45,
    tables.theme = theme_cleantable(),
    pval=TRUE,
    pval.coord = c(1,0.65),
    ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since 2-days before positive PCR",
    ylab = "Hospitalization probability",
    #break.time.by = 10,
    xlim = c(0, STUDY_MAX_PCR_TO_HOSP - STUDY_MIN_PCR_TO_HOSP + 1)
  )}

# Aged 55 - 65

plt.hosp_km55_65 <- study.hosp_km %>%
  # Leave only patient Aged 55-65 and above
  filter(between(age, 55, 64)) %>%
  {ggsurvplot(
    fit=survfit(
      Surv(days_first_positive_result_to_hospitalization, ind_hospitalized) ~ case,
      data=  .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    risk.table = TRUE,
    risk.table.title = "No. at Risk",
    table.heigth=.45,
    tables.theme = theme_cleantable(),
    pval=TRUE,
    pval.coord = c(1,0.65),
    ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since 2-days before positive PCR",
    ylab = "Hospitalization probability",
    #break.time.by = 10,
    xlim = c(0, STUDY_MAX_PCR_TO_HOSP - STUDY_MIN_PCR_TO_HOSP + 1)
  )}

# By groups + gender

plt.hosp_km_ngrp_gen <- study.hosp_km %>%
  mutate(age_group_narrow=factor(.$age_group_narrow, levels = levels(.$age_group_narrow), labels = sub("^(.*)", "Age group: \\1", levels(.$age_group_narrow)))) %>%
  {ggsurvplot_facet(
    short.panel.labs = TRUE,
    facet.by = c('age_group_narrow', 'gender'),
    fit = survfit(
      Surv(days_first_positive_result_to_hospitalization, ind_hospitalized) ~ case,
      data= .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    pval=TRUE,
    #pval.coord = c(1,0.65),
    #ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since 2-days before positive PCR",
    ylab = "Hospitalization probability",
    break.time.by = 10
  ) + theme(panel.spacing = unit(1, "lines"), plot.margin = margin(20, 20, 20, 20))}

plt.hosp_km_wgrp_gen <- study.hosp_km %>%
  mutate(age_group=factor(.$age_group, levels = levels(.$age_group), labels = sub("^(.*)", "Age group: \\1", levels(.$age_group)))) %>%
  {ggsurvplot_facet(
    short.panel.labs = TRUE,
    facet.by = c('age_group', 'gender'),
    fit = survfit(
      Surv(days_first_positive_result_to_hospitalization, ind_hospitalized) ~ case,
      data= .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    pval=TRUE,
    #pval.coord = c(1,0.65),
    #ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since 2-days before positive PCR",
    ylab = "Hospitalization probability",
    break.time.by = 10
  ) + theme(panel.spacing = unit(1, "lines"), plot.margin = margin(20, 20, 20, 20))}

# By groups

plt.hosp_km_ngrp <- study.hosp_km %>%
  mutate(age_group_narrow=factor(.$age_group_narrow, levels = levels(.$age_group_narrow), labels = sub("^(.*)", "Age group: \\1", levels(.$age_group_narrow)))) %>%
  {ggsurvplot_facet(
    short.panel.labs = TRUE,
    facet.by = c('age_group_narrow'),
    fit = survfit(
      Surv(days_first_positive_result_to_hospitalization, ind_hospitalized) ~ case,
      data= .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    pval=TRUE,
    #pval.coord = c(1,0.65),
    #ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since 2-days before positive PCR",
    ylab = "Hospitalization probability",
    break.time.by = 10
  ) + theme(panel.spacing = unit(1, "lines"), plot.margin = margin(20, 20, 20, 20))}

plt.hosp_km_wgrp <- study.hosp_km %>%
  mutate(age_group=factor(.$age_group, levels = levels(.$age_group), labels = sub("^(.*)", "Age group: \\1", levels(.$age_group)))) %>%
  {ggsurvplot_facet(
    # ncol = 1,
    short.panel.labs = TRUE,
    facet.by = c('age_group'),
    fit = survfit(
      Surv(days_first_positive_result_to_hospitalization, ind_hospitalized) ~ case,
      data= .,
      se.fit=TRUE,
      conf.int=.95
    ),
    data= .,
    palette=c("black","darkred"),
    legend.title = "",
    legend.labs = c("Unvaccinated", "Recently injected"),
    pval=TRUE,
    #pval.coord = c(1,0.65),
    #ylim=c(0.6,1.01),
    axes.offset = FALSE,
    xlab = "Days since 2-days before positive PCR",
    ylab = "Hospitalization probability",
    break.time.by = 10
  ) + theme(panel.spacing = unit(1, "lines"), plot.margin = margin(20, 20, 20, 20))}

## Hospitalization ----
#.............................

# Here we test the correlation between vaccination and the probability to hospitalize in corona dedicated division

anl.hosp <- bind_rows(
  
  # Narrow age group
  study %>%
    # Remove unknow hospitalization status
    filter(ind_hospitalized %in% c('0', '1')) %>%
    # Update ind_hospitalized as factor
    mutate(ind_hospitalized = as.factor(ind_hospitalized)) %>%
    mutate(description='Narrow age group', age_group=age_group_narrow) %>%
    group_by(description, age_group, gender) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.hosp.n=sum(case=='1' & ind_hospitalized=='1'),
      control.hosp.n=sum(case=='0' & ind_hospitalized=='1'),
      case.hosp.p=sum(case=='1' & ind_hospitalized=='1')/sum(case=='1'),
      control.hosp.p=sum(case=='0' & ind_hospitalized=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, as.character(ind_hospitalized))$p.value,
      fisher.test.or=fisher.test(case, as.character(ind_hospitalized))$estimate,
      fisher.test.conf.int=paste(fisher.test(case, as.character(ind_hospitalized))$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Wider age group
  study %>%
    # Remove unknow hospitalization status
    filter(ind_hospitalized %in% c('0', '1')) %>%
    # Update ind_hospitalized as factor
    mutate(ind_hospitalized = as.factor(ind_hospitalized)) %>%
    mutate(description='Wider age group') %>%
    group_by(description, age_group, gender) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.hosp.n=sum(case=='1' & ind_hospitalized=='1'),
      control.hosp.n=sum(case=='0' & ind_hospitalized=='1'),
      case.hosp.p=sum(case=='1' & ind_hospitalized=='1')/sum(case=='1'),
      control.hosp.p=sum(case=='0' & ind_hospitalized=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, as.character(ind_hospitalized))$p.value,
      fisher.test.or=fisher.test(case, as.character(ind_hospitalized))$estimate,
      fisher.test.conf.int=paste(fisher.test(case, as.character(ind_hospitalized))$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Wider age group no gender
  study %>%
    # Remove unknow hospitalization status
    filter(ind_hospitalized %in% c('0', '1')) %>%
    # Update ind_hospitalized as factor
    mutate(ind_hospitalized = as.factor(ind_hospitalized)) %>%
    mutate(description='Wider age group') %>%
    group_by(description, age_group) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.hosp.n=sum(case=='1' & ind_hospitalized=='1'),
      control.hosp.n=sum(case=='0' & ind_hospitalized=='1'),
      case.hosp.p=sum(case=='1' & ind_hospitalized=='1')/sum(case=='1'),
      control.hosp.p=sum(case=='0' & ind_hospitalized=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, as.character(ind_hospitalized))$p.value,
      fisher.test.or=fisher.test(case, as.character(ind_hospitalized))$estimate,
      fisher.test.conf.int=paste(fisher.test(case, as.character(ind_hospitalized))$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Aged 55 and above by gender
  study %>%
    # Remove unknow hospitalization status
    filter(ind_hospitalized %in% c('0', '1')) %>%
    # Update ind_hospitalized as factor
    mutate(ind_hospitalized = as.factor(ind_hospitalized)) %>%
    # Leave only patient Aged 55 and above
    filter(age >= 55) %>%
    mutate(description='Patient aged 55 and above') %>%
    group_by(description, gender) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.hosp.n=sum(case=='1' & ind_hospitalized=='1'),
      control.hosp.n=sum(case=='0' & ind_hospitalized=='1'),
      case.hosp.p=sum(case=='1' & ind_hospitalized=='1')/sum(case=='1'),
      control.hosp.p=sum(case=='0' & ind_hospitalized=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, as.character(ind_hospitalized))$p.value,
      fisher.test.or=fisher.test(case, as.character(ind_hospitalized))$estimate,
      fisher.test.conf.int=paste(fisher.test(case, as.character(ind_hospitalized))$conf.int, collapse = '-'),
      .groups = 'drop'),
  
  # Aged 55 and above
  study %>%
    # Remove unknow hospitalization status
    filter(ind_hospitalized %in% c('0', '1')) %>%
    # Update ind_hospitalized as factor
    mutate(ind_hospitalized = as.factor(ind_hospitalized)) %>%
    # Leave only patient Aged 55 and above
    filter(age >= 55) %>%
    mutate(description='Patient aged 55 and above') %>%
    group_by(description) %>% 
    dplyr::summarise(
      case.n=sum(case=='1'),
      control.n=sum(case=='0'),
      case.hosp.n=sum(case=='1' & ind_hospitalized=='1'),
      control.hosp.n=sum(case=='0' & ind_hospitalized=='1'),
      case.hosp.p=sum(case=='1' & ind_hospitalized=='1')/sum(case=='1'),
      control.hosp.p=sum(case=='0' & ind_hospitalized=='1')/sum(case=='0'),
      fisher.test.p.value=fisher.test(case, as.character(ind_hospitalized))$p.value,
      fisher.test.or=fisher.test(case, as.character(ind_hospitalized))$estimate,
      fisher.test.conf.int=paste(fisher.test(case, as.character(ind_hospitalized))$conf.int, collapse = '-'),
    .groups = 'drop')
)

anl.hosp$description <- as.factor(anl.hosp$description)

## Make parameters table ----
#.............................

anl.params <- data.frame(description=c(), value=c())

anl.params <- rbind(anl.params, data.frame(
  description='Study start',
  value=STUDY_START
))

anl.params <- rbind(anl.params, data.frame(
  description='Study end',
  value=STUDY_END
))

anl.params <- rbind(anl.params, data.frame(
  description='Study maximum positive pcr to hospitalization',
  value=STUDY_MAX_PCR_TO_HOSP
))

anl.params <- rbind(anl.params, data.frame(
  description='Study minimum positive pcr to hospitalization',
  value=STUDY_MIN_PCR_TO_HOSP
))

anl.params <- rbind(anl.params, data.frame(
  description='Study maximum positive pcr to death',
  value=STUDY_MAX_PCR_TO_DEATH
))

anl.params <- rbind(anl.params, data.frame(
  description='Study minimum days to death',
  value=STUDY_MIN_TO_DEATH
))

anl.params <- rbind(anl.params, data.frame(
  description='Study minimum age',
  value=STUDY_MIN_AGE
))

anl.params <- rbind(anl.params, data.frame(
  description='Rows in master source',
  value=nrow(source.master)
))

anl.params <- rbind(anl.params, data.frame(
  description='Rows in master without duplicate complete rows',
  value=nrow(unique.master)
))

anl.params <- rbind(anl.params, data.frame(
  description='Rows in master with complete non-confilicted cases',
  value=nrow(data.master %>% filter(err == ''))
))

anl.params <- anl.params %>% mutate(description=as.factor(description))

## Make demographic table ----
#.............................

anl.demographic <- data.frame(description=c(), general=c(), case=c(), control=c())

anl.demographic <- rbind(
  anl.demographic,
  data.comb %>%
    filter(between(first_positive_result_test_date, as.Date(STUDY_START), as.Date(CUTTOF_3T4WAVE) - 1)) %>%
    filter(STUDY_MIN_AGE <= age) %>%
    {data.frame(
      description=sprintf('Complete non-confilicted cases aged %d+ during %s - %s', STUDY_MIN_AGE, STUDY_START, as.Date(CUTTOF_3T4WAVE) - 1),
      general=sprintf('%d', nrow(.) ),
      case=sprintf('%d(%.2f%%)', nrow(filter(.,case == '1')), nrow(filter(.,case == '1')) / nrow(.) * 100),
      control=sprintf('%d(%.2f%%)', nrow(filter(.,case == '0')), nrow(filter(.,case == '0')) / nrow(.) * 100)
    )}
)

anl.demographic <- rbind(
  anl.demographic,
  data.comb %>%
    filter(between(first_positive_result_test_date, as.Date(CUTTOF_3T4WAVE), as.Date(STUDY_END))) %>%
    filter(STUDY_MIN_AGE <= age) %>%
    {data.frame(
      description=sprintf('Complete non-confilicted cases aged %d+ during %s - %s', STUDY_MIN_AGE, CUTTOF_3T4WAVE, STUDY_END),
      general=sprintf('%d', nrow(.) ),
      case=sprintf('%d(%.2f%%)', nrow(filter(.,case == '1')), nrow(filter(.,case == '1')) / nrow(.) * 100),
      control=sprintf('%d(%.2f%%)', nrow(filter(.,case == '0')), nrow(filter(.,case == '0')) / nrow(.) * 100)
    )}
)

anl.demographic <- rbind(
  anl.demographic,
  data.comb %>%
    filter(between(first_positive_result_test_date, as.Date(STUDY_START), as.Date(STUDY_END))) %>%
    filter(STUDY_MIN_AGE <= age) %>%
    {data.frame(
       description=sprintf('Complete non-confilicted cases aged %d+ during %s - %s', STUDY_MIN_AGE, STUDY_START, STUDY_END),
       general=sprintf('%d', nrow(.) ),
       case=sprintf('%d(%.2f%%)', nrow(filter(.,case == '1')), nrow(filter(.,case == '1')) / nrow(.) * 100),
       control=sprintf('%d(%.2f%%)', nrow(filter(.,case == '0')), nrow(filter(.,case == '0')) / nrow(.) * 100)
    )}
)

anl.demographic <- rbind(
  anl.demographic,
  study.pre.match %>%
    {data.frame(
      description='In the study group brfore match',
      general=sprintf(
        '%d',
        nrow(.)
      ),
      case=sprintf(
        '%d(%.2f%%)',
        nrow(filter(.,case == '1')),
        nrow(filter(.,case == '1')) / nrow(.) * 100
      ),
      control=sprintf(
        '%d(%.2f%%)',
        nrow(filter(.,case == '0')),
        nrow(filter(.,case == '0')) / nrow(.) * 100
      )
    )}
)

anl.demographic <- rbind(
  anl.demographic,
  study %>%
    {data.frame(
      description='In the study group after match (% of the same group before match)',
      general=sprintf(
        '%d',
        nrow(.)
      ),
      case=sprintf(
        '%d(%.2f%%)',
        nrow(filter(.,case == '1')),
        nrow(filter(.,case == '1')) / nrow(filter(study.pre.match,case == '1')) * 100
      ),
      control=sprintf(
        '%d(%.2f%%)',
        nrow(filter(.,case == '0')),
        nrow(filter(.,case == '0')) / nrow(filter(study.pre.match,case == '0')) * 100
      )
    )}
)

anl.demographic <- rbind(
  anl.demographic,
  study %>%
    {data.frame(
      description=sprintf('Study median age(lower quartile - upper quartile)'),
      general=sprintf(
        '%d(%.2f-%.2f)',
        median(.$age),
        quantile(.$age)[2],
        quantile(.$age)[4]
      ),
      case=sprintf(
        '%d(%.2f-%.2f)',
        median(filter(.,case == '1')$age),
        quantile(filter(.,case == '1')$age)[2],
        quantile(filter(.,case == '1')$age)[4]
      ),
      control=sprintf(
        '%d(%.2f-%.2f)',
        median(filter(.,case == '0')$age),
        quantile(filter(.,case == '0')$age)[2],
        quantile(filter(.,case == '0')$age)[4]
      )
    )}
)

anl.demographic <- rbind(
  anl.demographic,
  study %>%
    group_by(age_group) %>%
    dplyr::summarise(
      general=sprintf(
        '%d(%.2f%%)',
        n(),
        n() / nrow(.) * 100
      ),
      control=sprintf(
        '%d(%.2f%%)',
        sum(case=='0'),
        sum(case=='0') / sum(.$case=='0') * 100
      ),
      case=sprintf(
        '%d(%.2f%%)',
        sum(case=='1'),
        sum(case=='1') / sum(.$case=='1') * 100
      )
    ) %>%
    mutate(
      description=sprintf('Study age group %s', age_group)
    ) %>%
    dplyr::select(description, general, case, control)
)

anl.demographic <- rbind(
  anl.demographic,
  study %>%
    group_by(age_group_narrow) %>%
    dplyr::summarise(
      general=sprintf(
        '%d(%.2f%%)',
        n(),
        n() / nrow(.) * 100
      ),
      control=sprintf(
        '%d(%.2f%%)',
        sum(case=='0'),
        sum(case=='0') / sum(.$case=='0') * 100
      ),
      case=sprintf(
        '%d(%.2f%%)',
        sum(case=='1'),
        sum(case=='1') / sum(.$case=='1') * 100
      )
    ) %>%
    mutate(
      description=sprintf('Study age group %s', age_group_narrow)
    ) %>%
    dplyr::select(description, general, case, control)
)

anl.demographic <- rbind(
  anl.demographic,
  study %>%
    group_by(gender) %>%
    dplyr::summarise(
      general=sprintf(
        '%d(%.2f%%)',
        n(),
        n() / nrow(.) * 100
      ),
      control=sprintf(
        '%d(%.2f%%)',
        sum(case=='0'),
        sum(case=='0') / sum(.$case=='0') * 100
      ),
      case=sprintf(
        '%d(%.2f%%)',
        sum(case=='1'),
        sum(case=='1') / sum(.$case=='1') * 100
      )
    ) %>%
    mutate(
      description=sprintf('Study sex %s', gender)
    ) %>%
    dplyr::select(description, general, case, control)
)

anl.demographic <- anl.demographic %>% mutate(description=as.factor(description))

## Hospitalization duration models----
#.............................

# Here we test the correlation between vaccination and the hospitalization duration in corona dedicated division

mdl.hosp.poisson <- study %>%
  # Remove rows with NA
  drop_na(duration_hospitalization) %>%
  # Remove dead patients
  filter(ind_death != '1') %>%
  # Make model
  {glm(as.numeric(duration_hospitalization) ~ case + as.factor(age_group) + as.factor(gender),
            family="poisson",data=.)}

mdl.hosp.poisson.all <- study %>%
  # Remove rows with NA
  drop_na(duration_hospitalization) %>%
  # Remove dead patients
  filter(ind_death != '1') %>%
  # Make model
  {glm(as.numeric(duration_hospitalization) ~ case * as.factor(age_group) * as.factor(gender),
               family="poisson",data=.)}


# Negative Binomial model

mdl.hosp.np<- study %>%
  # Remove rows with NA
  drop_na(duration_hospitalization) %>%
  # Remove dead patients
  filter(ind_death != '1') %>%
  # Make model
  {glm.nb(as.numeric(duration_hospitalization) ~ case + as.factor(age_group) + as.factor(gender),
                  data=.)}

mdl.hosp.nb.all <- study %>%
  # Remove rows with NA
  drop_na(duration_hospitalization) %>%
  # Remove dead patients
  filter(ind_death != '1') %>%
  # Make model
  {glm.nb(as.numeric(duration_hospitalization) ~ case * as.factor(age_group) * as.factor(gender),
                  data=.)}

# MAKE DATA FOR EXPORT ----
#+++++++++++++++++++++++++++++

exp_df_death <- function(df){
  return(df %>%
    mutate(
      case = as.character(case_when(
        0 < case.death.n & case.death.n < MIN_N & MIN_N <= case.n ~ sprintf(
          '<%d/%d(<%.2f%%)',
          MIN_N,
          case.n,
          100 * MIN_N / case.n
        ),
        TRUE ~ sprintf(
          '%s/%s(%.2f%%)',
          case_when(0 < case.death.n & case.death.n < MIN_N ~ sprintf("<%d", MIN_N), TRUE ~ sprintf("%d", case.death.n)),
          case_when(0 < case.n & case.n < MIN_N ~ sprintf("<%d", MIN_N), TRUE ~ sprintf("%d", case.n)),
          100 * case.death.p
        )
      )),
      control = as.character(case_when(
        0 < control.death.n & control.death.n < MIN_N & MIN_N <= control.n ~ sprintf(
          '<%d/%d(<%.2f%%)',
          MIN_N,
          control.n,
          100 * MIN_N / control.n
        ),
        TRUE ~ sprintf(
          '%s/%s(%.2f%%)',
          case_when(0 < control.death.n & control.death.n < MIN_N ~ sprintf("<%d", MIN_N), TRUE ~ sprintf("%d", control.death.n)),
          case_when(0 < control.n & control.n < MIN_N ~ sprintf("<%d", MIN_N), TRUE ~ sprintf("%d", control.n)),
          100 * control.death.p
        )
      )),
      odds.ratio=sprintf(
        '%.2f(%.2f-%.2f)',
        fisher.test.or,
        as.numeric(str_split_fixed(.$fisher.test.conf.int, '-', Inf)[,1]),
        as.numeric(str_split_fixed(.$fisher.test.conf.int, '-', Inf)[,2])
      ),
      p.value=as.character(case_when(
        fisher.test.p.value < 0.001 ~ '<0.001',
        TRUE ~ sprintf('%.3f',fisher.test.p.value)
      ))
    ) %>%
    dplyr::select(-case.n, -case.death.n, -case.death.p, -control.n, -control.death.n, -control.death.p, -fisher.test.or, -fisher.test.p.value, -fisher.test.conf.int)
  )
}

exp_df_hosp <- function(df){
  return(df %>%
     mutate(
       case = as.character(case_when(
         0 < case.hosp.n & case.hosp.n < MIN_N & MIN_N <= case.n ~ sprintf(
           '<%d/%d(<%.2f%%)',
           MIN_N,
           case.n,
           100 * MIN_N / case.n
         ),
         TRUE ~ sprintf(
           '%s/%s(%.2f%%)',
           case_when(0 < case.hosp.n & case.hosp.n < MIN_N ~ sprintf("<%d", MIN_N), TRUE ~ sprintf("%d", case.hosp.n)),
           case_when(0 < case.n & case.n < MIN_N ~ sprintf("<%d", MIN_N), TRUE ~ sprintf("%d", case.n)),
           100 * case.hosp.p
         )
       )),
       control = as.character(case_when(
         0 < control.hosp.n & control.hosp.n < MIN_N & MIN_N <= control.n ~ sprintf(
           '<%d/%d(<%.2f%%)',
           MIN_N,
           control.n,
           100 * MIN_N / control.n
         ),
         TRUE ~ sprintf(
           '%s/%s(%.2f%%)',
           case_when(0 < control.hosp.n & control.hosp.n < MIN_N ~ sprintf("<%d", MIN_N), TRUE ~ sprintf("%d", control.hosp.n)),
           case_when(0 < control.n & control.n < MIN_N ~ sprintf("<%d", MIN_N), TRUE ~ sprintf("%d", control.n)),
           100 * control.hosp.p
         )
       )),
       odds.ratio=sprintf(
         '%.2f(%.2f-%.2f)',
         fisher.test.or,
         as.numeric(str_split_fixed(.$fisher.test.conf.int, '-', Inf)[,1]),
         as.numeric(str_split_fixed(.$fisher.test.conf.int, '-', Inf)[,2])
       ),
       p.value=as.character(case_when(
         fisher.test.p.value < 0.001 ~ '<0.001',
         TRUE ~ sprintf('%.3f',fisher.test.p.value)
       ))
     ) %>%
     dplyr::select(-case.n, -case.hosp.n, -case.hosp.p, -control.n, -control.hosp.n, -control.hosp.p, -fisher.test.or, -fisher.test.p.value, -fisher.test.conf.int)
  )
}

exp.death <- exp_df_death(anl.death)
exp.death2d <- exp_df_death(anl.death2d)
exp.death2t5d <- exp_df_death(anl.death2t5d)
exp.hosp <- exp_df_hosp(anl.hosp)
exp.params <- anl.params
exp.demographic <- anl.demographic

# SAVE FILES ----
#+++++++++++++++++++++++++++++

pdf(paste0(RESDIR, 'plt.covariate.pdf'))
print(plt.covariate, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.first_positive_result_test_month.pdf'))
print(plt.first_positive_result_test_month, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.days_first_positive_result_to_hospitalization.pdf'))
print(plt.days_first_positive_result_to_hospitalization, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.days_first_positive_result_to_hospitalization_study.pdf'))
print(plt.days_first_positive_result_to_hospitalization_study, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.survival.pdf'))
print(plt.survival, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.survival55p.pdf'))
print(plt.survival55p, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.survival65p.pdf'))
print(plt.survival65p, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.survival55_65.pdf'))
print(plt.survival55_65, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.survival_ngrp_gen.pdf'))
print(plt.survival_ngrp_gen, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.survival_wgrp_gen.pdf'))
print(plt.survival_wgrp_gen, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.survival_ngrp.pdf'))
print(plt.survival_ngrp, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.survival_wgrp.pdf'))
print(plt.survival_wgrp, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.hosp_km.pdf'))
print(plt.hosp_km, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.hosp_km55p.pdf'))
print(plt.hosp_km55p, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.hosp_km65p.pdf'))
print(plt.hosp_km65p, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.hosp_km55_65.pdf'))
print(plt.hosp_km55_65, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.hosp_km_ngrp_gen.pdf'))
print(plt.hosp_km_ngrp_gen, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.hosp_km_wgrp_gen.pdf'))
print(plt.hosp_km_wgrp_gen, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.hosp_km_ngrp.pdf'))
print(plt.hosp_km_ngrp, newpage = FALSE)
invisible(dev.off())

pdf(paste0(RESDIR, 'plt.hosp_km_wgrp.pdf'))
print(plt.hosp_km_wgrp, newpage = FALSE)
invisible(dev.off())

cat("", capture.output(summary(mdl.hosp.poisson)), file = paste0(RESDIR, 'mdl.hosp.poisson.txt'), sep = "\n")
cat("", capture.output(summary(mdl.hosp.poisson.all)), file = paste0(RESDIR, 'mdl.hosp.poisson.all.txt'), sep = "\n")
cat("", capture.output(summary(mdl.hosp.np)), file = paste0(RESDIR, 'mdl.hosp.np.txt'), sep = "\n")
cat("", capture.output(summary(mdl.hosp.nb.all)), file = paste0(RESDIR, 'mdl.hosp.nb.all.txt'), sep = "\n")

write.csv(exp.death, paste0(RESDIR, 'death.csv'), row.names = FALSE, na = "")
write.csv(exp.death2d, paste0(RESDIR, 'death2d.csv'), row.names = FALSE, na = "")
write.csv(exp.death2t5d, paste0(RESDIR, 'death2t5d.csv'), row.names = FALSE, na = "")
write.csv(exp.hosp, paste0(RESDIR, 'hosp.csv'), row.names = FALSE, na = "")
write.csv(exp.params, paste0(RESDIR, 'params.csv'), row.names = FALSE, na = "")
write.csv(exp.demographic, paste0(RESDIR, 'demographic.csv'), row.names = FALSE, na = "")

# Set output back to terminal
sink(NULL, type = "output")
sink(NULL, type = "message")
close(ofile)

# Put code in file
ofile <- file(paste0(RESDIR, 'code.R'), open = "wt")
sink(ofile, type = "output")
cat(readChar(
  rstudioapi::getSourceEditorContext()$path,
  file.info(rstudioapi::getSourceEditorContext()$path)$size)
)
sink(NULL, type = "output")
close(ofile)

# SHOW RESULTS ----
#+++++++++++++++++++++++++++++

plt.covariate
plt.first_positive_result_test_month
plt.days_first_positive_result_to_hospitalization
plt.days_first_positive_result_to_hospitalization_study
plt.survival
plt.survival55p
plt.survival65p
plt.survival55_65
plt.survival_ngrp_gen
plt.survival_wgrp_gen
plt.survival_ngrp
plt.survival_wgrp
plt.hosp_km
plt.hosp_km55p
plt.hosp_km65p
plt.hosp_km55_65
plt.hosp_km_ngrp_gen
plt.hosp_km_wgrp_gen
plt.hosp_km_ngrp
plt.hosp_km_wgrp

summary(mdl.hosp.poisson)
summary(mdl.hosp.poisson.all)
summary(mdl.hosp.np)
summary(mdl.hosp.nb.all)

View(anl.death)
View(anl.death2d)
View(anl.death2t5d)
View(anl.hosp)
View(anl.params)
View(anl.demographic)
