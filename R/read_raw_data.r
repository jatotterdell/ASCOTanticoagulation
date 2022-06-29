# Functions to read in and make minor fixes/formatting
# to the raw data tables.

# PACKAGES ----

library(tidyverse)

# FUNCTIONS ----

# read_dictionary <- function() {
#   read_csv(file.path(ASCOT_DATA, "dictionary", "2021-11-23", "ASCOT Dictionary Fields.csv"))
# }


read_dictionary <- function() {
  read_csv(list.files(file.path(ASCOT_DATA, "dictionary", "2022-05-17"), full.names = T)[-9])
}


read_eligibility_file <- function(fn) {
  eligibility <- read_csv(
    fn,
    na = c("", "NA", "n/a"),
    col_types = cols(
      .default = col_character(),
      EL_DOB = col_date("%d-%b-%Y"),
      EL_YOB = col_integer(),
      EL_AgeAtEntry = col_integer(),
      EL_AdmittedToHospital = col_date("%d-%b-%Y"),
      EL_Referral = col_date("%d-%b-%Y"),
      EL_Screened = col_date("%d-%b-%Y"),
      EL_BloodPlateletTestValue = col_double(),
      EL_BloodPlateletTestAs_x10_9_L = col_double(),
      EL_SerumCreatinineBlood = col_double(),
      EL_SerumCreatinine_umolL = col_double(),
      EL_SerumPotassium = col_double(),
      EL_SerumSodium = col_double(),
      EL_eGFR = col_double(),
      EL_FirstPositiveTest = col_date("%d-%b-%Y"),
      EL_FirstEnteredLocalTime = col_datetime("%d/%m/%Y %H:%M:%S"),
      EL_FirstEnteredUTC = col_datetime("%d/%m/%Y %H:%M:%S"),
      EL_EligibilityLastUpdated = col_datetime("%d/%m/%Y %H:%M:%S"),
      EL_EligibilityLastUpdatedUTC = col_datetime("%d/%m/%Y %H:%M:%S"),
      EL_Con_ConsentDate = col_datetime("%d/%m/%Y %H:%M:%S"),
      EL_PregnancyDueDate = col_date("%d-%b-%Y"),
      FirstSymptoms = col_date("%d-%b-%Y")
    )
  ) %>%
    rename(EL_FirstSymptoms = FirstSymptoms) %>%
    filter(EL_ProtocolVersion != "1.0") %>%
    mutate(EL_rec = 1)

  # Note that there appears to be a duplicated column: "EL_OralTherapeuticAnticoagAgents"
  # The two columns appear to be exact duplicates, so one is removed here
  eligibility <- eligibility %>%
    rename(EL_OralTherapeuticAnticoagAgents = EL_OralTherapeuticAnticoagAgents...38) %>%
    select(-EL_OralTherapeuticAnticoagAgents...50) %>%
    # As requested, exclude ineligible participant
    filter(StudyPatientID != "BLK00001" | is.na(StudyPatientID))
  return(eligibility)
}


read_enrolled_file <- function(fn) {
  enrolled <- read_csv(
    fn,
    na = c("", "NA", "n/a"),
    col_types = cols(
      .default = col_character(),
      AgeAtEntry = col_double(),
      PT_DOD = col_date("%d-%b-%Y"),
      PT_YOB = col_double()
    )
  ) %>%
    mutate(
      RandomisedLocal = parse_datetime(str_replace_all(
        RandomisedLocal, c("a.m." = "AM", "p.m." = "PM")
      ), "%d/%m/%Y %H:%M:%S %p"),
      DOB = parse_date(str_replace_all(DOB, c(" 12:00:00 a.m." = "", " 12:00:00 p.m." = "")), "%d/%m/%Y"),
      RandDate = as.Date(RandomisedLocal)
    ) %>%
    rename(StudyPatientID = StudyID) %>%
    mutate(ENR_rec = 1)
  return(enrolled)
}


read_consent_file <- function(fn) {
  consent <- read_csv(
    fn,
    na = c("", "NA", "n/a"),
    col_types = cols(
      .default = col_character(),
      CON_ConsentDate = col_datetime("%d/%m/%Y %H:%M:%S")
    )
  ) %>%
    filter(PT_ProtocolVersion != "1.0") %>%
    select(-PT_ProtocolVersion) %>%
    mutate(CON_rec = 1)
}


read_withdrawal_file <- function(fn) {
  withdrawal <- read_csv(
    fn,
    na = c("", "NA", "n/a"),
    col_types = cols(
      .default = col_character(),
      Con_WithdrawnDate = col_date("%d-%b-%Y"),
      CON_DatePatientIneligible = col_date("%d-%b-%Y")
    )
  ) %>%
    filter(PT_ProtocolVersion != "1.0") %>%
    # Fix some inconsistent naming capitalisation
    rename(
      CON_WithdrawnBy = Con_WithdrawnBy,
      CON_WithdrawnDate = Con_WithdrawnDate,
      CON_WithdrawalReason = Con_WithdrawalReason,
      CON_WithdrawalReason_Other = Con_WithdrawalReason_Other,
      CON_WithdrawalClincianReason = Con_WithdrawalClincianReason,
      CON_WithdrawnDomainA = Con_WithdrawnDomainA
    ) %>%
    select(-PT_ProtocolVersion) %>%
    mutate(WTH_rec = 1)
}


read_baseline_file <- function(fn) {
  baseline <- read_csv(
    fn,
    na = c("", "NA", "n/a"),
    col_types = cols(
      .default = col_character(),
      BAS_Weight = col_double(),
      BAS_VaccineDosesReceived = col_double(),
      BAS_DateLastVaccineDose = col_datetime("%d/%m/%Y %H:%M:%S"),
      BAS_PeripheralOxygen = col_double(),
      BAS_RespRateHighest = col_double(),
      BAS_UreaEntered = col_double(),
      BAS_UreaResult = col_double(),
      BAS_CRPEntered = col_double(),
      BAS_CRPResult = col_double(),
      BAS_DateRespiratoryTest = col_date("%d-%b-%Y"),
      BAS_CycleThresholdValue = col_double(),
      BAS_DDimerEntered = col_double(),
      BAS_DDimerResult = col_double(),
      BAS_APTT = col_double(),
      BAS_INR = col_double(),
      BAS_FibrinogenTest = col_double(),
      BAS_FibrinogenResult = col_double(),
      BAS_ProthrombinTime = col_double()
    )
  ) %>%
    filter(PT_ProtocolVersion != "1.0") %>%
    select(-PT_ProtocolVersion) %>%
    mutate(BAS_rec = 1)
  return(baseline)
}


read_discharge_file <- function(fn) {
  discharge <- read_csv(
    fn,
    na = c("", "NA", "n/a"),
    col_types = cols(
      .default = col_character(),
      DIS_DateOfDischarge = col_date("%d-%b-%Y"),
      DIS_DateOfDeath = col_date("%d-%b-%Y")
    )
  ) %>%
    filter(PT_ProtocolVersion != "1.0") %>%
    select(-PT_ProtocolVersion) %>%
    mutate(DIS_rec = 1)
  return(discharge)
}


read_daily_file <- function(fn) {
  daily <- read_csv(
    fn,
    na = c("", "NA", "n/a"),
    col_types = cols(
      .default = col_character(),
      DD_Date = col_date(format = "%d-%b-%Y"),
      DD_StudyDay = col_double(),
      DD_NafamostatDailyDose = col_double(),
      DD_NafamostatDuration = col_double(),
      DD_Potassium = col_double(),
      DD_Sodium = col_double(),
      DD_ALTLabs = col_double(),
      DD_ASTLabs = col_double(),
      DD_RespiratorySampleDate = col_date(format = "%d-%b-%Y")
    )
  ) %>%
    group_by(StudyPatientID) %>%
    mutate(DD_n = max(DD_StudyDay)) %>%
    ungroup() %>%
    filter(PT_ProtocolVersion != "1.0") %>%
    select(-PT_ProtocolVersion) %>%
    mutate(DD_rec = 1)
}


read_d28_file <- function(fn) {
  d28 <- read_csv(
    fn,
    na = c("", "NA", "n/a"),
    col_types = cols(
      .default = col_character(),
      D28_DateOfDeath = col_date("%d-%b-%Y"),
      D28_DateOfFollowUp = col_date("%d-%b-%Y"),
      D28_LastKnownDateAlive = col_date("%d-%b-%Y"),
      D28_OutcomeTotalDaysHospitalised = col_double(),
      D28_OutcomeDaysFreeOfVentilation = col_double(),
      D28_BreathScale = col_double(),
      D28_EQMobility = col_double(),
      D28_EQPersonalCare = col_double(),
      D28_EQUsualActivities = col_double(),
      D28_EQPainDiscomfort = col_double(),
      D28_EQAnxietyDepression = col_double(),
      D28_EQOverallHealthScore = col_double()
    )
  ) %>%
    filter(PT_ProtocolVersion != "1.0") %>%
    select(-PT_ProtocolVersion) %>%
    mutate(D28_rec = 1)
  return(d28)
}


read_d90_file <- function(fn) {
  d90 <- read_csv(
    fn,
    na = c("", "NA", "n/a"),
    col_types = cols(
      D90_DateOfFollowUp = col_date(format = "%d-%b-%Y"),
      PT_ProtocolVersion = col_character()
    )
  ) %>%
    filter(PT_ProtocolVersion != "1.0") %>%
    select(-PT_ProtocolVersion) %>%
    mutate(D90_rec = 1)
  return(d90)
}


read_deviation_file <- function(fn) {
  deviation <- read_csv(
    fn,
    na = c("", "NA", "n/a"),
  ) %>%
    filter(PT_ProtocolVersion != "1.0") %>%
    select(-PT_ProtocolVersion)
  return(deviation)
}


read_site_activation_file <- function() {
  active <- read_csv(file.path(
    ASCOT_DATA_RAW, "ASCOTADAPT_SiteActivationDates_20220620.csv"
  ),
  col_types = cols(
    `Site activated on protocol version 3.0` = col_date(format = "%d-%b-%y"),
    `Site activated on protocol version 5.0 (Therapeutic AC)` = col_date(format = "%d-%b-%y")
  )) %>%
    rename(
      state = STATE,
      site = `Site Code`,
      health_service = `Health service`,
      hospname_hospaddress = `Hospital name/address`,
      active_domains = `Domains Active`,
      active_v3 = `Site activated on protocol version 3.0`,
      active_v5 = `Site activated on protocol version 5.0 (Therapeutic AC)`
    ) %>%
    mutate(
      country = case_when(
        state == "India" ~ "IN",
        state == "Nepal" ~ "NP",
        state == "NZ" ~ "NZ",
        TRUE ~ "AU"
      )
    )
  return(active)
}


## Manually created files from trial management group ----


read_per_protocol_list <- function() {
  fn <- file.path(ASCOT_DATA_RAW, "Patients Not Per Protocol_31MAY2022.csv")
  out <- read_csv(fn) %>%
    mutate(Reason = str_replace(Reason, "Protol", "Protocol")) %>%
    rename(StudyPatientID = `Patient ID`)
  return(out)
}


read_protocol_deviations <- function() {
  fn <- file.path(ASCOT_DATA_RAW, "PDs_31MAY2022.csv")
  out <- read_csv(
    fn,
    col_types = cols(
      DateofDeviation = col_date(format = "%d/%m/%Y")
    ))
  return(out)
}


read_sars <- function() {
  fn <- file.path(ASCOT_DATA_RAW, "SARs_31MAY2022.csv")
  out <- read_csv(
    fn,
    na = c("", "NA", "n/a", "N/A"),
    col_types = cols(
      DateofSAEOnset = col_date(format = "%d/%m/%Y"),
      DateLoggedOnDatabase = col_date(format = "%d/%m/%Y"),
      ProtocolTreatmentStartDate = col_date(format = "%d/%m/%Y"),
      ProtocolTreatmentEndDate = col_date(format = "%d/%m/%Y"),
      DateOfResolution = col_date(format = "%d/%m/%Y")
    ))
  return(out)
}


read_saes <- function() {
  fn <- file.path(ASCOT_DATA_RAW, "SAEs_31MAY2022.csv")
  out <- read_csv(
    fn,
    na = c("", "NA", "n/a", "N/A"),
    col_types = cols(
      DateofSAEOnset = col_date(format = "%d/%m/%Y"),
      DateLoggedOnDatabase = col_date(format = "%d/%m/%Y"),
      ProtocolTreatmentStartDate = col_date(format = "%d/%m/%Y"),
      ProtocolTreatmentEndDate = col_date(format = "%d/%m/%Y"),
      DateofResolution = col_date(format = "%d/%m/%Y")
    ))
  return(out)
}


## Helpers to read all and save derived dataset ----

#' @title read_raw_extracts(fn)
#' @description
#' Read the raw data for extract "fn".
#'
#' @param fn The name of the extract as a string, or a vector of extract names
#' @return Returns nothing, but adds the extracted data to the global environment.
#' Always reads the most recent data (based on ordering of file name).
read_raw_extracts <- function(fn) {
  dirs <- file.path(ASCOT_DATA_RAW, fn)
  versions <- sapply(dirs, function(x) tail(list.files(x), 1))
  paths <- normalizePath(file.path(dirs, versions))
  for (i in 1:length(fn)) {
    assign(fn[i], get(paste0("read_", fn[i], "_file"))(paths[i]), envir = .GlobalEnv)
  }
}


read_all_raw_extracts <- function() {
  fns <- c("eligibility", "consent", "enrolled", "baseline", "withdrawal", "discharge", "daily", "d28")
  read_raw_extracts(fns)
}


#' @title save_derived_dataset(dat, fn)
#' @description
#' Save a derived dataset to the appropriate directory (ANTICOAG_DATA)
#'
#' @param dat The dataset
#' @param fn The file name which will be appended to `ANTICOAG_DATA`.
#' @return Nothing
save_derived_dataset <- function(dat, fn) {
  path <- file.path(ANTICOAG_DATA, fn)
  saveRDS(dat, path)
}
