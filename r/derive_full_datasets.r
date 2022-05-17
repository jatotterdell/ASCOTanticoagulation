# Functions to derive "full" datasets.
# Analysis scripts will filter/select on these for specific summaries/analyses.
# This should only be run "once" and then other scripts to use
# the full derived datasets.

# DEPS ----

source("r/read_raw_data.r")

# PACKAGES ----

library(tidyverse)
library(labelled)

# FUNCTIONS ----

findfirst <- function(x, v = NA) {
  j <- which(x)
  if(length(j)) min(j) else v
}

## CREATED DERIVED DATA ----

### FORMAT RAW DATA ----


#' @title add_database_corrections
#' @description
#' Some corrections to the data base are required, but cannot be made
#' to the database directly. These corrections are listed in
#' "ASCOT_ADAPT_DatabaseCorrectionsErrors.xlsx" and are here applied
#' manually to the data extracts.
#' @param dat Dataset containing relevant variables
#' @return `dat` but with fields corrected
add_database_corrections <- function(dat) {
  platID <- paste0("LUD000", formatC(1:22, width = 2, flag = "0"))
  creaID <- c(paste0("PUN000", formatC(1:32, width = 2, flag = "0")), paste0("PUN000", c(34, 55, 60, 66)))
  dat %>%
    mutate(
      D28_OutcomeDaysFreeOfVentilation = if_else(StudyPatientID == "MID00006", 28, D28_OutcomeDaysFreeOfVentilation),
      # Fix platelet units
      EL_BloodPlateletTestValueUnits = if_else(StudyPatientID %in% platID, "x 10<sup>9</sup>/L", EL_BloodPlateletTestValueUnits),
      EL_BloodPlateletTestAs_x10_9_L = if_else(StudyPatientID %in% platID, EL_BloodPlateletTestValue, EL_BloodPlateletTestAs_x10_9_L),
      # Fix wrong units for Creatinine units, note 1 mg/dL -> 88.42 umol/L
      EL_SerumCreatinineUnits = if_else(StudyPatientID %in% creaID, "mg/dL", EL_SerumCreatinineUnits),
      EL_SerumCreatinine_umolL = if_else(StudyPatientID %in% creaID, EL_SerumCreatinineBlood * 88.42, EL_SerumCreatinine_umolL)
    )
}


#' @title add_withdrawn_followup
#'
#' Add withdrawal flag to participant records.
#'
#' @param dat Dataset containing relevant withdrawal variables
#' @return Returns dat with a new indicator variable for whether
#' participant withdrew from follow-up.
add_withdrawn_followup <- function(dat) {
  dat %>%
    mutate(WTH_FU = if_else(
      WTH_rec == 0, 0L, as.integer(
        CON_WithdrawnContact28 == "Yes" |
          CON_WithdrawnDailyCollection == "Withdrawn" |
          # Covers withdrawal due to ineligibility
          is.na(CON_WithdrawnContact28)
      )
    ))
}


#' @title add_primary_outcome(dat, daily_dat)
#' @description Add the derived primary outcome for each participant.
#' @param dat Dataset with one record per participant with their outcomes
#' @param daily_dat Daily data set with multiple records per participant
#' @return Returns `dat` but with primary outcome fields appended.
add_primary_outcome <- function(dat, daily_dat) {

}


format_eligibility_data <- function(el) {
  el %>%
    select(-EligibilityID) %>%
    mutate(

      # Ineligible for Nafamostat
      EL_inelg_a2 = labelled(as.integer(
        EL_TherapeuticAnticoagBleeding == "Yes" |
          EL_HyperNafamostat == "Yes" |
          EL_ReceivedNafamostat == "Yes" |
          EL_HeartRenalDialysis == "Yes" |
          as.numeric(EL_SerumPotassium) > 5.5 |
          as.numeric(EL_SerumSodium) < 120
      ), label = "Ineligible for Nafamostat (A2)"),

      # Given only two interventions in domain, if
      # ineligible for A2 then ineligible for the domain
      EL_inelg_a = labelled(EL_inelg_a2, label = "Ineligible for domain A"),

      # Ineligible for the anticoagulation domain
      EL_inelg_c = labelled(as.integer(
        EL_TherapeuticAnticoag == "Yes" |
          EL_DualAntiplateletTherapy == "Yes" |
          EL_ContraHeparinReact == "Yes" |
          EL_IntracranialHaemorrhage == "Yes" |
          EL_BleedingConditionThrombo == "Yes" |
          as.numeric(EL_BloodPlateletTestAs_x10_9_L) < 30
      ), label = "Ineligible for domain C"),

      # Ineligible for the intervention C3
      # Note, once C3 had been removed, ceased to be assessed,
      # protocol v5 participants eligibility is unknown...
      EL_inelg_c3 = labelled(as.integer(
        EL_ReceivingAntiplatelet == "Yes" |
          EL_HyperAspirin == "Yes"
      ), label = "Ineligible for LWMH + Aspirin (C3)")
    )
}


format_consent_data <- function(con) {
  con %>%
    select(
      StudyPatientID,
      CON_Consent,
      CON_DomainA,
      CON_DomainB,
      CON_DomainC,
      CON_rec
    )
}


format_enrolled_data <- function(enr) {
  enr %>%
    select(-PT_YOB, -PT_DOD, -EID, -ECode) %>%
    mutate(
      ENR_regimen = paste0(AAssignment, BAssignment, CAssignment),
      FAS_ITT = 1L,
      ACS_ITT = if_else(CAssignment != "C0", 1L, 0L),
      AVS_ITT = if_else(AAssignment != "A0", 1L, 0L),
      agegte60 = labelled(as.integer(AgeAtEntry >= 60), label = "Age >= 60")
    )
}


format_baseline_data <- function(bas) {
  bas %>%
    select(-SDV, -FormLock, -BAS_PatientInitials) %>%
    # Summarise ethnicity data
    relocate(BAS_EthnicityUnknown, .before = "BAS_EthnicityAboriginal") %>%
    rowwise() %>%
    mutate(
      BAS_eth_all_na = all(is.na(c_across(BAS_EthnicityUnknown:BAS_EthnicityOther))),
      BAS_eth_all_but_unknown_na = all(is.na(c_across(BAS_EthnicityAboriginal:BAS_EthnicityOther))),
      BAS_eth_count = sum(c_across(BAS_EthnicityAboriginal:BAS_EthnicityOther) == "Yes", na.rm = T)
    ) %>%
    ungroup() %>%
    mutate(
      # Replace NA with "No"
      across(BAS_EthnicityUnknown:BAS_EthnicityOther, ~ if_else(is.na(.x), "No", .x))
    ) %>%
    # Summarise co-morbidities data
    relocate(BAS_Comorbidities_None, .before = "BAS_ChonicCardiacDisease") %>%
    rowwise() %>%
    mutate(
      BAS_com_all_na = all(is.na(c_across(BAS_Comorbidities_None:BAS_IatrogenicImmuno))),
      BAS_com_all_but_none_na = all(is.na(c_across(BAS_ChonicCardiacDisease:BAS_IatrogenicImmuno))),
      BAS_com_count = sum(c_across(BAS_ChonicCardiacDisease:BAS_IatrogenicImmuno) == "Yes", na.rm = T)
    ) %>%
    ungroup() %>%
    mutate(
      across(BAS_Comorbidities_None:BAS_IatrogenicImmuno, ~ if_else(is.na(.x), "No", .x))
    )
}


format_withdrawal_data <- function(wth) {
  wth %>%
    add_withdrawn_followup()
}


format_discharge_data <- function(dis) {
  dis %>%
    mutate(
      DIS_Outcome = factor(DIS_Outcome, levels = c(
        "Death",
        "Discharged alive (including hospital in the home)",
        "Discharged alive but against medical advice",
        "Unknown"
      )),
      DIS_Death = case_when(
        DIS_Outcome == "Death" ~ 1L,
        is.na(DIS_Outcome) | DIS_Outcome == "Unknown" ~ NA_integer_,
        TRUE ~ 0L
      ),
      DIS_DAMA = case_when(
        DIS_Outcome == "Discharged alive but against medical advice" ~ 1L,
        is.na(DIS_Outcome) | DIS_Outcome == "Unknown" ~ NA_integer_,
        TRUE ~ 0L
      ),
      DIS_DAMAlikelytodie = case_when(
        DIS_DAMA == 1 & DIS_LikelyToDie28 == "Yes" ~ 1L,
        DIS_DAMA == 1 & DIS_LikelyToDie28 == "No" ~ 0L,
        DIS_DAMA == 0 ~ 0L,
        is.na(DIS_DAMA) | DIS_DAMA == 1 & (is.na(DIS_LikelyToDie28) | DIS_LikelyToDie28 == "Unknown") ~ NA_integer_
      )
    ) %>%
    # Check antiviral and immouno use
    relocate(DIS_ReceivedNone, .before = "DIS_CamostatReceived") %>%
    relocate(Dis_ImmunoNone, .before = "DIS_ImmunoAnakinra") %>%
    rowwise() %>%
    mutate(
      DIS_av_all_na = all(is.na(c_across(DIS_ReceivedNone:DIS_ReceivedOther))),
      DIS_av_count = sum(c_across(DIS_CamostatReceived:DIS_ReceivedOther) == 1, na.rm = TRUE),
      DIS_im_all_na = all(is.na(c_across(Dis_ImmunoNone:DIS_IummunoOther))),
      DIS_im_count = sum(c_across(DIS_ImmunoAnakinra:DIS_IummunoOther) == 1, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      # Replace NA with "No"
      across(DIS_ReceivedNone:DIS_ReceivedOther, ~ if_else(is.na(.x), "No", .x)),
      across(Dis_ImmunoNone:DIS_IummunoOther, ~ if_else(is.na(.x), "No", .x))
    )
}


format_d28_data <- function(d28) {
  d28 %>%
    mutate(
      D28_who = as.integer(substr(D28_Status, 1, 1)),
      D28_who2 = as.integer(substr(D28_StatusVentilation, 1, 1)),
      # Required respiratory support AT D28
      D28_resp = case_when(
        D28_who < 6 ~ 0L,
        D28_who == 7 | (D28_who == 6 & D28_who2 %in% c(2, 4)) ~ 1L,
        is.na(D28_who) ~ NA_integer_,
        TRUE ~ NA_integer_
      ),
      # Days free of ventilation up to day 28
      D28_dfv = pmin(28, D28_OutcomeDaysFreeOfVentilation),
      # Days free of hospital up to day 28
      D28_dfh = pmax(0, 28 - D28_OutcomeTotalDaysHospitalised)
    )
}


#' @title format_daily_data
#' @description
#' Apply formatting to the raw daily data extract.
#' @param dd Raw daily data extract
format_daily_data <- function(dd) {
  dd %>%
    mutate(
      DD_who = as.integer(substr(DD_ParticipantDailyStatus, 1, 1)),
      DD_who2 = as.integer(substr(DD_O2, 1, 1))
    )
}


#' @title summarise_daily_data
#' @description
#' Summarise the daily data.
#' Requires fields from other baseline extracts, so
#' assumes all necessary variables are included in `dd`.
#' @param dd Daily data extract with additional variables as required.
summarise_daily_data <- function(dd) {
  dd %>%
    filter(DD_StudyDay <= 28 | is.na(DD_StudyDay)) %>%
    group_by(StudyPatientID) %>%
    summarise(
      DD_total_rec = n(),
      DD_total_days = max(DD_StudyDay),
      DD_who_missing = sum(is.na(DD_who)),
      DD_who_worst = max(DD_who, na.rm = TRUE),
      DD_who_best = min(DD_who, na.rm = TRUE),
      DD_resp = as.integer(any(DD_who2 %in% 1:2)),
      DD_vaso = as.integer(any(DD_O2VasopressorsInotropes == "Yes")),
      DD_death = as.integer(any(DD_who == 8)),
      DD_rec = as.integer(any(DD_who <= 3) | (DD_total_days < 28 & DD_who_worst < 8)),
      DD_ttr = if_else(DD_rec == 1, min(DD_total_days, findfirst(DD_who <= 3, v = 29)), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      DD_missing = pmin(28, DD_total_days) - DD_total_rec,
      DD_any_missing = if_else(is.na(DD_total_days) | DD_missing > 0, 1, 0),
    )
}


### JOIN DATASETS ADD FIELDS ----

#' @title create_fulldata_no_daily()
#' @description
#' Basically joins all the datasets together, except for the daily records.
#' Perform additional formatting here, in sub-steps for each raw data table.
#' Assumes all raw data extracts are in the global environment, e.g.
#' via `read_all_raw_extracts()`.
create_fulldata_no_daily <- function() {
  sub_eligibility <- format_eligibility_data(eligibility)
  sub_consent <- format_consent_data(consent)
  sub_enrolled <- format_enrolled_data(enrolled)
  sub_baseline <- format_baseline_data(baseline)
  sub_withdrawal <- format_withdrawal_data(withdrawal)
  sub_discharge <- format_discharge_data(discharge)
  sub_d28 <- format_d28_data(d28)

  sub_eligibility %>%
    left_join(sub_enrolled, by = "StudyPatientID") %>%
    left_join(sub_consent, by = "StudyPatientID") %>%
    left_join(sub_baseline, by = "StudyPatientID") %>%
    left_join(sub_withdrawal, by = "StudyPatientID") %>%
    left_join(sub_discharge, by = "StudyPatientID") %>%
    left_join(sub_d28, by = "StudyPatientID") %>%
    add_database_corrections() %>%
    mutate(
      ENR_rec = if_else(is.na(ENR_rec), 0, ENR_rec),
      CON_rec = if_else(is.na(CON_rec), 0, CON_rec),
      BAS_rec = if_else(is.na(BAS_rec), 0, BAS_rec),
      WTH_rec = if_else(is.na(WTH_rec), 0, WTH_rec),
      DIS_rec = if_else(is.na(DIS_rec), 0, DIS_rec),
      D28_rec = if_else(is.na(D28_rec), 0, D28_rec),
      WTH_FU = if_else(is.na(WTH_FU), 0L, WTH_FU),
      WTH_day = as.integer(CON_WithdrawnDate - RandDate + 1),
      DIS_day = as.numeric(DIS_DateOfDischarge - RandDate + 1),
      DIS_dday = as.numeric(DIS_DateOfDeath - RandDate + 1),
      DIS_deathlt28 = labelled(
        as.integer(DIS_Death == 1 & DIS_day <= 28),
        label = "Discharge: death within 28 days"
      ),
      D28_death = labelled(
        case_when(
          D28_PatientStatusDay28 == "Alive" ~ 0L,
          D28_PatientStatusDay28 == "Dead" ~ 1L,
          D28_PatientStatusDay28 == "Unknown" ~ NA_integer_,
          DIS_deathlt28 == 1 ~ 1L,
          TRUE ~ NA_integer_
        ),
        label = "All cause mortality at 28 days"
      )
    )
}


create_fulldata_add_daily <- function(dat) {
  dat %>%
    full_join(
      daily %>%
        format_daily_data(),
      by = "StudyPatientID"
    ) %>%
    mutate(
      ENR_rec = if_else(is.na(ENR_rec), 0, ENR_rec),
      CON_rec = if_else(is.na(CON_rec), 0, CON_rec),
      BAS_rec = if_else(is.na(BAS_rec), 0, BAS_rec),
      WTH_rec = if_else(is.na(WTH_rec), 0, WTH_rec),
      DIS_rec = if_else(is.na(DIS_rec), 0, DIS_rec),
      D28_rec = if_else(is.na(D28_rec), 0, D28_rec),
      DD_rec = if_else(is.na(DD_rec), 0, DD_rec)
    )
}

## SAVE DERIVED DATA ----

create_and_save_derived_data <- function() {
  read_all_raw_extracts()
  all_dat <- create_fulldata_no_daily()
  all_dat_daily <- create_fulldata_add_daily(all_dat)
  save_derived_dataset(all_dat, "all_data.rds")
  save_derived_dataset(all_dat_daily, "all_daily_data.rds")
}

## READ DERIVED DATA ----

read_all_no_daily <- function() {
  readRDS(file.path(ANTICOAG_DATA, "all_data.rds"))
}


read_all_daily <- function() {
  readRDS(file.path(ANTICOAG_DATA, "all_daily_data.rds"))
}
