# These are intended to be the "full" data with all relevant
# fields from all relevant data tables.
# Analysis scripts will filter/select on these for specific summaries/analyses.
# This should only be run "once" and then other scripts to use
# the full derived datasets.

# DEPS ----

source("r/read_raw_data.r")

# PACKAGES ----

library(tidyverse)

# FUNCTIONS ----


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


#' @title create_fulldata_no_daily()
#' @description
#' Basically joins all the datasets together, except for the daily records.
#' Perform additional formatting here, in sub-steps for each raw data table.
#' Assumes all raw data extracts are in the global environment, e.g.
#' via `read_all_raw_extracts()`.
create_fulldata_no_daily <- function() {

  sub_eligibility <- eligibility %>%
    select(-EligibilityID, -EligibilityCode)

  sub_enrolled <- enrolled %>%
    select(-PT_YOB, -PT_DOD, -EID, -ECode) %>%
    # Flag set membership
    mutate(
      FAS_ITT = 1L,
      ACS_ITT = if_else(CAssignment != "C0", 1L, 0L),
      AVS_ITT = if_else(AAssignment != "A0", 1L, 0L)
    )

  sub_baseline <- baseline %>%
    select(-SDV, -FormLock, -BAS_PatientInitials)

  sub_withdrawal <- withdrawal %>%
    add_withdrawn_followup()

  sub_discharge <- discharge

  sub_d28 <- d28

  sub_eligibility %>%
    left_join(sub_enrolled, by = "StudyPatientID") %>%
    left_join(sub_baseline, by = "StudyPatientID") %>%
    left_join(sub_withdrawal, by = "StudyPatientID") %>%
    left_join(sub_discharge, by = "StudyPatientID") %>%
    left_join(sub_d28, by = "StudyPatientID") %>%
    mutate(
      ENR_rec = if_else(is.na(ENR_rec), 0, ENR_rec),
      BAS_rec = if_else(is.na(BAS_rec), 0, BAS_rec),
      WTH_rec = if_else(is.na(WTH_rec), 0, WTH_rec),
      DIS_rec = if_else(is.na(DIS_rec), 0, DIS_rec),
      D28_rec = if_else(is.na(D28_rec), 0, D28_rec)
    ) %>%
    mutate(
      WTH_FU = if_else(is.na(WTH_FU), 0L, WTH_FU),
      WTH_Day = as.integer(CON_WithdrawnDate - RandDate)
    )
}


create_fulldata_add_daily <- function(dat) {
  dat %>%
    full_join(daily, by = "StudyPatientID") %>%
    mutate(
      ENR_rec = if_else(is.na(ENR_rec), 0, ENR_rec),
      BAS_rec = if_else(is.na(BAS_rec), 0, BAS_rec),
      WTH_rec = if_else(is.na(WTH_rec), 0, WTH_rec),
      DIS_rec = if_else(is.na(DIS_rec), 0, DIS_rec),
      D28_rec = if_else(is.na(D28_rec), 0, D28_rec)
    )
}


# PROCESSING ----

if (!interactive()) {
  read_all_raw_extracts()
  all_dat <- create_fulldata_no_daily()
  all_dat_daily <- create_fulldata_add_daily(all_dat)
  save_derived_dataset(all_dat, "all_data.rds")
  save_derived_dataset(all_dat_daily, "all_daily_data.rds")
}

#
