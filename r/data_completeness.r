
library(tidyverse)

# FUNCTIONS ----

#' Generate completeness for the daily data of participants who
#' are 28 days past their date of randomisation.
#'
#' Each participant should have:
#' - one baseline record
#' - one discharge record
#' - one daily record per day up until the date of their discharge
#' - one day 28 record
#'
#' @param ddat The current daily dataset
#' @param rand_ref The reference date for participant inclusion in completeness summaries
#' Default is to include everyone
#' @return A tibble of data completeness
generate_completeness_data <- function(ddat) {
  # Filter to those at least 28 days post randomisation
  out <- ddat %>%
    select(
      StudyPatientID,
      Country,
      PT_CountryName,
      Location,
      PT_LocationName,
      RandDate,
      DIS_DateOfDischarge,
      DIS_Outcome,
      DD_n,
      BAS_rec,
      DD_rec,
      DIS_rec,
      D28_rec,
      WTH_FU
    ) %>%
    unique() %>%
    arrange(RandDate) %>%
    mutate(
      DaysToDischarge = DIS_DateOfDischarge - RandDate + 1,
      DD_expect = pmin(28, DaysToDischarge),
      DD_complete = as.integer(DD_n == DD_expect),
      FU_complete = as.integer(DD_complete == 1 & DIS_rec == 1 & D28_rec == 1),
      FU_complete = if_else(is.na(FU_complete) | FU_complete == 0, 0, 1)
    )
  return(out)
}


#' Summarise data completeness
#'
#' @param cp_dat Output from `generate_completeness_data`.
#' @return A number of tibbles summarising data completeness which can be used
#' with kableExtra to generate tables fro interim reports.
summarise_completeness_data <- function(cp_dat) {
  overall <- cp_dat %>%
    summarise(
      `Enrolled` = n(),
      `Withdrew` = sum(WTH_FU == 1),
      `Expect follow-up` = n() - sum(WTH_FU == 1),
      # Ignore participants who withdrew consent for follow-up
      `Expected daily` = if_else(is.na(sum(DD_expect[WTH_FU == 0])), "Unknown", as.character(sum(DD_expect[WTH_FU == 0]))),
      `Baseline` = sprintf("%i (%.2f)", sum(BAS_rec), mean(BAS_rec)),
      `Daily` = sprintf("%i (%.2f)", sum(DD_n[WTH_FU == 0], na.rm = T), mean(if_else(is.na(DD_complete[WTH_FU == 0]), 0L, DD_complete[WTH_FU == 0]))),
      `Discharge` = sprintf("%i (%.2f)", sum(DIS_rec[WTH_FU == 0]), mean(DIS_rec[WTH_FU == 0])),
      `Day 28` = sprintf("%i (%.2f)", sum(D28_rec[WTH_FU == 0]), sum(D28_rec[WTH_FU == 0]) / (n() - sum(WTH_FU == 1))),
      `Complete follow-up` = sprintf("%i (%.2f)", sum(FU_complete[WTH_FU == 0]), mean(FU_complete[WTH_FU == 0])),
    )
  by_site <- cp_dat %>%
    mutate(Site = paste0(Country, ": ", Location)) %>%
    group_by(Site) %>%
    summarise(
      `Enrolled` = n(),
      `Withdrew` = sum(WTH_FU == 1),
      `Expect follow-up` = n() - sum(WTH_FU == 1),
      # Ignore participants who withdrew consent for follow-up
      `Expected daily` = if_else(is.na(sum(DD_expect[WTH_FU == 0])), "Unknown", as.character(sum(DD_expect[WTH_FU == 0]))),
      `Baseline` = sprintf("%i (%.2f)", sum(BAS_rec), mean(BAS_rec)),
      `Daily` = sprintf("%i (%.2f)", sum(DD_n[WTH_FU == 0], na.rm = T), mean(if_else(is.na(DD_complete[WTH_FU == 0]), 0L, DD_complete[WTH_FU == 0]))),
      `Discharge` = sprintf("%i (%.2f)", sum(DIS_rec[WTH_FU == 0]), mean(DIS_rec[WTH_FU == 0])),
      `Day 28` = sprintf("%i (%.2f)", sum(D28_rec[WTH_FU == 0]), sum(D28_rec[WTH_FU == 0]) / (n() - sum(WTH_FU == 1))),
      `Complete follow-up` = sprintf("%i (%.2f)", sum(FU_complete[WTH_FU == 0]), mean(FU_complete[WTH_FU == 0])),
    )
  combined <- bind_rows(by_site, overall %>% mutate(Site = "Total"))
  crosstab <- cp_dat %>%
    count(BAS_rec, WTH_FU, DD_complete, DIS_rec, D28_rec, FU_complete)

  return(list(overall = overall, by_site = by_site, combined = combined, crosstab = crosstab))
}
