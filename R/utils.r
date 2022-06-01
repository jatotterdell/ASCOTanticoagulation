#' @title read_all_no_daily
#' @description
#' Read in the derived dataset with everything except daily data
#' @references The full dataset
read_all_no_daily <- function() {
  readRDS(file.path(ANTICOAG_DATA, "all_data.rds"))
}


#' @title read_all_daily
#' @description
#' Read in the derived dataset with everything including daily data
#' @references The full dataset with daily
read_all_daily <- function() {
  readRDS(file.path(ANTICOAG_DATA, "all_daily_data.rds"))
}


#' @title creatinine_clearance
#' @description
#' Serum creatinine clearance by Cockcroft-Gault formula
#' https://www.eviq.org.au/clinical-resources/eviq-calculators/3200-creatinine-clearance-calculator
#' @param sex Either "Male" or "Female"
#' @param age Age in years
#' @param weight Weight in kg
#' @param creatinine Serum creatinine in umol/L
#' @return Creatinine clearance rate in mL/min.
creatinine_clearance <- function(sex, age, weight, creatinine) {
  out <- (140 - age) * weight / (0.814 * 72 * creatinine)
  out[sex == "Female"] <- 0.85 * out[sex == "Female"]
  return(out)
}


#' @title transform_coding
#' @description
#' Convert from one contrast coding to another
#' @param cod_from The coding to transform from
#' @param cod_to The coding to transform to
#' @return The transformation matrix
transform_coding <- function(cod_from, cod_to) {
  MASS::ginv(cod_to) %*% cod_from
}


intervention_labels <- function() {
  list(
    AAssignment = c(
      "A0" = "Not randomised to antiviral",
      "A1" = "No specific antiviral",
      "A2" = "Nafamostat"
    ),
    CAssignment = c(
      "C0" = "Not randomised to anticoagulation",
      "C1" = "Standard dose",
      "C2" = "Intermediate dose",
      "C3" = "Standard dose plus aspirin",
      "C4" = "Therapeutic dose"
    )
  )
}


get_intervention_dates <- function() {
  tribble(
    ~ Domain, ~ Intervention, ~ stdate, ~ endate,
    "Anti-coagulation", "C1", as.Date("2021-02-18"), as.Date("2022-04-08"),
    "Anti-coagulation", "C2", as.Date("2021-02-18"), as.Date("2022-04-08"),
    "Anti-coagulation", "C3", as.Date("2021-02-18"), as.Date("2021-09-10"),
    "Anti-coagulation", "C4", as.Date("2021-10-14"), as.Date("2022-04-08"),
    "Anti-viral", "A1", as.Date("2021-06-10"), Sys.Date(),
    "Anti-viral", "A2", as.Date("2021-06-10"), Sys.Date()
  ) %>%
    mutate(Intervention = labelled(
      Intervention,
      labels = c(
        C1 = "Standard dose", C2 = "Intermediate dose", C3 = "Standard dose + aspirin", C4 = "Therapeutic dose",
        A1 = "No specific antiviral", A2 = "Nafamostat"
      ),
      label = "Domain intervention"))
}


get_interim_dates <- function() {
  tribble(
    ~ meet_num, ~ meet_date,
    1, as.Date("2021-07-21"),
    2, as.Date("2021-09-15"),
    3, as.Date("2021-12-01"),
    4, as.Date("2022-02-22")
  )
}


intervention_strata <- function() {
  tribble(
    ~ stdate, ~ endate, ~ strata_int, ~ strata_lab,
    as.Date("2021-02-18"), as.Date("2021-06-09"), 1, "Add A",
    as.Date("2021-06-10"), as.Date("2021-09-10"), 2, "Drop C3",
    as.Date("2021-09-11"), as.Date("2021-10-13"), 3, "Add C4",
    as.Date("2021-10-14"), as.Date("2022-04-08"), 4, "Drop C",
  )
}

#' @title filter_fas_itt
#' @description
#' Filter the full dataset down to the FAS-ITT,
#' that is, those who were enrolled and did not
#' withdraw from follow-up.
#' @param dat Data from `read_all_no_daily`
filter_fas_itt <- function(dat) {
  filter(
    dat,
    ENR_rec == 1,
    WTH_FU == 0
  )
}


#' @title filter_acs_itt
#' @description
#' Filter the full dataset down to the ACS-ITT,
#' that is, those who were enrolled and did not
#' withdraw from follow-up and were randomised to
#' intervention in domain C.
#' @param dat Data from `read_all_no_daily`
filter_acs_itt <- function(dat) {
  filter(
    dat,
    ENR_rec == 1,
    ACS_ITT == 1,
    WTH_FU == 0
  )
}


transmute_model_cols <- function(dat) {
  site_counts <- dat %>%
    dplyr::count(Country = factor(Country, levels = c("IN", "AU", "NP", "NZ")), Location)
  merge_aus <- site_counts %>%
    filter(Country == "AU", n < 5) %>%
    pull(Location)
  merge_nz <- site_counts %>%
    filter(Country == "NZ", n < 5) %>%
    pull(Location)
  dat <- dat %>%
    mutate(
      ctry = factor(Country, levels = c("IN", "AU", "NP", "NZ")),
      # group sites with < 5 counts
      site = fct_collapse(
        factor(Location, levels = site_counts$Location),
        `AU other` = merge_aus,
        `NZ other` = merge_nz
      )
    )
  region_site <- dat %>%
    dplyr::count(ctry, site) %>%
    mutate(
      ctry_num = as.numeric(ctry),
      site_num = as.numeric(fct_inorder(site))
    )
  dat <- dat %>%
    left_join(region_site %>% select(-n), by = c("ctry", "site")) %>%
    transmute(
      StudyPatientID,
      Sex,
      AAssignment= droplevels(factor(
        AAssignment, levels = c("A1", "A0", "A2"))),
      CAssignment = droplevels(factor(
        CAssignment, levels = c("C1", "C0", "C2", "C3", "C4"))),
      RandDate,
      PO,
      AgeAtEntry,
      weight = BAS_Weight,
      weightgt120 = as.numeric(weight > 120),
      oxygen_sat = if_else(BAS_PeripheralOxygen < 10, NA_real_, BAS_PeripheralOxygen),
      dsfs = as.numeric(RandDate - EL_FirstSymptoms),
      dsfsgt7 = as.numeric(dsfs > 7),
      age_c = AgeAtEntry - mean(AgeAtEntry),
      agegte60,
      agegte60_c = agegte60 - mean(agegte60),
      country = factor(Country, levels = c("IN", "AU", "NP", "NZ")),
      inelgc3 = if_else(EL_inelg_c3 == 1 | is.na(EL_inelg_c3), 1, 0),
      ctry = factor(Country, levels = c("IN", "AU", "NP", "NZ")),
      ctry_num,
      # group sites with < 5 counts
      site = fct_collapse(
        factor(Location, levels = site_counts$Location),
        `AU other` = merge_aus,
        `NZ other` = merge_nz
      ),
      site_num,
      relRandDate = as.numeric(max(RandDate) - RandDate),
      epoch = floor(relRandDate / 28),
      # Manual merge after check of count(epoch)
      epoch = case_when(
        epoch %in% 0:1 ~ 2,
        epoch == 14 ~ 13,
        TRUE ~ epoch
      ) - 1,
    ) %>%
    group_by(epoch) %>%
    mutate(epoch_lab = paste(
      format(min(RandDate), "%d%b%y"),
      format(max(RandDate), "%d%b%y"),
      sep = "-")
    ) %>%
    ungroup()
  return(dat)
}

# Model helpers ----

logit <- function(x) log(x) - log(1 - x)

expit <- function(x) 1 / (1 + exp(-x))

ordered_logit <- function(x) {
  c(
    1 - expit(-x[1]),
    expit(-x[1:(length(x)-1)]) - expit(-x[2:(length(x))]),
    expit(-tail(x, 1))
  )
}

# Save outputs ----


#' @title save_tex_table
#' @description
#' Save a latex table
#' @param tab The latex table
#' @param fn The file name (without .tex extension)
#' @return Nothing, but writes the table to outputs/tables/<fn>.tex
save_tex_table <- function(tab, fn) {
  writeLines(tab, con = file.path("outputs", "tables", paste0(fn, ".tex")))
}



save_cmdstanr_model <- function(mod, fn) {
  sf <- paste0(fn, ".rds")
  mod$save_object(file.path("outputs", "models", sf))
}