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
      "C1" = "Standard<br>dose",
      "C2" = "Intermediate<br>dose",
      "C3" = "Standard dose<br>plus aspirin",
      "C4" = "Therapeutic<br>dose"
    )
  )
}

intervention_labels2 <- function() {
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
    "Anticoagulation", "C1", as.Date("2021-02-18"), as.Date("2022-04-08"),
    "Anticoagulation", "C2", as.Date("2021-02-18"), as.Date("2022-04-08"),
    "Anticoagulation", "C3", as.Date("2021-02-18"), as.Date("2021-09-10"),
    "Anticoagulation", "C4", as.Date("2021-10-14"), as.Date("2022-04-08"),
    "Antiviral", "A1", as.Date("2021-06-10"), as.Date("2022-04-09"),
    "Antiviral", "A2", as.Date("2021-06-10"), as.Date("2022-04-09")
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

# Filter analysis sets ----

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

# Filter concurrent controls ----

filter_concurrent_intermediate <- function(dat) {
  dat %>%
    # Restrict to participants randomised to C1 or C2
    filter(CAssignment %in% c("C1", "C2")) %>%
    mutate(
      CAssignment = droplevels(CAssignment),
      randC  = droplevels(randC)
    )
}


filter_concurrent_std_aspirin <- function(dat) {
  dat %>%
    # Patients randomised before closure of C3
    filter(RandDate < get_intervention_dates()$endate[3]) %>%
    # Patients ineligible for aspirin arm
    filter(inelgc3 == 0) %>%
    mutate(CAssignment = droplevels(CAssignment),
           randC  = droplevels(randC),
           ctry = droplevels(ctry))
}


filter_concurrent_therapeutic <- function(dat) {
  dat %>%
    # Patients randomised after opening of C4
    filter(RandDate >= get_intervention_dates()$stdate[4] &
           RandDate < get_intervention_dates()$endate[4]) %>%
    filter(country %in% c("AU", "NP", "NZ") | site == "BCM") %>%
    mutate(CAssignment = droplevels(CAssignment),
           randC  = droplevels(randC),
           ctry = relevel(ctry, ref = "NP"))
}


# Transmute model columns ----

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
      randA = factor(AAssignment, levels = c("A1", "A2")),
      randC = factor(CAssignment, levels = c("C1", "C2", "C3", "C4")),
      PO,
      AgeAtEntry,
      aspirin = if_else(BAS_PatientTakingAspirin == "Yes", 1, 0),
      weight = BAS_Weight,
      weightgt120 = as.numeric(weight > 120),
      oxygen_sat = if_else(BAS_PeripheralOxygen < 10, NA_real_, BAS_PeripheralOxygen),
      dsfs = as.numeric(RandDate - EL_FirstSymptoms),
      dsfsgt7 = as.numeric(dsfs > 7),
      rec_steroids = if_else(DIS_ImmunoCorticosteroids == "Yes", 1, 0),
      age_c = AgeAtEntry - mean(AgeAtEntry),
      agegte60,
      agegte60_c = agegte60 - mean(agegte60),
      country = factor(Country, levels = c("IN", "AU", "NP", "NZ")),
      inelgc3 = if_else(EL_inelg_c3 == 0 | is.na(EL_inelg_c3), 0, 1),
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
      epoch_raw = floor(relRandDate / 28),
      # Manual merge after check of count(epoch)
      epoch = case_when(
        epoch_raw %in% 0:1 ~ 2,
        epoch_raw == 14 ~ 13,
        TRUE ~ epoch_raw
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


#' transmute_model_cols_grp_aus_nz
#'
#' Same as transmute_model_cols, but join Australia and new zealand
#' together as one "region" with nested sites
transmute_model_cols_grp_aus_nz <- function(dat) {
  site_counts <- dat %>%
    dplyr::count(
      Region = fct_collapse(
        factor(Country, levels = c("IN", "AU", "NP", "NZ")),
        "AU/NZ" = c("AU", "NZ")),
      Location)
  merge_ausnz <- site_counts %>%
    filter(Region == "AU/NZ", n < 5) %>%
    pull(Location)
  dat <- dat %>%
    mutate(
      ctry = fct_collapse(
        factor(Country, levels = c("IN", "AU", "NP", "NZ")),
        "AU/NZ" = c("AU", "NZ")),
      # group sites with < 5 counts
      site = fct_collapse(
        factor(Location, levels = site_counts$Location),
        `AU/NZ other` = merge_ausnz
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
        AAssignment, levels = c("A0", "A1", "A2"))),
      CAssignment = droplevels(factor(
        CAssignment, levels = c("C1", "C0", "C2", "C3", "C4"))),
      RandDate,
      randA = factor(AAssignment, levels = c("A1", "A2")),
      randC = factor(CAssignment, levels = c("C1", "C2", "C3", "C4")),
      PO,
      AgeAtEntry,
      aspirin = if_else(BAS_PatientTakingAspirin == "Yes", 1, 0),
      weight = BAS_Weight,
      weightgt120 = as.numeric(weight > 120),
      oxygen_sat = if_else(BAS_PeripheralOxygen < 10, NA_real_, BAS_PeripheralOxygen),
      dsfs = as.numeric(RandDate - EL_FirstSymptoms),
      dsfsgt7 = as.numeric(dsfs > 7),
      rec_steroids = if_else(DIS_ImmunoCorticosteroids == "Yes", 1, 0),
      age_c = AgeAtEntry - mean(AgeAtEntry),
      agegte60,
      agegte60_c = agegte60 - mean(agegte60),
      country = factor(Country, levels = c("IN", "AU", "NP", "NZ")),
      inelgc3 = if_else(EL_inelg_c3 == 0 | is.na(EL_inelg_c3), 0, 1),
      ctry = fct_collapse(
        factor(Country, levels = c("IN", "AU", "NP", "NZ")),
        "AU/NZ" = c("AU", "NZ")),
      ctry_num,
      # group sites with < 5 counts
      site = fct_collapse(
        factor(Location, levels = site_counts$Location),
        `AU/NZ other` = merge_ausnz
      ),
      site_num,
      relRandDate = as.numeric(max(RandDate) - RandDate),
      epoch_raw = floor(relRandDate / 28),
      # Manual merge after check of count(epoch)
      epoch = case_when(
        epoch_raw %in% 0:1 ~ 2,
        epoch_raw == 14 ~ 13,
        TRUE ~ epoch_raw
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

logit <- function(x) {
  log(x) - log(1 - x)
}

expit <- function(x) {
  1 / (1 + exp(-x))
}

ordered_logit <- function(x) {
  c(
    1 - expit(-x[1]),
    expit(-x[1:(length(x)-1)]) - expit(-x[2:(length(x))]),
    expit(-tail(x, 1))
  )
}


#' make_domA_design
#'
#' Create design matrix where domain has an intercept which
#' indicates whether a participant was randomised to it or not
#' and (by default) uses orthonormal contrasts for available
#' interventions.
#'
#' This keeps the "not randomised" group out of the overall
#' mean around which interventions are centred.
#'
#' @param dat A dataset with factor variable "randA" for which NA means not randomised.
#' @param ctr Contrast, by defualt orthonormal
#' @return Design matrix
make_domA_design <- function(dat, ctr = contr.orthonorm) {
  XA <- model.matrix(
    ~ randA,
    model.frame(~ randA, dat, na.action = na.pass),
    contrasts = list(randA = ctr)
  )
  # If not randomised to A
  XA[is.na(XA[, 2]), ] <- 0
  colnames(XA)[1] <- "randA"
  return(XA)
}


#' make_domC_design
#'
#' Create design matrix where domain has an intercept which
#' indicates whether a participant was randomised to it or not
#' and (by default) uses orthonormal contrasts for available
#' interventions.
#'
#' This keeps the "not randomised" group out of the overall
#' mean around which interventions are centred.
#'
#' @param dat A dataset with factor variable "randC" for which NA means not randomised.
#' @param ctr Contrast, by defualt orthonormal
#' @return Design matrix
make_domC_design <- function(dat, ctr = contr.orthonorm) {
  XC <- model.matrix(
    ~ randC,
    model.frame(~ randC, dat, na.action = na.pass),
    contrasts = list(randC = ctr)
  )
  # If not randomised to A
  XC[is.na(XC[, 2]), ] <- 0
  colnames(XC)[1] <- "randC"
  return(XC)
}


make_X_design <- function(
    dat,
    vars = NULL,
    ctr = contr.orthonorm,
    includeA = TRUE) {

  XA <- make_domA_design(dat, ctr)
  XC <- make_domC_design(dat, ctr)

  if (!is.null(vars)) {
    Xother <- model.matrix(
      as.formula(paste(" ~ ", paste(vars, collapse = " + "))),
      data = dat
    ) [, -1]
  } else {
    Xother <- NULL
  }

  if(includeA) {
    X <- cbind(XC, XA, Xother)
    attributes(X)$contrasts <- list(
      "randA" = attr(XA, "contrasts")$randA,
      "randC" = attr(XC, "contrasts")$randC
    )
  } else {
    X <- cbind(XC, Xother)
    attributes(X)$contrasts <- list(
      "randC" = attr(XC, "contrasts")$randC
    )
  }
  return(X)
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
