read_all_no_daily <- function() {
  readRDS(file.path(ANTICOAG_DATA, "all_data.rds"))
}


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


filter_fas_itt <- function(dat) {
  filter(
    dat,
    ENR_rec == 1,
    WTH_FU == 0
  )
}


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
      AAssignment= factor(
        AAssignment, levels = c("A1", "A0", "A2")),
      CAssignment = factor(
        CAssignment, levels = c("C1", "C0", "C2", "C3", "C4")),
      RandDate,
      PO,
      age_c = AgeAtEntry - mean(AgeAtEntry),
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
      ) - 1
    )
  return(dat)
}

# Save outputs ----

save_tex_table <- function(tab, fn) {
  writeLines(tab, con = file.path("outputs", "tables", paste0(fn, ".tex")))
}

save_cmdstanr_model <- function(mod, fn) {
  sf <- paste0(fn, ".rds")
  mod$save_object(file.path("outputs", "models", sf))
}
