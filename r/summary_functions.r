library(tidyverse)
library(kableExtra)

# Baseline summary - demographics  ----

#' Generate baseline demographics summary by a grouping variable
#'
#' @param data The dataset. assumed to have baseline data included
#' @param grpvar The grouping variable
generate_baseline_demographics_by <- function(dat, grpvar = NULL) {
  grpvar <- enquo(grpvar)
  no_baseline <- sum(dat$BAS_rec == 0)
  tab <- dat %>%
    group_by(!!grpvar) %>%
    summarise(

      `Age (years), Median (IQR)` = sprintf("%.0f (%.0f, %.0f)", median(AgeAtEntry), quantile(AgeAtEntry, 0.25), quantile(AgeAtEntry, 0.75)),

      `Sex_Male, n (\\%)` = sprintf("%i (%.0f)", sum(Sex == "Male", na.rm = TRUE), 100 * sum(Sex == "Male", na.rm = TRUE) / n()),
      `Sex_Female, n (\\%)` = sprintf("%i (%.0f)", sum(Sex == "Female", na.rm = TRUE), 100 * sum(Sex == "Female", na.rm = TRUE) / n()),

      `Weight_Median, (IQR)` = sprintf("%.0f (%.0f, %.0f)", median(BAS_Weight, na.rm = T), quantile(BAS_Weight, 0.25, na.rm = T), quantile(BAS_Weight, 0.75, na.rm = T)),
      `Weight_Missing, n (\\%)` = sprintf("%i (%.0f)", sum(is.na(BAS_Weight)), 100 * sum(is.na(BAS_Weight)) / n()),

      `Vaccinated_Yes, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_PatientVaccinated == "Yes", na.rm = TRUE), 100 * sum(BAS_PatientVaccinated == "Yes", na.rm = TRUE) / n()),
      `Vaccinated_Missing, n (\\%)` = sprintf("%i (%.0f)", sum(is.na(BAS_PatientVaccinated)), 100 * sum(is.na(BAS_PatientVaccinated)) / n()),

      `Ethnicity_Indian, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityIndian == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityIndian == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Asian, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityAsian == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityAsian == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_European, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityEuropean == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityEuropean == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Maori, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityMaori == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityMaori == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Pacific Islander, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityPacificIslander == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityPacificIslander == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_African, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityAfrican == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityAfrican == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Aboriginal, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityAboriginal == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityAboriginal == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Latin American, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityLatinAmerican == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityLatinAmerican == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Middle Eastern, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityMiddleEastern == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityMiddleEastern == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Other, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityOther == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityOther == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Unknown, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityUnknown == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityUnknown == "Yes", na.rm = TRUE) / n()),

      `Smoking_Current, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_Smoking == "Current", na.rm = TRUE), 100 * sum(BAS_Smoking == "Current", na.rm = TRUE) / n()),
      `Smoking_Former, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_Smoking == "Former", na.rm = TRUE), 100 * sum(BAS_Smoking == "Former", na.rm = TRUE) / n()),
      `Smoking_Never, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_Smoking == "Never", na.rm = TRUE), 100 * sum(BAS_Smoking == "Never", na.rm = TRUE) / n()),
      `Smoking_Missing, n (\\%)` = sprintf("%i (%.0f)", sum(is.na(BAS_Smoking)), 100 * sum(is.na(BAS_Smoking)) / n())

    ) %>%
    gather(Variable, value, -!!grpvar, factor_key = TRUE) %>%
    spread(!!grpvar, value)
  colnames(tab)[-1] <- dat %>%
    count(!!grpvar) %>%
    mutate(lab = paste0(!!grpvar, "<br>(n = ", n, ")")) %>%
    pull(lab)
  return(tab)
}


#' Generate baseline demographics summary overall
#'
#' @param data The dataset. assumed to have baseline data included
generate_baseline_demographics <- function(dat) {
  no_baseline <- sum(dat$BAS_rec == 0)
  overall <- dat %>%
    summarise(
      `Age (years), Median (IQR)` = sprintf("%.0f (%.0f, %.0f)", median(AgeAtEntry), quantile(AgeAtEntry, 0.25), quantile(AgeAtEntry, 0.75)),

      `Sex_Male, n (\\%)` = sprintf("%i (%.0f)", sum(Sex == "Male", na.rm = TRUE), 100 * sum(Sex == "Male", na.rm = TRUE) / n()),
      `Sex_Female, n (\\%)` = sprintf("%i (%.0f)", sum(Sex == "Female", na.rm = TRUE), 100 * sum(Sex == "Female", na.rm = TRUE) / n()),

      `Weight_Median, (IQR)` = sprintf("%.0f (%.0f, %.0f)", median(BAS_Weight, na.rm = T), quantile(BAS_Weight, 0.25, na.rm = T), quantile(BAS_Weight, 0.75, na.rm = T)),
      `Weight_Missing, n (\\%)` = sprintf("%i (%.0f)", sum(is.na(BAS_Weight)), 100 * sum(is.na(BAS_Weight)) / n()),

      `Vaccinated_Yes, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_PatientVaccinated == "Yes", na.rm = TRUE), 100 * sum(BAS_PatientVaccinated == "Yes", na.rm = TRUE) / n()),
      `Vaccinated_Missing, n (\\%)` = sprintf("%i (%.0f)", sum(is.na(BAS_PatientVaccinated)), 100 * sum(is.na(BAS_PatientVaccinated)) / n()),

      `Ethnicity_Indian, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityIndian == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityIndian == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Asian, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityAsian == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityAsian == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_European, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityEuropean == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityEuropean == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Maori, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityMaori == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityMaori == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Pacific Islander, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityPacificIslander == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityPacificIslander == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_African, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityAfrican == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityAfrican == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Aboriginal, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityAboriginal == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityAboriginal == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Latin American, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityLatinAmerican == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityLatinAmerican == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Middle Eastern, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityMiddleEastern == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityMiddleEastern == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Other, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityOther == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityOther == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Unknown, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityUnknown == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityUnknown == "Yes", na.rm = TRUE) / n()),

      `Smoking_Current, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_Smoking == "Current", na.rm = TRUE), 100 * sum(BAS_Smoking == "Current", na.rm = TRUE) / n()),
      `Smoking_Former, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_Smoking == "Former", na.rm = TRUE), 100 * sum(BAS_Smoking == "Former", na.rm = TRUE) / n()),
      `Smoking_Never, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_Smoking == "Never", na.rm = TRUE), 100 * sum(BAS_Smoking == "Never", na.rm = TRUE) / n()),
      `Smoking_Missing, n (\\%)` = sprintf("%i (%.0f)", sum(is.na(BAS_Smoking)), 100 * sum(is.na(BAS_Smoking)) / n())

    ) %>%
    gather(Variable, Overall, factor_key = TRUE)
  colnames(overall)[2] <- paste0("Overall<br>(n = ", nrow(dat), ")")
  return(overall)
}


#' Generate baseline demographics tables
#'
#' @param data The dataset. assumed to have baseline data included
#' @param closed If FALSE, calculate overall only for open report.
#' If FALSE produce one table per domain.
#' @return Either a single table or a list of tables, one for each domain..
generate_baseline_demographics_table <- function(dat, format = "html") {
    byCgrp <- generate_baseline_demographics_by(dat %>% filter(CAssignment != "C0"), CAssignment)
    ovrC   <- generate_baseline_demographics(dat %>% filter(CAssignment != "C0"))
    tabC <- left_join(byCgrp, ovrC, by = "Variable") %>%
      mutate(Variable = str_replace(Variable, "[A-z]*_", ""))
    if(format == "latex") {
      colnames(tabC) <- linebreak(colnames(tabC), linebreaker = "<br>", align = "c")
    }
    outC <- kable(
      tabC,
      format = format,
      booktabs = T,
      caption = "Baseline demographics for participants randomised into domain C.",
      escape = F,
      align = "lrrrrr") %>%
      kable_styling(
        bootstrap_options = "striped",
        font_size = 11,
        latex_options = "HOLD_position") %>%
      group_rows("Sex", 2, 3) %>%
      group_rows("Weight (kg)", 4, 5) %>%
      group_rows("Vaccinated<sup>1</sup>", 6, 7, escape = FALSE) %>%
      group_rows("Ethnicity", 8, 18) %>%
      group_rows("Smoking", 19, 22) %>%
      add_header_above(c(" " = 1, "Anticoagulation" = ncol(byCgrp) - 1, " " = 1))  %>%
      row_spec(0, align = "c") %>%
      add_footnote("Site LUD does not have ethics approval for collection of vaccination status and accounts for most missingness", notation = "number")
    return(outC)
}

# Baseline summary - co-morbidities  ----

#' Generates baseline co-morbidity summary across all participants
#' @param dat
#' @return A tibble giving the summary
generate_baseline_comorbidities <- function(dat) {
  basdat <- dat %>%
    select(
      BAS_rec,
      BAS_Comorbidities_None,
      BAS_ChonicCardiacDisease,
      BAS_Hypertension,
      BAS_Obesity,
      BAS_ChronicLungDisease,
      BAS_ObsSleepApnoea,
      BAS_Asthma,
      BAS_Diabetes,
      BAS_ChronicKidneyDisease,
      BAS_Dialysis,
      BAS_ModSevLiverDisease,
      BAS_Dementia,
      BAS_MalignantNeoplasm,
      BAS_HIVInfection,
      BAS_IatrogenicImmuno)
  ovr_missing <- basdat %>%
    summarise(name = "Missing, n (\\%)", value = sprintf("%i (%3.1f)", sum(BAS_rec == 0), 100 * sum(BAS_rec == 0) / n()))
  overall <- basdat %>%
    select(-BAS_rec) %>%
    summarise_all(., list(`.n` = ~ sum(.x == "Yes", na.rm = T), `.p` = ~ sum(.x == "Yes", na.rm = T) / n())) %>%
    pivot_longer(everything(), names_to = c("name", "measure"), names_sep = "_\\.") %>%
    pivot_wider(id_cols = name, names_from = measure, values_from = value) %>%
    arrange(-n) %>%
    mutate(value = sprintf("%i (%3.1f)", n, 100 * p)) %>%
    select(-n, -p) %>%
    mutate(name = str_replace(name, "BAS_", ""),
           name = str_replace(name, "_Variable", ""),
           name = str_replace(name, "_", ""),
           name = str_replace(name, "ObsS", "ObstructiveS"),
           name = str_replace(name, "Immuno", "Immunosuppression"),
           name = str_replace(name, "ModSev", "ModerateOrSevere"),
           name = str_replace(name, "Comorbidities", ""),
           name = fct_inorder(gsub("([[:upper:]]*)([[:upper:]][[:lower:]]+)", "\\1 \\2", name)),
           name = trimws(stringr::str_to_sentence(name)),
           name = str_replace(name, "Hiv", "HIV"),
           name = str_replace(name, "Chonic", "Chronic"),
           name = paste0(name, ", n (\\%)")) %>%
    add_row(ovr_missing) %>%
    mutate(name = fct_inorder(name))
  colnames(overall)[1] <- "Comorbidity"
  colnames(overall) <- c("Comorbidity", paste0("Overall<br>(n = ", nrow(basdat), ")"))
  return(overall)
}


#' Generates baseline co-morbidity summary across grouping variable
#' @param dat The dataset
#' @param grpvar The grouping variable
#' @return A tibble giving the summary
generate_baseline_comorbidities_by <- function(dat, grpvar = NULL) {
  grpvar <- enquo(grpvar)
  basdat <- dat %>%
    group_by(!!grpvar) %>%
    select(
      BAS_rec,
      !!grpvar,
      BAS_Comorbidities_None,
      BAS_ChonicCardiacDisease,
      BAS_Hypertension,
      BAS_Obesity,
      BAS_ChronicLungDisease,
      BAS_ObsSleepApnoea,
      BAS_Asthma,
      BAS_Diabetes,
      BAS_ChronicKidneyDisease,
      BAS_Dialysis,
      BAS_ModSevLiverDisease,
      BAS_Dementia,
      BAS_MalignantNeoplasm,
      BAS_HIVInfection,
      BAS_IatrogenicImmuno)
  missing <- basdat %>%
    group_by(!!grpvar) %>%
    summarise(name = "Missing, n (\\%)", value = sprintf("%i (%3.1f)", sum(BAS_rec == 0), 100 * sum(BAS_rec == 0) / n())) %>%
    spread(!!grpvar, value)
  tab <- basdat %>%
    select(-BAS_rec) %>%
    summarise_all(., list(`.n` = ~ sum(.x == "Yes", na.rm = T),
                          `.p` = ~ sum(.x == "Yes", na.rm = T) / n())) %>%
    pivot_longer(-!!grpvar, names_to = c("name", "measure"), names_sep = "_\\.") %>%
    pivot_wider(id_cols = !!grpvar:name, names_from = measure, values_from = value) %>%
    mutate(value = sprintf("%i (%3.1f)", n, 100 * p)) %>%
    select(-n, -p) %>%
    spread(!!grpvar, value) %>%
    mutate(name = str_replace(name, "BAS_", ""),
           name = str_replace(name, "_Variable", ""),
           name = str_replace(name, "_", ""),
           name = str_replace(name, "ObsS", "ObstructiveS"),
           name = str_replace(name, "Immuno", "Immunosuppression"),
           name = str_replace(name, "ModSev", "ModerateOrSevere"),
           name = str_replace(name, "Comorbidities", ""),
           name = fct_inorder(gsub("([[:upper:]]*)([[:upper:]][[:lower:]]+)", "\\1 \\2", name)),
           name = trimws(stringr::str_to_sentence(name)),
           name = str_replace(name, "Hiv", "HIV"),
           name = str_replace(name, "Chonic", "Chronic"),
           name = paste0(name, ", n (\\%)")) %>%
    add_row(missing)
  colnames(tab)[1] <- "Comorbidity"
  colnames(tab) <- c(
    "Comorbidity",
    basdat %>%
      count(!!grpvar) %>%
      mutate(lab = paste0(!!grpvar, "<br>(n = ", n, ")")) %>%
      pull(lab))
  return(tab)
}



#' Generates baseline co-morbidity summary tables for each domain
#' @param dat The dataset
#' @param closed If TRUE, one table per domain, otherwise just aggregated table
#' @return A list of tibbles giving the summary tables
generate_baseline_comorbidities_table <- function(dat, format = "html") {
  ovr   <- generate_baseline_comorbidities(dat)
  ovrC   <- generate_baseline_comorbidities(dat %>% filter(CAssignment != "C0"))
  byAgrp <- generate_baseline_comorbidities_by(dat %>% filter(AAssignment != "A0"), AAssignment)
  byCgrp <- generate_baseline_comorbidities_by(dat %>% filter(CAssignment != "C0"), CAssignment)
  tabC <- left_join(ovrC, byCgrp, by = "Comorbidity")[, c(1,(2 + 1:(ncol(byCgrp) - 1)), 2)]
  if(format == "latex") {
    colnames(tabC) <- linebreak(colnames(tabC), linebreaker = "<br>", align = "c")
  }
  outC <- kable(
    tabC,
    format = format,
    booktabs = T,
    linesep = "",
    caption = "Baseline co-morbidities for participants randomised into domain C.",
    align = "lrrrrrrrr",
    escape = F) %>%
    kable_styling(
      bootstrap_options = "striped",
      font_size = 11,
      latex_options = "HOLD_position") %>%
    add_header_above(c(" " = 1, "Anticoagulation" = ncol(byCgrp) - 1, " " = 1)) %>%
    row_spec(0, align = "c")
  return(outC)
}


# Baseline summary - prognostics  ----


#' Generates baseline prognostics summary across grouping variable
#' @param dat The dataset
#' @param grpvar The grouping variable
#' @return A tibble giving the summary
generate_baseline_prognostics_by <- function(dat, grpvar = NULL) {
  grpvar <- enquo(grpvar)
  tab <- dat %>%
    group_by(!!grpvar) %>%
    summarise(

      "On Room Air_Yes, n (\\%)" = sprintf("%i (%.0f)",
                                           sum(BAS_OnRoomAir24hrs == "Yes", na.rm = TRUE),
                                           100 * sum(BAS_OnRoomAir24hrs == "Yes", na.rm = TRUE) / n()),
      "On Room Air_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                               sum(is.na(BAS_OnRoomAir24hrs)),
                                               100 * sum(is.na(BAS_OnRoomAir24hrs)) / n()),

      "GCS < 15_Yes, n (\\%)" = sprintf("%i (%.0f)",
                                        sum(BAS_PatientGCS == "Yes", na.rm = TRUE),
                                        100*sum(BAS_PatientGCS == "Yes", na.rm = TRUE)/n()),
      "GCS < 15_Missing, n (\\%)" = sprintf("%i (%.0f)", sum(is.na(BAS_PatientGCS)), 100*sum(is.na(BAS_PatientGCS)) / n()),

      "Peripheral oxygen saturation (SpO2)_Median (IQR)" = sprintf("%.0f (%.0f, %.0f)",
                                                                   median(BAS_PeripheralOxygen, na.rm = TRUE),
                                                                   quantile(BAS_PeripheralOxygen, na.rm = TRUE, prob = 0.25),
                                                                   quantile(BAS_PeripheralOxygen, na.rm = TRUE, prob = 0.75)),
      "Peripheral oxygen saturation (SpO2)_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                                                       sum(is.na(BAS_PeripheralOxygen)),
                                                                       100*sum(is.na(BAS_PeripheralOxygen) / n())),

      "Highest respiratory rate (breaths/minute)_Median (IQR)" = sprintf("%.0f (%.0f, %.0f)",
                                                                         median(BAS_RespRateHighest, na.rm = TRUE),
                                                                         quantile(BAS_RespRateHighest, prob = 0.25, na.rm = TRUE),
                                                                         quantile(BAS_RespRateHighest, prob = 0.75, na.rm = TRUE)),
      "Highest respiratory rate_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                                            sum(is.na(BAS_RespRateHighest)),
                                                            100*sum(is.na(BAS_RespRateHighest) / n())),

      "Highest recorded Urea (mmol/L)_Median (IQR)" = sprintf("%.0f (%.0f, %.0f)",
                                                              median(BAS_UreaResult, na.rm = TRUE),
                                                              quantile(BAS_UreaResult, prob = 0.25, na.rm = TRUE),
                                                              quantile(BAS_UreaResult, prob = 0.75, na.rm = TRUE)),
      "Highest recorded Urea (mmol/L)_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                                                  sum(is.na(BAS_UreaResult)),
                                                                  100*sum(is.na(BAS_UreaResult) / n())),

      "Highest recorded CRP (mg/L)_Median (IQR)" = sprintf("%.0f (%.0f, %.0f)",
                                                           median(BAS_CRPResult, na.rm = TRUE),
                                                           quantile(BAS_CRPResult, prob = 0.25, na.rm = TRUE),
                                                           quantile(BAS_CRPResult, prob = 0.75, na.rm = TRUE)),
      "Highest recorded CRP (mg/L)_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                                               sum(is.na(BAS_CRPResult)),
                                                               100*sum(is.na(BAS_CRPResult) / n())),
      "APTT_Median (IQR)" = sprintf("%.0f (%.0f, %.0f)", median(BAS_APTT, na.rm = TRUE),
                                    quantile(BAS_APTT, prob = 0.25, na.rm = TRUE),
                                    quantile(BAS_APTT, prob = 0.75, na.rm = TRUE)),
      "APTT_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                        sum(is.na(BAS_APTT)),
                                        100*sum(is.na(BAS_APTT) / n())),
      "INR_Mean (SD)" = sprintf("%.2f (%.2f)", mean(BAS_INR, na.rm = TRUE), sd(BAS_INR, na.rm = TRUE)),
      "INR_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                       sum(is.na(BAS_INR)),
                                       100*sum(is.na(BAS_INR) / n())),
      "Fibrinogen (g/L)_Mean (SD)" = sprintf("%.2f (%.2f)", mean(BAS_FibrinogenResult, na.rm = TRUE), sd(BAS_FibrinogenResult, na.rm = TRUE)),
      "Fibrinogen (g/L)_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                                    sum(is.na(BAS_FibrinogenResult)),
                                                    100*sum(is.na(BAS_FibrinogenResult) / n())),
      "Prothrombin time (sec)_Median (IQR)" = sprintf("%.0f (%.0f, %.0f)", median(BAS_ProthrombinTime, na.rm = TRUE),
                                                      quantile(BAS_ProthrombinTime, prob = 0.25, na.rm = TRUE),
                                                      quantile(BAS_ProthrombinTime, prob = 0.75, na.rm = TRUE)),
      "Prothrombin time (sec)_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                                          sum(is.na(BAS_ProthrombinTime)),
                                                          100*sum(is.na(BAS_ProthrombinTime) / n())),
      "Taking aspirin_Yes, n (\\%)" = sprintf("%i (%.0f)",
                                              sum(BAS_PatientTakingAspirin == "Yes", na.rm = TRUE),
                                              100 * sum(BAS_PatientTakingAspirin == "Yes", na.rm = TRUE) / n()),
      "Taking aspirin_Missing, n (\\%)" = sprintf("%i (%.0f)", sum(is.na(BAS_PatientTakingAspirin)), 100 * sum(is.na(BAS_PatientTakingAspirin)) / n())
    ) %>%
    gather(Variable, value, -!!grpvar, factor_key = TRUE) %>%
    spread(!!grpvar, value)
  colnames(tab)[-1] <- dat %>%
    count(!!grpvar) %>%
    mutate(lab = paste0(!!grpvar, "<br>(n = ", n, ")")) %>%
    pull(lab)
  return(tab)
}


#' Generates baseline prognostics summary overall
#' @param dat The dataset
#' @return A tibble giving the summary
generate_baseline_prognostics <- function(dat) {
  tab <- dat %>%
    summarise(
      "On Room Air_Yes, n (\\%)" = sprintf("%i (%.0f)",
                                           sum(BAS_OnRoomAir24hrs == "Yes", na.rm = TRUE),
                                           100 * sum(BAS_OnRoomAir24hrs == "Yes", na.rm = TRUE) / n()),
      "On Room Air_Missing, n (\\%)" = sprintf("%i (%.0f)", sum(is.na(BAS_OnRoomAir24hrs)), 100 * sum(is.na(BAS_OnRoomAir24hrs)) / n()),

      "GCS < 15_Yes, n (\\%)" = sprintf("%i (%.0f)",
                                        sum(BAS_PatientGCS == "Yes", na.rm = TRUE),
                                        100*sum(BAS_PatientGCS == "Yes", na.rm = TRUE)/n()),
      "GCS < 15_Missing, n (\\%)" = sprintf("%i (%.0f)", sum(is.na(BAS_PatientGCS)), 100*sum(is.na(BAS_PatientGCS)) / n()),

      "Peripheral oxygen saturation (SpO2)_Median (IQR)" = sprintf("%.0f (%.0f, %.0f)",
                                                                   median(BAS_PeripheralOxygen, na.rm = TRUE),
                                                                   quantile(BAS_PeripheralOxygen, na.rm = TRUE, prob = 0.25),
                                                                   quantile(BAS_PeripheralOxygen, na.rm = TRUE, prob = 0.75)),
      "Peripheral oxygen saturation (SpO2)_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                                                       sum(is.na(BAS_PeripheralOxygen)),
                                                                       100*sum(is.na(BAS_PeripheralOxygen) / n())),

      "Highest respiratory rate (breaths/minute)_Median (IQR)" = sprintf("%.0f (%.0f, %.0f)",
                                                                         median(BAS_RespRateHighest, na.rm = TRUE),
                                                                         quantile(BAS_RespRateHighest, prob = 0.25, na.rm = TRUE),
                                                                         quantile(BAS_RespRateHighest, prob = 0.75, na.rm = TRUE)),
      "Highest respiratory rate_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                                            sum(is.na(BAS_RespRateHighest)),
                                                            100*sum(is.na(BAS_RespRateHighest) / n())),

      "Highest recorded Urea (mmol/L)_Median (IQR)" = sprintf("%.0f (%.0f, %.0f)",
                                                              median(BAS_UreaResult, na.rm = TRUE),
                                                              quantile(BAS_UreaResult, prob = 0.25, na.rm = TRUE),
                                                              quantile(BAS_UreaResult, prob = 0.75, na.rm = TRUE)),
      "Highest recorded Urea (mmol/L)_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                                                  sum(is.na(BAS_UreaResult)),
                                                                  100*sum(is.na(BAS_UreaResult) / n())),

      "Highest recorded CRP (mg/L)_Median (IQR)" = sprintf("%.0f (%.0f, %.0f)",
                                                           median(BAS_CRPResult, na.rm = TRUE),
                                                           quantile(BAS_CRPResult, prob = 0.25, na.rm = TRUE),
                                                           quantile(BAS_CRPResult, prob = 0.75, na.rm = TRUE)),
      "Highest recorded CRP (mg/L)_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                                               sum(is.na(BAS_CRPResult)),
                                                               100*sum(is.na(BAS_CRPResult) / n())),
      "APTT_Median (IQR)" = sprintf("%.0f (%.0f, %.0f)", median(BAS_APTT, na.rm = TRUE),
                                    quantile(BAS_APTT, prob = 0.25, na.rm = TRUE),
                                    quantile(BAS_APTT, prob = 0.75, na.rm = TRUE)),
      "APTT_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                        sum(is.na(BAS_APTT)),
                                        100*sum(is.na(BAS_APTT) / n())),
      "INR_Mean (SD)" = sprintf("%.2f (%.2f)", mean(BAS_INR, na.rm = TRUE), sd(BAS_INR, na.rm = TRUE)),
      "INR_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                       sum(is.na(BAS_INR)),
                                       100*sum(is.na(BAS_INR) / n())),
      "Fibrinogen (g/L)_Mean (SD)" = sprintf("%.2f (%.2f)", mean(BAS_FibrinogenResult, na.rm = TRUE), sd(BAS_FibrinogenResult, na.rm = TRUE)),
      "Fibrinogen (g/L)_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                                    sum(is.na(BAS_FibrinogenResult)),
                                                    100*sum(is.na(BAS_FibrinogenResult) / n())),
      "Prothrombin time (sec)_Median (IQR)" = sprintf("%.0f (%.0f, %.0f)", median(BAS_ProthrombinTime, na.rm = TRUE),
                                                      quantile(BAS_ProthrombinTime, prob = 0.25, na.rm = TRUE),
                                                      quantile(BAS_ProthrombinTime, prob = 0.75, na.rm = TRUE)),
      "Prothrombin time (sec)_Missing, n (\\%)" = sprintf("%i (%.0f)",
                                                          sum(is.na(BAS_ProthrombinTime)),
                                                          100*sum(is.na(BAS_ProthrombinTime) / n())),
      "Taking aspirin_Yes, n (\\%)" = sprintf("%i (%.0f)",
                                              sum(BAS_PatientTakingAspirin == "Yes", na.rm = TRUE),
                                              100 * sum(BAS_PatientTakingAspirin == "Yes", na.rm = TRUE) / n()),
      "Taking aspirin_Missing, n (\\%)" = sprintf("%i (%.0f)", sum(is.na(BAS_PatientTakingAspirin)), 100 * sum(is.na(BAS_PatientTakingAspirin)) / n())

    ) %>%
    gather(Variable, value, factor_key = TRUE)
  colnames(tab)[2] <- paste0("Overall<br>(n = ", nrow(dat), ")")
  return(tab)
}


#' Generates baseline prognostics summary tables
#' @param dat The dataset
#' @param closed If TRUE, generates summary only overall, if FALSE generates by domain
#' @return A tibble giving the summary
generate_baseline_prognostics_table <- function(dat, format = "html") {
    ovrC   <- generate_baseline_prognostics(dat %>% filter(CAssignment != "C0"))
    byCgrp <- generate_baseline_prognostics_by(dat %>% filter(CAssignment != "C0"), CAssignment)
    tabC <- left_join(ovrC, byCgrp, by = "Variable")[, c(1,(2 + 1:(ncol(byCgrp) - 1)), 2)]
    if(format == "latex") {
      colnames(tabC) <- linebreak(colnames(tabC), linebreaker = "<br>", align = "c")
    }
    outC <- tabC %>%
      mutate(Variable = str_replace_all(Variable, ".*_", "")) %>%
      kable(
        format = format,
        booktabs = TRUE,
        escape = F,
        linesep = "",
        caption = "Baseline prognostic variables for participants randomised into domain C.",
        align = "lrrrrr") %>%
      kable_styling(
        bootstrap_options = "striped",
        font_size = 11,
        latex_options = "HOLD_position") %>%
      group_rows("Was the patient on room air for any of the preceding 24 hours?", 1, 2) %>%
      group_rows("Was the patient's GCS < 15?", 3, 4) %>%
      group_rows("Peripheral oxygen saturation (SpO2) on room air (Lowest)", 5, 6) %>%
      group_rows("Highest respiratory rate (breaths/minute)", 7, 8) %>%
      group_rows("Highest recorded Urea in the last 24 hours (mmol/L)", 9, 10) %>%
      group_rows("Highest recorded CRP in the last 24 hours (mg/L)", 11, 12) %>%
      group_rows("APTT\\\\textsuperscript{1}", 13, 14, escape = F) %>%
      group_rows("INR\\\\textsuperscript{1}", 15, 16, escape = F) %>%
      group_rows("Fibrinogen\\\\textsuperscript{1} (g/L)", 17, 18, escape = F) %>%
      group_rows("Prothrombin time\\\\textsuperscript{1} (sec)", 19, 20, escape = F) %>%
      group_rows("Taking aspirin", 21, 22) %>%
      add_header_above(c(" " = 1, "Anticoagulation" = ncol(byCgrp) - 1, " " = 1)) %>%
      row_spec(0, align = "c") %>%
      add_footnote("For APTT, INR, Fibrinogen, and Prothrombin only at least one required.",
                   notation = "number")
    return(outC)
}
