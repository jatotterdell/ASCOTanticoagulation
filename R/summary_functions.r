library(tidyverse)
library(kableExtra)

# Outcome summary functions - model outputs ----

#' odds_ratio_summary_table
#'
#' @param OR A vector of rvars which represent posterior odds ratio draws
#' @param format Either html or latex
#' @param fn Optional file name if saving the latex table.
#'  Note that this uses `save_tex_table` which prepends `outputs/tables`
#'  to the file path and appends .tex  e.g. should just use fn = "outcomes/primary/filename"
#'  to save a file to outputs/tables/outcomes/primary/filename.tex
#' @return A kable object
odds_ratio_summary_table <- function(OR, format = "html", fn = NULL) {
  out <- tibble(
    Parameter = names(OR),
    Median = median(OR),
    `95% CrI` = apply(
      quantile(OR, c(0.025, 0.975)), 2,
      function(x) sprintf("(%.2f, %.2f)", x[1], x[2])),
    `Mean (SD)` = sprintf("%.2f (%.2f)", E(OR), sd(OR)),
    `Pr(OR < 1)` = Pr(OR < 1),
  ) %>%
    kable(
      format = format,
      digits = 2,
      align = "lrrrr",
      linesep = "",
      booktabs = TRUE) %>%
    kable_styling(
      font_size = 9,
      bootstrap_options = "striped",
      latex_options = "HOLD_position")
  if (!is.null(fn) & format == "latex") {
    save_tex_table(out, fn)
  } else {
    return(out)
  }
}

#' odds_ratio_summary_table_rev
#'
#' As for odds_ratio_summary_table but usin Pr(OR > 1) instead of Pr(OR < 1)
odds_ratio_summary_table_rev <- function(OR, format = "html", fn = NULL) {
  out <- tibble(
    Parameter = names(OR),
    Median = median(OR),
    `95% CrI` = apply(
      quantile(OR, c(0.025, 0.975)), 2,
      function(x) sprintf("(%.2f, %.2f)", x[1], x[2])),
    `Mean (SD)` = sprintf("%.2f (%.2f)", E(OR), sd(OR)),
    `Pr(OR > 1)` = Pr(OR > 1),
  ) %>%
    kable(
      format = format,
      digits = 2,
      align = "lrrrr",
      linesep = "",
      booktabs = TRUE) %>%
    kable_styling(
      font_size = 9,
      bootstrap_options = "striped",
      latex_options = "HOLD_position")
  if (!is.null(fn) & format == "latex") {
    save_tex_table(out, fn)
  } else {
    return(out)
  }
}


plot_or_densities <- function(rvs) {
  tibble(Contrast = fct_inorder(names(rvs)), RV = rvs) %>%
    ggplot(., aes(y = Contrast, xdist = RV)) +
    stat_halfeye(
      aes(fill =
            stat(cut_cdf_qi(
              cdf,
              .width = c(.5, .8, .95, 0.99),
              labels = scales::percent_format()))),
      adjust = 1, n = 1001, .width = c(0.5, 0.8, 0.95)
    ) +
    scale_fill_brewer(
      palette = "Reds",
      direction = -1,
      na.translate = FALSE) +
    labs(
      x = "Odds ratio contrast",
      fill = "Interval"
    ) +
    scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10)) +
    geom_vline(xintercept = 1)
}


plot_epoch_terms <- function(rvs_epoch) {
  orsdat <- tibble(
    Group = "Epoch",
    Parameter = fct_inorder(names(rvs_epoch)),
    posterior = rvs_epoch
  )
  p_epoch <- ggplot(orsdat, aes(xdist = posterior, y = Parameter)) +
    stat_pointinterval(.width = c(0.75, 0.95), fatten_point = 1.5) +
    geom_vline(xintercept = 1, linetype = 2) +
    scale_x_log10("Odds ratio (log scale)") +
    labs(y = "Epoch")
  return(p_epoch)
}


plot_site_terms <- function(rvs_site, region) {
  orsdat <- tibble(
    Group = "Site",
    Country = region,
    Parameter = fct_inorder(names(rvs_site)),
    posterior = rvs_site
  )
  p_site <- ggplot(orsdat, aes(xdist = posterior, y = Parameter)) +
    facet_grid(Country ~ ., scales = "free_y", space = "free_y") +
    stat_pointinterval(.width = c(0.75, 0.95), fatten_point = 1.5) +
    geom_vline(xintercept = 1, linetype = 2) +
    scale_x_log10("Odds ratio (log scale)") +
    labs(y = "Site") +
    theme(panel.border = element_rect(fill = NA))
}


plot_epoch_site_terms <- function(rvs_epoch, rvs_site, region) {
  p_epoch <- plot_epoch_terms(rvs_epoch)
  p_site <- plot_site_terms(rvs_site, region)
  p <- p_epoch | p_site
  p
}

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

      `Country_India, n (\\%)` = sprintf("%i (%.0f)", sum(Country == "IN", na.rm = TRUE), 100 * sum(Country == "IN", na.rm = TRUE) / n()),
      `Country_Australia, n (\\%)` = sprintf("%i (%.0f)", sum(Country == "AU", na.rm = TRUE), 100 * sum(Country == "AU", na.rm = TRUE) / n()),
      `Country_Nepal, n (\\%)` = sprintf("%i (%.0f)", sum(Country == "NP", na.rm = TRUE), 100 * sum(Country == "NP", na.rm = TRUE) / n()),
      `Country_New Zealand, n (\\%)` = sprintf("%i (%.0f)", sum(Country == "NZ", na.rm = TRUE), 100 * sum(Country == "NZ", na.rm = TRUE) / n()),

      `Sex_Male, n (\\%)` = sprintf("%i (%.0f)", sum(Sex == "Male", na.rm = TRUE), 100 * sum(Sex == "Male", na.rm = TRUE) / n()),
      `Sex_Female, n (\\%)` = sprintf("%i (%.0f)", sum(Sex == "Female", na.rm = TRUE), 100 * sum(Sex == "Female", na.rm = TRUE) / n()),

      `Weight_Median, (IQR)` = sprintf("%.0f (%.0f, %.0f)", median(BAS_Weight, na.rm = T), quantile(BAS_Weight, 0.25, na.rm = T), quantile(BAS_Weight, 0.75, na.rm = T)),
      `Weight_Missing, n (\\%)` = sprintf("%i (%.0f)", sum(is.na(BAS_Weight)), 100 * sum(is.na(BAS_Weight)) / n()),

      `Vaccinated_Yes, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_PatientVaccinated == "Yes", na.rm = TRUE), 100 * sum(BAS_PatientVaccinated == "Yes", na.rm = TRUE) / n()),
      `Vaccinated_Missing, n (\\%)` = sprintf("%i (%.0f)", sum(is.na(BAS_PatientVaccinated)), 100 * sum(is.na(BAS_PatientVaccinated)) / n()),

      `Ethnicity_Indian, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityIndian == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityIndian == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_European, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityEuropean == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityEuropean == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Asian, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityAsian == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityAsian == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Pacific Islander, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityPacificIslander == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityPacificIslander == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Middle Eastern, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityMiddleEastern == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityMiddleEastern == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Maori, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityMaori == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityMaori == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_African, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityAfrican == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityAfrican == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Aboriginal, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityAboriginal == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityAboriginal == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Latin American, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityLatinAmerican == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityLatinAmerican == "Yes", na.rm = TRUE) / n()),
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
    dplyr::count(!!grpvar) %>%
    mutate(lab = paste0(!!grpvar, "<br><br>(n = ", n, ")")) %>%
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

      `Country_India, n (\\%)` = sprintf("%i (%.0f)", sum(Country == "IN", na.rm = TRUE), 100 * sum(Country == "IN", na.rm = TRUE) / n()),
      `Country_Australia, n (\\%)` = sprintf("%i (%.0f)", sum(Country == "AU", na.rm = TRUE), 100 * sum(Country == "AU", na.rm = TRUE) / n()),
      `Country_Nepal, n (\\%)` = sprintf("%i (%.0f)", sum(Country == "NP", na.rm = TRUE), 100 * sum(Country == "NP", na.rm = TRUE) / n()),
      `Country_New Zealand, n (\\%)` = sprintf("%i (%.0f)", sum(Country == "NZ", na.rm = TRUE), 100 * sum(Country == "NZ", na.rm = TRUE) / n()),

      `Sex_Male, n (\\%)` = sprintf("%i (%.0f)", sum(Sex == "Male", na.rm = TRUE), 100 * sum(Sex == "Male", na.rm = TRUE) / n()),
      `Sex_Female, n (\\%)` = sprintf("%i (%.0f)", sum(Sex == "Female", na.rm = TRUE), 100 * sum(Sex == "Female", na.rm = TRUE) / n()),

      `Weight_Median, (IQR)` = sprintf("%.0f (%.0f, %.0f)", median(BAS_Weight, na.rm = T), quantile(BAS_Weight, 0.25, na.rm = T), quantile(BAS_Weight, 0.75, na.rm = T)),
      `Weight_Missing, n (\\%)` = sprintf("%i (%.0f)", sum(is.na(BAS_Weight)), 100 * sum(is.na(BAS_Weight)) / n()),

      `Vaccinated_Yes, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_PatientVaccinated == "Yes", na.rm = TRUE), 100 * sum(BAS_PatientVaccinated == "Yes", na.rm = TRUE) / n()),
      `Vaccinated_Missing, n (\\%)` = sprintf("%i (%.0f)", sum(is.na(BAS_PatientVaccinated)), 100 * sum(is.na(BAS_PatientVaccinated)) / n()),

      `Ethnicity_Indian, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityIndian == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityIndian == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_European, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityEuropean == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityEuropean == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Asian, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityAsian == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityAsian == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Pacific Islander, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityPacificIslander == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityPacificIslander == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Middle Eastern, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityMiddleEastern == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityMiddleEastern == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Maori, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityMaori == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityMaori == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_African, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityAfrican == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityAfrican == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Aboriginal, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityAboriginal == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityAboriginal == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Latin American, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityLatinAmerican == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityLatinAmerican == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Other, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityOther == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityOther == "Yes", na.rm = TRUE) / n()),
      `Ethnicity_Unknown, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_EthnicityUnknown == "Yes", na.rm = TRUE), 100 * sum(BAS_EthnicityUnknown == "Yes", na.rm = TRUE) / n()),

      `Smoking_Current, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_Smoking == "Current", na.rm = TRUE), 100 * sum(BAS_Smoking == "Current", na.rm = TRUE) / n()),
      `Smoking_Former, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_Smoking == "Former", na.rm = TRUE), 100 * sum(BAS_Smoking == "Former", na.rm = TRUE) / n()),
      `Smoking_Never, n (\\%)` = sprintf("%i (%.0f)", sum(BAS_Smoking == "Never", na.rm = TRUE), 100 * sum(BAS_Smoking == "Never", na.rm = TRUE) / n()),
      `Smoking_Missing, n (\\%)` = sprintf("%i (%.0f)", sum(is.na(BAS_Smoking)), 100 * sum(is.na(BAS_Smoking)) / n())

    ) %>%
    gather(Variable, Overall, factor_key = TRUE)
  colnames(overall)[2] <- paste0("Overall<br><br>(n = ", nrow(dat), ")")
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
    fsize <- 12
    if(format == "latex") {
      fsize <- 9
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
        latex_options = "HOLD_position",
        font_size = fsize) %>%
      group_rows("Country", 2, 5) %>%
      group_rows("Sex", 6, 7) %>%
      group_rows("Weight (kg)", 8, 9) %>%
      group_rows("Ethnicity", 12, 22) %>%
      group_rows("Smoking", 23, 26) %>%
      add_header_above(c(" " = 1, "Anticoagulation" = ncol(byCgrp) - 1, " " = 1))  %>%
      row_spec(0, align = "c") %>%
      footnote(number = "Site LUD did not have ethics approval for collection of vaccination status.",
               fixed_small_size = TRUE)
    if(format == "latex") {
      outC <- outC %>% group_rows(., "Vaccinated\\\\textsuperscript{1}", 10, 11, escape = FALSE)
    } else {
      outC <- outC %>% group_rows("Vaccinated<sup>1</sup>", 10, 11, escape = FALSE)
    }
    return(outC)
}

# Baseline summary - co-morbidities  ----

#' Generates baseline co-morbidity summary across all participants
#' @param dat Dataset with baseline variables
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
    mutate(value = sprintf("%i (%3.0f)", n, 100 * p)) %>%
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
    summarise(name = "Missing, n (\\%)", value = sprintf("%i (%3.0f)", sum(BAS_rec == 0), 100 * sum(BAS_rec == 0) / n())) %>%
    spread(!!grpvar, value)
  tab <- basdat %>%
    select(-BAS_rec) %>%
    summarise_all(., list(`.n` = ~ sum(.x == "Yes", na.rm = T),
                          `.p` = ~ sum(.x == "Yes", na.rm = T) / n())) %>%
    pivot_longer(-!!grpvar, names_to = c("name", "measure"), names_sep = "_\\.") %>%
    pivot_wider(id_cols = !!grpvar:name, names_from = measure, values_from = value) %>%
    mutate(value = sprintf("%i (%3.0f)", n, 100 * p)) %>%
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
  fsize <- 12
  if(format == "latex") {
    fsize <- 9
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
      font_size = fsize,
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
      "Taking aspirin_Missing, n (\\%)" = sprintf("%i (%.0f)", sum(is.na(BAS_PatientTakingAspirin)), 100 * sum(is.na(BAS_PatientTakingAspirin)) / n()),
      "Time from onset of symptoms to hospitalisation_Median (IQR)" = sprintf(
        "%.0f (%.0f, %.0f)",
        median(as.numeric(EL_AdmittedToHospital - EL_FirstSymptoms)),
        quantile(as.numeric(EL_AdmittedToHospital - EL_FirstSymptoms), 0.25),
        quantile(as.numeric(EL_AdmittedToHospital - EL_FirstSymptoms), 0.75)),
      "Time from hospitalisation to randomisation_Median (IQR)" = sprintf(
        "%.0f (%.0f, %.0f)",
        median(as.numeric(RandDate - EL_AdmittedToHospital)),
        quantile(as.numeric(RandDate - EL_AdmittedToHospital), 0.25),
        quantile(as.numeric(RandDate - EL_AdmittedToHospital), 0.75))
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
      "Taking aspirin_Missing, n (\\%)" = sprintf("%i (%.0f)", sum(is.na(BAS_PatientTakingAspirin)), 100 * sum(is.na(BAS_PatientTakingAspirin)) / n()),
      "Time from onset of symptoms to hospitalisation_Median (IQR)" = sprintf(
        "%.0f (%.0f, %.0f)",
        median(as.numeric(EL_AdmittedToHospital - EL_FirstSymptoms)),
        quantile(as.numeric(EL_AdmittedToHospital - EL_FirstSymptoms), 0.25),
        quantile(as.numeric(EL_AdmittedToHospital - EL_FirstSymptoms), 0.75)),
      "Time from hospitalisation to randomisation_Median (IQR)" = sprintf(
        "%.0f (%.0f, %.0f)",
        median(as.numeric(RandDate - EL_AdmittedToHospital)),
        quantile(as.numeric(RandDate - EL_AdmittedToHospital), 0.25),
        quantile(as.numeric(RandDate - EL_AdmittedToHospital), 0.75))

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
    fsize <- 12
    if(format == "latex") {
      fsize <- 9
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
        font_size = fsize,
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
      group_rows("Time from onset of symptoms to hospitalisation", 23, 23) %>%
      group_rows("Time from hospitalisation to randomisation", 24, 24) %>%
      add_header_above(c(" " = 1, "Anticoagulation" = ncol(byCgrp) - 1, " " = 1)) %>%
      row_spec(0, align = "c") %>%
      footnote(number = "For APTT, INR, Fibrinogen, and Prothrombin only at least one required.")
    return(outC)
}

# Discharge summary - other drugs ----

#' Generates overall summary of drugs used during hospital stay,
#' as recorded on discharge
#'
#' @param dat The dataset
#' @return A tibble giving the summary
generate_discharge_drugs <- function(dat) {
  tab <- dat %>%
    filter(DIS_rec == 1) %>%
    select(
      DIS_ReceivedAntibacterialDrugs,
      DIS_NoAntiviral = DIS_ReceivedNone,
      DIS_CamostatReceived,
      DIS_FavipiravirReceived,
      DIS_DoxycyclineReceived,
      DIS_IvermectinReceived,
      DIS_RemdesivirReceived,
      DIS_OtherAntiviral = DIS_ReceivedOther,
      DIS_NoImmunomodulatory = Dis_ImmunoNone,
      DIS_ImmunoAnakinra,
      DIS_ImmunoCorticosteroids,
      DIS_ImmunoSarilumab,
      DIS_ImmunoAzithromycin,
      DIS_ImmunoTocilizumab,
      DIS_ImmunoBaricitinib,
      DIS_ImmunoRuxolitinib,
      DIS_ImmunoTofacitinib,
      DIS_ImmunoZinc,
      DIS_OtherImmunomodulatory = DIS_IummunoOther
    ) %>%
    summarise_all(., list(Variable = ~ sprintf("%i (%.0f)", sum(.x == "Yes", na.rm = TRUE), 100 * sum(.x == "Yes", na.rm = TRUE) / n()))) %>%
    pivot_longer(everything()) %>%
    mutate(name = str_replace(name, "DIS_Received", ""),
           name = str_replace(name, "DIS_Immuno", ""),
           name = str_replace(name, "DIS_Iummuno", ""),
           name = str_replace(name, "Dis_Immuno", ""),
           name = str_replace(name, "DIS_", ""),
           name = str_replace(name, "Received_Variable", ""),
           name = str_replace(name, "_Variable", ""),
           name = fct_inorder(gsub("([[:upper:]]*)([[:upper:]][[:lower:]]+)", "\\1 \\2", name)),
           name = str_to_sentence(trimws(name)),
           name = paste0(name, ", n (\\%)"))
  colnames(tab) <- c(linebreak("Drug received"), linebreak(paste0("Overall\n(n = ", nrow(dat %>% filter(DIS_rec == 1)), ")"), align = "c"))
  return(tab)
}


#' Generates summary of drugs used during hospital stay,
#' as recorded on discharge, by grouping variable
#'
#' @param dat The dataset
#' @param grpvar The grouping variable
#' @return A tibble giving the summary
generate_discharge_drugs_by <- function(dat, grpvar = NULL) {
  grpvar <- enquo(grpvar)
  tab <- dat %>%
    filter(DIS_rec == 1) %>%
    group_by(!!grpvar) %>%
    select(
      !!grpvar,
      DIS_ReceivedAntibacterialDrugs,
      DIS_NoAntiviral = DIS_ReceivedNone,
      DIS_CamostatReceived,
      DIS_FavipiravirReceived,
      DIS_DoxycyclineReceived,
      DIS_IvermectinReceived,
      DIS_RemdesivirReceived,
      DIS_OtherAntiviral = DIS_ReceivedOther,
      DIS_NoImmunomodulatory = Dis_ImmunoNone,
      DIS_ImmunoAnakinra,
      DIS_ImmunoCorticosteroids,
      DIS_ImmunoSarilumab,
      DIS_ImmunoAzithromycin,
      DIS_ImmunoTocilizumab,
      DIS_ImmunoBaricitinib,
      DIS_ImmunoRuxolitinib,
      DIS_ImmunoTofacitinib,
      DIS_ImmunoZinc,
      DIS_OtherImmunomodulatory = DIS_IummunoOther
    ) %>%
    summarise_all(., list(Variable = ~ sprintf("%i (%.0f)", sum(.x == "Yes", na.rm = TRUE), 100 * sum(.x == "Yes", na.rm = TRUE) / n()))) %>%
    gather(name, value, -!!grpvar, factor_key = TRUE) %>%
    spread(!!grpvar, value) %>%
    mutate(name = str_replace(name, "DIS_Received", ""),
           name = str_replace(name, "DIS_Immuno", ""),
           name = str_replace(name, "DIS_Iummuno", ""),
           name = str_replace(name, "Dis_Immuno", ""),
           name = str_replace(name, "DIS_", ""),
           name = str_replace(name, "Received_Variable", ""),
           name = str_replace(name, "_Variable", ""),
           name = fct_inorder(gsub("([[:upper:]]*)([[:upper:]][[:lower:]]+)", "\\1 \\2", name)),
           name = str_to_sentence(trimws(name)),
           name = paste0(name, ", n (\\%)"))
  colnames(tab) <- c(linebreak("Drug received"), dat  %>% filter(DIS_rec == 1) %>% count(!!grpvar) %>% mutate(lab = linebreak(paste0(!!grpvar, "\n(n = ", n, ")"), align = "c")) %>% pull(lab))
  return(tab)
}


#' Drugs used during hospital stay
#'
#' Generates summary table of drugs used during hospital stay,
#' as recorded on discharge, by domain or overall.
#'
#' @param dat The dataset
#' @param closed If TRUE, generate overall table only, if FALSE one table per domain
#' @return A tibble giving the summary
generate_discharge_drugs_table <- function(dat, format = "html") {
  ovrC   <- generate_discharge_drugs(dat %>% filter(CAssignment != "C0"))
  bygrpC <- generate_discharge_drugs_by(dat %>% filter(CAssignment != "C0"), CAssignment)
  tabC <- left_join(ovrC, bygrpC, by = "Drug received")[, c(1,(2 + 1:(ncol(bygrpC) - 1)), 2)]
  fsize <- 12
  if(format == "latex") {
    fsize <- 9
    colnames(tabC) <- linebreak(colnames(tabC), linebreaker = "<br>", align = "c")
  }
  outC <- kable(
    tabC,
    format = format,
    booktabs = T,
    caption = "Drugs received during hospital stay, domain C.",
    align = "lrrrrrrrr",
    escape = F) %>%
    kable_styling(
      bootstrap_options = "striped",
      font_size = fsize,
      latex_options = "HOLD_position") %>%
    group_rows("Antivirals", 2, 8) %>%
    group_rows("Immunomodulatory", 9, 19) %>%
    add_header_above(c(" " = 1, "Anticoagulation" = ncol(bygrpC) - 1, " " = 1))
  return(outC)
}
