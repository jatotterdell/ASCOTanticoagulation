# Adherence analysis ------------------------------------------------------
# Christine Sommerville's adherence scripts

# Note for therapeutic dose
# Enoxaparin: PP_Dose is 1 x BAS_Weight twice per day, or 1.5 x BAS_Weight once
# per day. PP_Dose is in milligrams (mg) at the frequency in the 'times' column.
# Tinzaparin: PP_Dose is IU/kg - Note only 1 patient received Tinzaparin and
# none received Dalteparin.

# common to protocol 3 (incl C3 arm) and protocol 5.
# Matches James Totterdell's original pp_lmwh_dose table that did not account
# for creatinine levels at baseline and assumes creatinine clearance >30mL/min.
pp_lwmh_dose <- function() {
  tribble(
    ~DD_TypeLMWH, ~times, ~weight_low, ~weight_high, ~weight_group, ~CAssignment, ~PP_Dose,
    "Enoxaparin", "Once", 0, 50, "< 50", "C1", 20,
    "Enoxaparin", "Once", 0, 50, "< 50", "C2", 40,
    "Enoxaparin", "Once", 0, 50, "< 50", "C3", 20,
    "Enoxaparin", "Twice", 0, 50, "< 50", "C4", 1,
    "Enoxaparin", "Once", 0, 50, "< 50", "C4", 1.5,
    "Enoxaparin", "Once", 50, 120, "50 - 120", "C1", 40,
    "Enoxaparin", "Twice", 50, 120, "50 - 120", "C2", 40,
    "Enoxaparin", "Once", 50, 120, "50 - 120", "C2", 80,
    "Enoxaparin", "Once", 50, 120, "50 - 120", "C3", 40,
    "Enoxaparin", "Twice", 50, 120, "50 - 120", "C4", 1,
    "Enoxaparin", "Once", 50, 120, "50 - 120", "C4", 1.5,
    "Enoxaparin", "Once", 120, Inf, "> 120", "C1", 60,
    "Enoxaparin", "Twice", 120, Inf, "> 120", "C2", 60,
    "Enoxaparin", "Once", 120, Inf, "> 120", "C2", 120,
    "Enoxaparin", "Once", 120, Inf, "> 120", "C3", 60,
    "Enoxaparin", "Twice", 120, Inf, "> 120", "C4", 1,
    "Enoxaparin", "Once", 120, Inf, "> 120", "C4", 1.5,
    "Tinzaparin", "Once", NA, NA, NA, "C1", 75,
    "Tinzaparin", "Once", NA, NA, NA, "C2", 125,
    "Tinzaparin", "Once", NA, NA, NA, "C3", 75,
    "Tinzaparin", "Once", NA, NA, NA, "C4", 175
  )
}

# Protocol 5 Dosing table is not applicable to C3 (standard dose + aspirin)
# Different means differs from creatinine over >30 dosing
pp_lwmh_dose_creat_under <- function() {
  tribble(
    ~DD_TypeLMWH, ~times, ~weight_low, ~weight_high, ~weight_group, ~CAssignment, ~PP_Dose,
    "Enoxaparin", "Once", 0, 50, "< 50", "C1", 20,
    "Enoxaparin", "Once", 0, 50, "< 50", "C2", 0.5, # different
    "Enoxaparin", "Once", 0, 50, "< 50", "C4", 1, # different
    "Enoxaparin", "Once", 50, 120, "50 - 120", "C1", 20,
    "Enoxaparin", "Once", 0, 50, "50 - 120", "C2", 0.5, # different
    "Enoxaparin", "Once", 50, 120, "50 - 120", "C4", 1,
    "Enoxaparin", "Once", 120, Inf, "> 120", "C1", 40,
    "Enoxaparin", "Once", 120, Inf, "> 120", "C2", 0.5,
    "Enoxaparin", "Once", 120, Inf, "> 120", "C4", 1,
    "Tinzaparin", "Once", NA, NA, NA, "C1", 75,
    "Tinzaparin", "Once", NA, NA, NA, "C2", 125,
    "Tinzaparin", "Once", NA, NA, NA, "C4", 175
  )
}

# Prepare LMWH dosing table -----

generate_lmwh_dosing_table <- function(ddat) {
  lwmhdat <- ddat %>%
    filter(grepl("[1-4]", CAssignment)) %>% # omits anyone not in Domain C
    select(
      StudyPatientID,
      CAssignment,
      BAS_Weight,
      DD_StudyDay,
      DD_Date,
      DIS_DateOfDischarge,
      DD_LMWHAdministered,
      DD_TypeLMWH,
      DD_DoseLMWH,
      DD_LMWHAdministeredToday,
      DD_AspirinAdministered
    ) %>%
    mutate(DD_DoseLMWH = as.numeric(DD_DoseLMWH)) %>%
    arrange(StudyPatientID, DD_StudyDay) %>%
    # May not dose on day of discharge
    filter(DD_Date < DIS_DateOfDischarge) %>%
    # Disregard dosing on day 1
    filter(DD_StudyDay > 1) %>%
    # The one occasion reported in IU rather than IU/kg, convert
    mutate(DD_DoseLMWH = if_else(DD_TypeLMWH == "Tinzaparin", DD_DoseLMWH / BAS_Weight, DD_DoseLMWH))
  n_by_CAssignment <- lwmhdat %>%
    group_by(CAssignment, DD_LMWHAdministered, DD_TypeLMWH, DD_LMWHAdministeredToday, DD_DoseLMWH) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(CAssignment) %>%
    mutate(value = sprintf("%s (%3.1f)", formatC(n, big.mark = ","), 100 * n / sum(n))) %>%
    ungroup() %>%
    select(-n) %>%
    spread(CAssignment, value, fill = "0 (0.0)")
  n_total <- lwmhdat %>%
    group_by(DD_LMWHAdministered = "Total days", CAssignment) %>%
    summarise(value = formatC(n(), big.mark = ","), .groups = "drop")
  tabdat <- n_by_CAssignment %>%
    rename(
      Administered = DD_LMWHAdministered,
      Type = DD_TypeLMWH,
      Occasions = DD_LMWHAdministeredToday,
      Dose = DD_DoseLMWH
    ) %>%
    mutate(Administered = factor(Administered, levels = c("Yes", "No"))) %>%
    arrange(Administered)
  colnames(tabdat)[-(1:4)] <- n_total %>%
    mutate(lab = linebreak(paste0(CAssignment, "\n(Days = ", value, ")"), align = "c")) %>%
    pull(lab)
  tab <- kable(tabdat,
               booktabs = T, linesep = "",
               caption = "Daily dosing of low molecular weight heparin participants randomised to domain C.",
               align = "lrrrrrrrr", escape = F, longtable = T
  ) %>%
    kable_styling(font_size = 9, latex_options = "HOLD_position") %>%
    collapse_rows(2:3) %>%
    add_header_above(c(" " = 4, "Anticoagulation" = 4))
  return(tab)
}

# Check dosings are per-protocol for patients with creatinine clearance over
# threshold 30ml/min

check_lmwh_dose <- function(ddat) {
  ddat %>%
    # filter(grepl("[1-4]", CAssignment))
    mutate(
      expected_dose = case_when(
        CAssignment %in% c("C1", "C3") & DD_TypeLMWH == "Enoxaparin" & BAS_Weight < 50 & DD_LMWHAdministeredToday == "Once" ~ 20,
        CAssignment %in% c("C1", "C3") & DD_TypeLMWH == "Enoxaparin" & (BAS_Weight >= 50 & BAS_Weight < 120) & DD_LMWHAdministeredToday == "Once" ~ 40,
        CAssignment %in% c("C1", "C3") & DD_TypeLMWH == "Enoxaparin" & BAS_Weight > 120 & DD_LMWHAdministeredToday == "Once" ~ 60,
        CAssignment %in% c("C1", "C3") & DD_TypeLMWH == "Tinzaparin" & DD_LMWHAdministeredToday == "Once" ~ 75 * BAS_Weight,

        # Invalid dosing as no "Twice per day" option specified in protocol
        CAssignment %in% c("C1", "C3") & DD_TypeLMWH == "Enoxaparin" & DD_LMWHAdministeredToday == "Twice" ~ -9999,

        # C2 per-protocol dosing check
        CAssignment == "C2" & DD_TypeLMWH == "Enoxaparin" & BAS_Weight < 50 & DD_LMWHAdministeredToday == "Once" ~ 40,
        CAssignment == "C2" & DD_TypeLMWH == "Enoxaparin" & (BAS_Weight >= 50 & BAS_Weight < 120) & DD_LMWHAdministeredToday == "Once" ~ 80,
        CAssignment == "C2" & DD_TypeLMWH == "Enoxaparin" & (BAS_Weight >= 50 & BAS_Weight < 120) & DD_LMWHAdministeredToday == "Twice" ~ 40,
        CAssignment == "C2" & DD_TypeLMWH == "Enoxaparin" & BAS_Weight > 120 & DD_LMWHAdministeredToday == "Once" ~ 120,
        CAssignment == "C2" & DD_TypeLMWH == "Enoxaparin" & BAS_Weight > 120 & DD_LMWHAdministeredToday == "Twice" ~ 60,
        CAssignment == "C2" & DD_TypeLMWH == "Tinzaparin" & DD_LMWHAdministeredToday == "Once" ~ 125 * BAS_Weight,

        # Invalid dosing as no "Twice per day" option specified in protocol
        CAssignment == "C2" & DD_TypeLMWH == "Enoxaparin" & BAS_Weight < 50 & DD_LMWHAdministeredToday == "Twice" ~ -9999,

        # C4 per-protocol dosing check
        CAssignment == "C4" & DD_TypeLMWH == "Enoxaparin" & DD_LMWHAdministeredToday == "Once" ~ 1.5 * BAS_Weight,
        CAssignment == "C4" & DD_TypeLMWH == "Enoxaparin" & DD_LMWHAdministeredToday == "Twice" ~ BAS_Weight,
        CAssignment == "C4" & DD_TypeLMWH == "Tinzaparin" & DD_LMWHAdministeredToday == "Once" ~ 175 * BAS_Weight,

        is.na(CAssignment) | is.na(DD_TypeLMWH) | is.na(DD_LMWHAdministeredToday) | is.na(DD_DoseLMWH) | is.na(BAS_Weight) ~ NA_real_,
        TRUE ~ NA_real_
      ),

      # actual_dose = if_else((CAssignment == "C0" & is.na(DD_LMWHAdministeredToday)), 0, DD_DoseLMWH)),
      dose_given = DD_LMWHAdministered,
      dose_equal_expected = expected_dose == DD_DoseLMWH,
      dose_close_expected = abs(expected_dose - DD_DoseLMWH) < 10
    )
}
