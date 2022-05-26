ASCOT_DATA <- Sys.getenv("ASCOT_DATA")
ASCOT_DATA_RAW <- file.path(ASCOT_DATA, "raw")
ANTICOAG_DATA <- file.path(ASCOT_DATA, "anticoagulation")

if (!dir.exists("outputs")) {
  dir.create(file.path("outputs", "tables"), recursive = TRUE)
  dir.create(file.path("outputs", "figures"), recursive = TRUE)
  dir.create(file.path("outputs", "models"), recursive = TRUE)
}
