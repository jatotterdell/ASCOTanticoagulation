# ASCOT Anticoagulation Domain Analyses

Analysis scripts for the close-out of the ASCOT anticoagulation domain.

## Project Setup

The project assumes that a `.Renviron` file exists in the project directory and defines a variable `ASCOT_DATA` which provides the top level directory of the RDS location. On my system this is:

```{r}
ASCOT_DATA = "/media/rds/PRJ-ascotsims"
```

The `.Rprofile` file then defines a few additional global constants.

Both are automatically loaded when opening `ASCOTanticoagulation.Rproj`.

## Raw Data

If desired, the raw datasets can be added to the global environment with

```{r}
source("r/read_raw_data.r")
read_raw_extracts("enrolled")
```

Otherwise, it is better to work from the derived datasets.
