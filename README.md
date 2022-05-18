# ASCOT Anticoagulation Domain Analyses

Analysis scripts for the close-out of the ASCOT anticoagulation domain.

## Project Setup

The project assumes that a `.Renviron` file exists in the project directory and defines a variable `ASCOT_DATA` which provides the top level directory of the RDS location. On my system this is:

```{r}
ASCOT_DATA = "/media/rds/PRJ-ascotsims"
```

The `.Rprofile` file then defines a few additional global constants.

Both are automatically loaded when opening `ASCOTanticoagulation.Rproj`.

### Quarto

Interactive notebooks are written using [Quarto](https://quarto.org/docs/get-started/).
To render/preview the notebooks, run them from the root project directory via

```
quarto preview <filename>.qmd --execute-dir .
```

### Directories

Planned structure

```
.
├── analyses
├── docs
├── eda
├── outputs
│   ├── figures
│   ├── models
│   └── tables
├── r
├── report
└── stan
```


## Raw Data

If desired, the raw datasets can be added to the global environment with, for example,

```{r}
source("r/read_raw_data.r")
read_raw_extract("enrolled")
read_all_raw_extracts()
```

Otherwise, it is better to work from the derived datasets.
