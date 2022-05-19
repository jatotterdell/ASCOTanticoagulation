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
To render the notebooks, either use `Build` in RStudio or use the CLI by running

```
quarto render .
```

to render all documents, or

```
quarto render analyses/<filename>.qmd
quarto preview analyses/<filename>.qmd
```

to render/preview a specific file only.

Alternatively, use 

Alternatively can render 

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
