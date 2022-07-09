# colab-r-libs

Compiling the following chunk on Google Colab on a new instance takes 13 minutes:
```r
if (!require("BiocManager", quietly = TRUE)) {
    install.packages(
        "BiocManager",
        quiet = TRUE
    )
}
packages <- c(
  "DESeq2",
  "apeglm",
  "EnhancedVolcano",
  "data.table",
  "tidyr"
)
BiocManager::install(
    packages,
    update = FALSE,
    quiet = TRUE,
    Ncpus = parallel::detectCores()
)
invisible(
    lapply(
      packages,
      library,
      character.only = TRUE
    )
)
```

Putting the content of the `/usr/local/lib/R/site-library/` in this repository, then cloning it on a new Colab instance allows to load the same R packages in less than 1 minute.

