# Biostatisztika webinár 2019-10-04 Budapest ÁOTE

Teendők:

Telepítsük a szükséges R csomagokat

``` R
pkgs <- c("bookdown", "detect", "remotes", 
  "intrval", "knitr", "mefa4", "deldir",
  "shiny", "sp", "unmarked", "Distance")
to_inst <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_inst))
  install.packages(to_inst, repos="https://cloud.r-project.org/")
remotes::install_github("psolymos/bSims")
still_missing <- setdiff(c(pkgs, "bSims"), 
  rownames(installed.packages()))
if (length(still_missing)) {
  cat("Az alábbi csomagok telepítése nem sikerült:\n",
    paste("\t-", still_missing, collapse="\n"), "\n")
} else {
  loading_problem <- character(0L)
  for (pkg in c(pkgs, "bSims")) {
    i <- try(library(pkg, character.only=TRUE), silent=TRUE)
    if (inherits(i, "try-error"))
      loading_problem <- c(loading_problem, pkg)
  }
  if (length(loading_problem)) {
    cat("Az alábbi csomagok betöltése nem sikerült:\n",
      paste("\t-", loading_problem, collapse="\n"), "\n")
  } else {
    cat("Csomagok telepítéses sikeresen befejeződött.\n")
  }
}
```

A GitHub-ról töltsük le a `qpad-book` repót (klónozás v. zip fájl): https://github.com/psolymos/qpad-book

Az előadás diái [itt](https://docs.google.com/presentation/d/1ankrezLL5KBCTkQahV9I857PmFuTtTV3GahFotPrjks/edit?usp=sharing) érhetők el.

