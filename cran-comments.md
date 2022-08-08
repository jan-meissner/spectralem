Possible Errors:
* Depends on `RcppFaddeeva`. Which is not updated to the current version of R
  * `remotes::install_version('RcppFaddeeva', '0.2.2')` works for github workflows.
* Is it ok to include the code of `RcppFaddeeva`? (ofcourse denote it in LICENSE)

Notes:
* Usage of `plotly` causes large vignettes (3.9MB), most likely causing `sub-directories of 1Mb or more: doc 5MB`. 
  Could use `ggplot2` instead but then can't zoom to distinguish the fit from the real signal.
