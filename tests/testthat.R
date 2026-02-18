library(testthat)
Sys.setenv(BRCASURV_SKIP_DATA_DOWNLOAD = "true")
library(brcasurv)

test_check("brcasurv")
