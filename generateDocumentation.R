setwd("/Users/annaluo/Wu/Robust_NLME")
library(sinew)

files <- list.files(
  "R",
  pattern = "\\.R$",
  full.names = TRUE
)

for (f in files) {
  makeOxyFile(f, overwrite = TRUE)
}
