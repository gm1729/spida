.onLoad <- function(x,y) {
  packageStartupMessage("\nspida: York Summer Program in Data Analysis\nVersion: 0.1-3 (2014 09 22)\nNote:\n'nlme' is no longer automatically loaded to allow users to load 'lme4' without conflicting with 'nlme'.\nUse '> library(nlme)' or '> library(lme4)' to load the package you want to use.")
}


# .onLoad() 