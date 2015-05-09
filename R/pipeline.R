## Pipeline Convenience functions ####

## Sort a data frame ####
# GM 2014 08 09

sortdf <- function(data, form) {
  xx <- as.list(model.frame(form, data, na.action=NULL))
  ord <- do.call(order, xx)
  data[ ord,]
}

## Assign in a pipeline

assn <- function(value,x, pos = 0, ...) assign(x=x,value=value, pos = pos, ...)
