library(htmlwidgets)
library(htmltools)
library(devtools)
dat <- data.frame(Peanut = runif(45), Apple = runif(45)+runif(45), Banana = c(rep("Cat", 25), rep("Bat", 10), rep("Dog",10)))

# run
install("~/Documents/Projects/mead/filterWidget/")
library(filterWidget)
filterWidget("Peanut", dat)
#categorical
filterWidget("Banana", dat)
