D <- "results"
q1 <- read.csv(file.path(D, "summary/q1_fit.csv"))
q2 <- read.csv(file.path(D, "summary/q2_tilt.csv"))
q3 <- read.csv(file.path(D, "summary/q3_minima.csv"))
hw <- read.csv(file.path(D, "summary/heywood.csv"))
pop <- readRDS(file.path(D, "population.rds"))
