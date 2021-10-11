library(rstan)
library(xlsx)
rstan_options(auto_write = TRUE)

#The following options are recommended when running a multicore CPU with excess RAM.  Comment out if this is not the case.
options(mc.cores = parallel::detectCores())
options(boot.parallel = "multicore")
#end comment

data.file.name <- "FramingEffectsDataSet.xlsx"
codes.file.name <- "FramingEffectsCodes.xlsx"
output.file.name <- "Outputs.xlsx"
model.file.name <- "ProspectTheoryAnalysis.stan"

chains <- 16
iter <- 6000
warmup <- 1000
thin <- 1
adapt.delta <- .95
max.treedepth <- 15

file.data <- read.xlsx(data.file.name, 1) 
codes.data <- read.xlsx(codes.file.name, 1)
ids <- file.data[1]
file.data[is.na(file.data)] <- 0 ## stan does not allow na, so rewrite as 0
file.data <- file.data[,3:ncol(file.data)] ##cut off first columns - they're just ids
stan.data <- list(N = nrow(file.data),
                  choices = file.data,
                  codes = codes.data)
control <- list(adapt_delta = adapt.delta,
                max_treedepth = max.treedepth)

fit <- stan(file = model.file.name, 
            data = stan.data,
            chains = chains, 
            iter = iter, 
            warmup = warmup, 
            control = control)

print(fit)
#traceplot(fit, pars = c("alpha_sd", "alphas[1]", "c_sd", "cs[1]", "delta_sd", "deltas[1]", "lambda_sd", "lambdas[1]"))
write.xlsx(summary(fit)$summary, "Outputs2.xlsx")