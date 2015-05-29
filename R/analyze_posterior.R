pp <- read.table(file = "../posterior_table", header = TRUE)

library(ggplot2)
library(tidyr)

ggplot(pp, aes(x = Pos.wo.LLR.omega,  y = Pos.wo.LLR.E.t)) +
  geom_point()

library(reshape2)

tau_col_idx <- grep("Pos.LLR.Tau.+", names(pp), value = TRUE)
pp_tau <- melt(pp[, tau_col_idx], measure.vars = tau_col_idx)
ggplot(pp_tau, aes(x = value)) + geom_density() + facet_wrap( ~ variable)

psi_col_idx <- grep("Psi", names(pp), value = TRUE)
pp_psi <- melt(pp[, psi_col_idx], measure.vars = psi_col_idx)
ggplot(pp_psi, aes(x = value)) + geom_density() + facet_wrap( ~ variable)
