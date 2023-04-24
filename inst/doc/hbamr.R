## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE)

load("data.rda")
library("dplyr")

## ----load_data, echo = TRUE, eval = FALSE, include = TRUE---------------------
#  library("hbamr")
#  data(LC1980)
#  LC1980[LC1980 == 0 | LC1980 == 8 | LC1980 == 9] <- NA
#  self <- LC1980[, 1]
#  stimuli <- LC1980[, -1]

## ----show_data, eval = FALSE--------------------------------------------------
#  head(stimuli)

## ---- echo = FALSE------------------------------------------------------------
head_stimuli 

## ----prep_data, echo = TRUE, eval = FALSE-------------------------------------
#  dat <- prep_data(self, stimuli)

## ----fithbam1, eval = FALSE---------------------------------------------------
#  fit_hbam <- hbam(self, stimuli)

## ----fitbam, eval = FALSE, cache = FALSE, message = FALSE, warning = FALSE----
#  fit_hbam_mini <- hbam(self, stimuli, model = "HBAM_MINI", allow_miss = 0)

## ----fithbam2, cache = FALSE, eval = FALSE, message = FALSE, warning = FALSE----
#  dat <- prep_data(self, stimuli, allow_miss = 0)
#  fit_hbam <- hbam(data = dat, prep_data = FALSE)

## ----prep_data3, include = FALSE, eval = FALSE--------------------------------
#  dat <- prep_data(self, stimuli, allow_miss = 0)

## ----fitfbam, eval = FALSE----------------------------------------------------
#  fit_fbam <- fbam(self, stimuli)

## ---- include = TRUE, eval = FALSE, fig.align = "center", fig.asp = .4, fig.width = 8, out.width = "100%"----
#  plot_stimuli(fit_hbam)

## ---- echo = FALSE, fig.align = "center", fig.asp = .4, fig.width = 8, out.width = "100%"----
knitr::include_graphics("p_stim.svg")

## ---- include = TRUE, eval = FALSE, fig.align = "center", fig.asp = .4, fig.width = 8, out.width = "100%"----
#  plot_respondents(fit_hbam, n_draws = 10)

## ---- echo = FALSE, fig.align = "center", fig.asp = .4, fig.width = 8, out.width = "100%"----
knitr::include_graphics("p_resp.svg")

## ----get_plot_dat, include = TRUE, eval = FALSE, include= FALSE, eval = FALSE----
#  get_plot_data(fit_hbam)

## ----plot_eta, include = TRUE, eval = FALSE, echo = TRUE, fig.asp = .67, fig.width = 4.25, fig.align = "center", out.width = "50%"----
#  plot_over_self(fit_hbam, dat, "eta")

## ---- echo = FALSE, fig.align = "center", fig.asp = .67, fig.width = 4.25, fig.align = "center", out.width = "50%"----
knitr::include_graphics("p_eta.svg")

## ----plot_alpha, include = TRUE, eval = FALSE, echo = TRUE, fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
#  plot_over_self(list(fit_hbam, fit_hbam_mini), dat, "alpha")

## ---- echo = FALSE, fig.align = "center", fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
knitr::include_graphics("p_alpha.svg")

## ----plot_beta, include = TRUE, eval = FALSE, echo = TRUE, fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
#  plot_over_self(list(fit_hbam, fit_hbam_mini), dat, "abs_beta")

## ---- echo = FALSE, fig.align = "center", fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
knitr::include_graphics("p_abs_beta.svg")

## ----plot_lambda, include = TRUE, eval = FALSE, echo = TRUE, fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
#  plot_over_self(list(fit_hbam, fit_hbam_mini), dat, "lambda")

## ---- echo = FALSE, fig.align = "center", fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
knitr::include_graphics("p_lambda.svg")

## ----plot_chi, include = TRUE, eval = FALSE, echo = TRUE, fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
#  plot_over_self(list(fit_hbam, fit_hbam_mini), dat, "chi")

## ---- echo = FALSE, fig.align = "center", fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
knitr::include_graphics("p_chi.svg")

## ----get_est_t, include = TRUE, eval = FALSE----------------------------------
#  get_est(fit_hbam, "theta")

## ---- echo = FALSE------------------------------------------------------------
est_theta

## ----get_est_v, include = TRUE, eval = FALSE----------------------------------
#  get_est(fit_hbam, "chi")

## ---- echo = FALSE------------------------------------------------------------
est_chi

## ----traceplot_code, eval = FALSE, include = TRUE-----------------------------
#  rstan::traceplot(fit_hbam, pars = "theta")

## ---- echo = FALSE, fig.align = "center", fig.asp = .5, fig.width = 8, fig.align = "center", out.width = "100%", warning = FALSE, message = FALSE----
knitr::include_graphics("p_trace_theta.svg")

## ---- include = FALSE---------------------------------------------------------
load("elpds.rda")

## ----run_cv_new, eval = FALSE-------------------------------------------------
#  elpd_hbam <- hbam_cv(self, stimuli, model = "HBAM", allow_miss = 0)
#  elpd_hbam_2 <- hbam_cv(self, stimuli, model = "HBAM_2", allow_miss = 0)
#  elpd_hbam_NE <- hbam_cv(self, stimuli, model = "HBAM_NE", allow_miss = 0)
#  elpd_hbam_0 <- hbam_cv(self, stimuli, model = "HBAM_0", allow_miss = 0)
#  elpd_bam <- hbam_cv(self, stimuli, model = "BAM", allow_miss = 0)

## ----elpd_ill, eval = FALSE, include = TRUE-----------------------------------
#  elpds <- rbind(elpd_bam, elpd_hbam, elpd_hbam0, elpd_hbam_ne, elpd_hbam_2)
#  elpds[order(elpds$ELPD), ]

## ---- echo = FALSE, include = TRUE--------------------------------------------
elpds <- rbind(elpd_bam, elpd_hbam, elpd_hbam0, elpd_hbam_ne, elpd_hbam_2)
elpds[order(elpds$ELPD), ]

## ---- eval = FALSE------------------------------------------------------------
#  elpd_hbam_HM <- hbam_cv(self, stimuli, model = "HBAM_HM", allow_miss = 0)
#  elpd_hbam_MINI <- hbam_cv(self, stimuli, model = "HBAM_MINI", allow_miss = 0)
#  rbind(elpd_hbam_hm, elpd_hbam_mini)

## ---- echo = FALSE------------------------------------------------------------
rbind(elpd_hbam_hm, elpd_hbam_mini)

