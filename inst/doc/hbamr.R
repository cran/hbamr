## ----include = FALSE----------------------------------------------------------
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

## ----echo = FALSE-------------------------------------------------------------
head_stimuli 

## ----prep_data, echo = TRUE, eval = FALSE-------------------------------------
#  dat <- prep_data(self, stimuli)

## ----eval = FALSE-------------------------------------------------------------
#  show_code("HBAM_MULTI")

## ----fithbam1, eval = FALSE---------------------------------------------------
#  fit_hbam <- hbam(self, stimuli)

## ----fitbam, eval = FALSE, cache = FALSE, message = FALSE, warning = FALSE----
#  fit_hbam_mini <- hbam(self, stimuli, model = "HBAM_MINI", allow_miss = 0)

## ----fithbam2, cache = FALSE, eval = FALSE, message = FALSE, warning = FALSE----
#  dat <- prep_data(self, stimuli, allow_miss = 0)
#  fit_hbam <- hbam(data = dat)

## ----prep_data3, include = FALSE, eval = FALSE--------------------------------
#  dat <- prep_data(self, stimuli, allow_miss = 0)

## ----eval = FALSE-------------------------------------------------------------
#  fit_fbam_multi <- hbam(self, stimuli, model = "FBAM_MULTI", group_id = self,
#                         sigma_alpha = .8, sigma_mu_alpha = .5,
#                         sigma_beta = .4, sigma_mu_beta = .25)

## ----fitfbam, eval = FALSE----------------------------------------------------
#  fit_fbam <- fbam(self, stimuli)

## ----include = TRUE, eval = FALSE, fig.align = "center", fig.asp = .4, fig.width = 8, out.width = "100%"----
#  plot_stimuli(fit_hbam)

## ----echo = FALSE, fig.align = "center", fig.asp = .4, fig.width = 8, out.width = "100%"----
knitr::include_graphics("p_stim.svg")

## ----include = TRUE, eval = FALSE, fig.align = "center", fig.asp = .4, fig.width = 8, out.width = "100%"----
#  plot_respondents(fit_hbam, n_draws = 10)

## ----echo = FALSE, fig.align = "center", fig.asp = .4, fig.width = 8, out.width = "100%"----
knitr::include_graphics("p_resp.svg")

## ----get_plot_dat, include = TRUE, eval = FALSE, include= FALSE, eval = FALSE----
#  get_plot_data(fit_hbam)

## ----plot_alpha, include = TRUE, eval = FALSE, echo = TRUE, fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
#  plot_over_self(list(fit_hbam, fit_hbam_mini), dat, "alpha")

## ----echo = FALSE, fig.align = "center", fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
knitr::include_graphics("p_alpha.svg")

## ----plot_beta, include = TRUE, eval = FALSE, echo = TRUE, fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
#  plot_over_self(list(fit_hbam, fit_hbam_mini), dat, "abs_beta")

## ----echo = FALSE, fig.align = "center", fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
knitr::include_graphics("p_abs_beta.svg")

## ----plot_lambda, include = TRUE, eval = FALSE, echo = TRUE, fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
#  plot_over_self(list(fit_hbam, fit_hbam_mini), dat, "lambda")

## ----echo = FALSE, fig.align = "center", fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
knitr::include_graphics("p_lambda.svg")

## ----plot_chi, include = TRUE, eval = FALSE, echo = TRUE, fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
#  plot_over_self(list(fit_hbam, fit_hbam_mini), dat, "chi")

## ----echo = FALSE, fig.align = "center", fig.asp = .37, fig.width = 8.5, fig.align = "center", out.width = "100%"----
knitr::include_graphics("p_chi.svg")

## ----plot_eta, include = TRUE, eval = FALSE, echo = TRUE, fig.asp = .67, fig.width = 4.25, fig.align = "center", out.width = "50%"----
#  fit_hbam <- hbam(data = dat, extra_pars = "eta")
#  plot_over_self(fit_hbam, dat, "eta")

## ----echo = FALSE, fig.align = "center", fig.asp = .67, fig.width = 4.25, fig.align = "center", out.width = "50%"----
knitr::include_graphics("p_eta.svg")

## ----get_est_t, include = TRUE, eval = FALSE----------------------------------
#  get_est(fit_hbam, "theta")

## ----echo = FALSE-------------------------------------------------------------
est_theta

## ----get_est_v, include = TRUE, eval = FALSE----------------------------------
#  get_est(fit_hbam, "chi")

## ----echo = FALSE-------------------------------------------------------------
est_chi

## ----include = FALSE----------------------------------------------------------
load("elpds.rda")

## ----run_cv_new, eval = FALSE-------------------------------------------------
#  elpd_hbam <- hbam_cv(self, stimuli, model = "HBAM")
#  elpd_hbam_nf <- hbam_cv(self, stimuli, model = "HBAM_NF")
#  elpd_bam <- hbam_cv(self, stimuli, model = "BAM")

## ----elpd_ill, eval = FALSE, include = TRUE-----------------------------------
#  elpds <- rbind(elpd_bam, elpd_hbam, elpd_hbam_nf)
#  elpds[order(elpds$ELPD), ]

## ----echo = FALSE, include = TRUE---------------------------------------------
elpds <- rbind(elpd_bam, elpd_hbam, elpd_hbam_nf)
elpds[order(elpds$ELPD), ]

## ----eval = FALSE-------------------------------------------------------------
#  elpd_hbam_mini <- hbam_cv(self, stimuli, model = "HBAM_MINI")
#  elpd_hbam_mini

## ----echo = FALSE-------------------------------------------------------------
elpd_hbam_mini

## ----traceplot_code, eval = FALSE, include = TRUE-----------------------------
#  rstan::traceplot(fit_hbam, pars = "theta")

## ----echo = FALSE, fig.align = "center", fig.asp = .5, fig.width = 8, fig.align = "center", out.width = "100%", warning = FALSE, message = FALSE----
knitr::include_graphics("p_trace_theta.svg")

