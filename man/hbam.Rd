% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hbam.R
\name{hbam}
\alias{hbam}
\title{Fit an HBAM model}
\usage{
hbam(
  self = NULL,
  stimuli = NULL,
  model = "HBAM",
  allow_miss = 2,
  req_valid = NA,
  req_unique = 2,
  prefs = NULL,
  group_id = NULL,
  data = NULL,
  pars = c("alpha", "beta", "chi", "lambda", "theta"),
  extra_pars = NULL,
  include = TRUE,
  chains = 4,
  cores = parallel::detectCores(logical = FALSE),
  warmup = 1000,
  iter = 2000,
  seed = sample.int(.Machine$integer.max, 1),
  control = list(max_treedepth = 7),
  sigma_alpha = NULL,
  sigma_beta = 0.3,
  sigma_mu_alpha = NULL,
  sigma_mu_beta = 0.2,
  ...
)
}
\arguments{
\item{self}{An optional numerical vector of N ideological self-placements. Any missing data must be coded as NA. If this argument is not supplied (either here or in a previous call to \code{prep_data()}), respondent positions will not be estimated. If the data have been prepared in advance via the \code{prep_data()} function, the argument supplied here will be ignored.}

\item{stimuli}{An N × J matrix of numerical stimulus placements, where J is the number of stimuli. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the \code{prep_data()} function.}

\item{model}{Character: Name of the model to be used. Defaults to HBAM. The available models are described under Details.}

\item{allow_miss}{Integer specifying how many missing stimulus positions to be accepted for an individual still to be included in the analysis. This argument will not be used if the data have been prepared in advance via the \code{prep_data()} function. Defaults to 2.}

\item{req_valid}{Integer specifying how many valid observations to require for a respondent to be included in the analysis. The default is \code{req_valid = J - allow_miss}, but if specified, \code{req_valid} takes precedence. This argument will not be used if the data have been prepared in advance via the \code{prep_data()} function.}

\item{req_unique}{Integer specifying how may unique positions on the ideological scale each respondent is required to have used when placing the stimuli in order to be included in the analysis. The default is \code{req_unique = 2}. This argument will not be used if the data have been prepared in advance via the \code{prep_data()} function.}

\item{prefs}{An N × J matrix of numerical stimulus ratings or preference scores. These data are only required by the HBAM_R_MINI model and will be ignored when fitting other models.}

\item{group_id}{Vector of length N identifying which group each respondent belongs to. The format can be factor, character, integer, or numeric. Respondents with NAs on \code{group_id} will be dropped when \code{group_id} is supplied. These data are only required by models with \code{"MULTI"} in their name and will be ignored when fitting other models.}

\item{data}{List of data that have been prepared in advance via the \code{prep_data()} function. Not required if the arguments \code{self} and \code{stimuli} are provided.}

\item{pars}{A vector of character strings specifying parameters of interest. If \code{include = TRUE}, only samples for parameters named in pars are stored in the fitted results. Conversely, if \code{include = FALSE}, samples for all parameters except those named in pars are stored in the fitted results. The default is to store results for "alpha", "beta", "chi", "lambda", and "theta".}

\item{extra_pars}{A vector of character strings specifying parameters or generated quantities to be added to \code{pars}. This makes it easy to add one or several additional parameters of interest, without having to repeat the default \code{pars} vector. The default is \code{NULL}.}

\item{include}{Logical scalar defaulting to \code{TRUE} indicating whether to include or exclude the parameters given by the pars argument. If \code{FALSE}, only entire multidimensional parameters can be excluded, rather than particular elements of them.}

\item{chains}{A positive integer specifying the number of Markov chains. Defaults to 4.}

\item{cores}{The number of cores to use when executing the Markov chains in parallel. By default, all detected physical cores will be used if \code{chains} is equal to or higher than the number of cores.}

\item{warmup}{A positive integer specifying the number of warmup (aka burn-in) iterations per chain. If step-size adaptation is on (which it is by default), this also controls the number of iterations for which adaptation is run (and hence these warmup samples should not be used for inference). The number of warmup iterations should be smaller than \code{iter}.}

\item{iter}{A positive integer specifying the number of iterations for each chain (including warmup).}

\item{seed}{A positive integer specifying an optional seed for reproducibility. If this argument is not supplied, a random seed will be generated and the function will produce slightly different results on each run.}

\item{control}{A named list of parameters to control the sampler's behavior. See the documentation for \code{rstan::stan} for more details.}

\item{sigma_alpha}{A positive numeric value specifying the standard deviation of the prior on the shift parameters in the FBAM model, or the standard deviation of the parameters' deviation from the group-means in FBAM_MULTI models. (This argument will be ignored by HBAM models.) Defaults to B / 5, where B measures the length of the survey scale as the number of possible placements on one side of the center.}

\item{sigma_beta}{A positive numeric value specifying the standard deviation of the prior on the logged stretch parameters in the FBAM model, or the standard deviation of the logged parameters' deviation from the group-means in FBAM_MULTI models. (This argument will be ignored by HBAM models.) Defaults to .3.}

\item{sigma_mu_alpha}{A positive numeric value specifying the standard deviation of the prior on the group-means of the shift parameters in MULTI-type models. Defaults to B / 10.}

\item{sigma_mu_beta}{A positive numeric value specifying the standard deviation of the prior on the group-means of the logged stretch parameters in MULTI-type models. Defaults to .2.}

\item{...}{Arguments passed to \code{rstan::sampling()}.}
}
\value{
An object of S4 class \code{stanfit}.
}
\description{
Fit a Hierarchical Bayesian Aldrich-McKelvey model using automatically tuned Hamiltonian Monte Carlo sampling (NUTS) via \code{rstan}.
}
\details{
This package provides several alternative models that can be selected using the names below. Users who are unsure which model to use are advised to use the default HBAM model. If speed or sampling diagnostics are an issue, HBAM_MINI may provide a useful alternative.

\strong{HBAM} is the default model, which allows for scale flipping and employs hierarchical priors on the shift and stretch parameters. It also models heteroskedastic errors that vary by both individual and stimuli. Compared to the model in Bølstad (2024), this version has been slightly revised to provide faster sampling. A key difference from the original model is that the respondent positions are not treated as parameters, but rather calculated as a function of self-placements, individual-level parameters, and simulated errors. This makes the model considerably faster, while yielding very similar results. The model simulates errors in the self-placements of the same magnitude as that with which the respondent in question places the stimulus with the smallest errors. All models in the package use this approach.

\strong{HBAM_MULTI} is a version that models differences between groups defined by the user. It requires a vector identifying the groups to be supplied as the argument \code{group_id}. The model gives each group separate hyperparameters for the locations of the prior distributions for the shift and stretch parameters. Rather than shrinking the estimates toward the mode for the whole dataset, this model shrinks the estimates toward the mode for the group. The vectors of hyperparameters are called \code{mu_alpha} and \code{mu_beta} and are constructed to have means of 0. The scales of the priors on these hyperparameters can be set by the user via the arguments \code{sigma_mu_alpha} and \code{sigma_mu_beta}. The default values are B / 10 and .2, respectively. (Here, B measures the length of the survey scale as the number of possible placements on one side of the center.) One potential use for this model is to supply self-placements as \code{group_id}, and thus give each self-placement group its own prior distribution for the shift and stretch parameters.

\strong{HBAM_NF} (formerly HBAM_0) is a version of the HBAM model that does not allow for scale flipping. This may be useful if there are truly zero cases of scale flipping in the data. Such scenarios can be created artificially, but may also arise in real data. For example, expert surveys appear unlikely to contain many instances of scale flipping. For data that contain zero cases of flipping, models that allow for flipping contain superfluous parameters that lead to inefficient sampling. Models that do not allow for flipping will sample faster and typically yield slightly more accurate estimates. Such models are therefore usually preferable when no flipping is present.

\strong{HBAM_MULTI_NF} is a version of the HBAM_MULTI model that does not allow for scale flipping.

\strong{HBAM_MINI} is a version of the HBAM model that assumes the prediction errors in the stimuli placements to be homoskedastic. This model tends to sample faster faster than the standard HBAM model while yielding very similar point estimates. For large datasets, this model may provide a reasonable compromise between model complexity and estimation speed.

\strong{FBAM} is a version of the HBAM model with fixed hyperparameters to allow fitting via optimization rather than MCMC -- which can be useful for large data sets. This model allows the user to specify the scales of the priors for the shift and (logged) stretch parameters via the arguments \code{sigma_alpha} and \code{sigma_beta}. The default values are B / 5 and .3, respectively. These defaults are intended to be realistic and moderately informative. Users who want to control the degree of shrinkage of the individual-level parameters may find it useful to fit this model -- or other FBAM models -- via either MCMC or optimization.

\strong{FBAM_MULTI} is a version of the FBAM model that shares the group-modeling features of the HBAM_MULTI model. It allows the user to set the scales of the priors for the shift and stretch parameters via the arguments \code{sigma_alpha} and \code{sigma_beta}, and set the scales of the priors on \code{mu_alpha} and \code{mu_beta} via the arguments \code{sigma_mu_alpha} and \code{sigma_mu_beta}.

\strong{FBAM_MULTI_NF} is a version of the FBAM_MULTI model that does not allow for scale flipping.

\strong{HBAM_R_MINI} is a version of the HBAM_MINI model that incorporates the rationalization component of the ISR model by Bølstad (2020). This model requires additional data to be supplied as the argument \code{pref}: An N × J matrix of stimuli ratings from the respondents. The rationalization part of the model is simplified relative to the original ISR model: The direction in which respondents move disfavored stimuli is estimated as a common expectation for each possible self-placement on the scale.

\strong{BAM} is an unpooled model with wide uniform priors on the shift and stretch parameters. It is similar to the JAGS version introduced by Hare et al. (2015), although the version included here has been adjusted to yield stretch parameters with an average of approximately one (and thus produce a scale similar to those of the other models in the package). Like the other models, this version of the BAM model also simulates errors in the self-placements to yield a realistic level of uncertainty. While this model is simple and fast, it tends to overfit the data and produce invalid posterior distributions for some respondent positions (see Bølstad 2024). However, it could potentially be useful as a baseline for model comparisons in situations where respondent positions are not of interest.

\strong{HBAM_2} has been replaced by the more general HBAM_MULTI model.

These models can also be used in situations where self-placements are not available and the only goal is to estimate stimulus positions or respondents' shift and stretch parameters. While the latent respondent positions will not be estimated, all other parameters are unaffected if the argument \code{self} is dropped when calling \code{prep_data()} or \code{hbam()}.

See the \code{hbamr} vignette for a table summarizing the key characteristics of the available models.
}
\examples{
\donttest{
# Loading and re-coding ANES 1980 data:
data(LC1980)
LC1980[LC1980 == 0 | LC1980 == 8 | LC1980 == 9] <- NA

# Making a small subset of the data for illustration:
self <- LC1980[1:100, 1]
stimuli <- LC1980[1:100, -1]

# Fitting the HBAM_MINI model, obtaining 1000 draws:
fit_hbam_mini <- hbam(self, stimuli, model = "HBAM_MINI",
                   warmup = 500, iter = 1000, chains = 2)

# Preparing the data before fitting, requiring complete responses:
dat <- prep_data(self, stimuli, allow_miss = 0)
fit_hbam_mini <- hbam(data = dat, model = "HBAM_MINI",
                   warmup = 500, iter = 1000, chains = 2)

# Obtaining posterior summaries for the latent stimulus positions:
theta_est <- get_est(fit_hbam_mini, par = "theta")

# Obtaining posterior summaries for the latent respondent positions
  # in a format matching the rows in the original dataset:
chi_est <- get_est(fit_hbam_mini, par = "chi", format_orig = TRUE)

# Fitting the FBAM_MULTI_NF model with self-placements as group_id:
fit_fbam_multi_nf <- hbam(self, stimuli, group_id = self, model = "FBAM_MULTI_NF",
                       warmup = 500, iter = 1000, chains = 2)
}
}
\references{
\itemize{
\item Bølstad, Jørgen (2024). Hierarchical Bayesian Aldrich-McKelvey Scaling. \emph{Political Analysis}. 32(1): 50–64. \doi{10.1017/pan.2023.18}.
\item Bølstad, Jørgen (2020). Capturing Rationalization Bias and Differential Item Functioning: A Unified Bayesian Scaling Approach. \emph{Political Analysis} 28(3): 340–355.
\item Hare, Christopher et al. (2015). Using Bayesian Aldrich-McKelvey Scaling to Study Citizens' Ideological Preferences and Perceptions. \emph{American Journal of Political Science} 59(3): 759–774.
}
}
