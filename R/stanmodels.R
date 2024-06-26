# Generated by rstantools.  Do not edit by hand.

# names of stan models
stanmodels <- c("BAM", "FBAM", "FBAM_MULTI", "FBAM_MULTI_NF", "HBAM", "HBAM_MINI", "HBAM_MULTI", "HBAM_MULTI_NF", "HBAM_NF", "HBAM_R_MINI")

# load each stan module
Rcpp::loadModule("stan_fit4BAM_mod", what = TRUE)
Rcpp::loadModule("stan_fit4FBAM_mod", what = TRUE)
Rcpp::loadModule("stan_fit4FBAM_MULTI_mod", what = TRUE)
Rcpp::loadModule("stan_fit4FBAM_MULTI_NF_mod", what = TRUE)
Rcpp::loadModule("stan_fit4HBAM_mod", what = TRUE)
Rcpp::loadModule("stan_fit4HBAM_MINI_mod", what = TRUE)
Rcpp::loadModule("stan_fit4HBAM_MULTI_mod", what = TRUE)
Rcpp::loadModule("stan_fit4HBAM_MULTI_NF_mod", what = TRUE)
Rcpp::loadModule("stan_fit4HBAM_NF_mod", what = TRUE)
Rcpp::loadModule("stan_fit4HBAM_R_MINI_mod", what = TRUE)

# instantiate each stanmodel object
stanmodels <- sapply(stanmodels, function(model_name) {
  # create C++ code for stan model
  stan_file <- if(dir.exists("stan")) "stan" else file.path("inst", "stan")
  stan_file <- file.path(stan_file, paste0(model_name, ".stan"))
  stanfit <- rstan::stanc_builder(stan_file,
                                  allow_undefined = TRUE,
                                  obfuscate_model_name = FALSE)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name,
                            model_cppcode = stanfit$cppcode)
  # create stanmodel object
  methods::new(Class = "stanmodel",
               model_name = stanfit$model_name,
               model_code = stanfit$model_code,
               model_cpp = stanfit$model_cpp,
               mk_cppmodule = function(x) get(paste0("rstantools_model_", model_name)))
})
