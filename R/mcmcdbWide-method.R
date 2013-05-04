#' @name McmcdbWide-methods
#' @rdname McmcdbWide-methods
#' @title Create \code{McmcdbWide} objects
#'
#' @description
#' Create \code{\linkS4class{McmcdbWide}} objects from \code{\linkS4class{stanfit}}
#' objects produced by \pkg{stan}.
#'
#' @param object An object for which a method is available.
#' @return An object of class \code{\linkS4class{McmcdbWide}}.
#' @examples
#' \dontrun{
#' # Convert stanfit object
#' library(rstan)
#' scode <- "
#'      parameters {
#'        real y[2]; 
#'      } 
#'      model {
#'        y[1] ~ normal(0, 1);
#'        y[2] ~ double_exponential(0, 2);
#'      } 
#'      "
#' fit1 <- stan(model_code = scode, iter = 10, verbose = FALSE)
#' fit2 <- McmcdbWide(fit1)
#' }
NULL

McmcdbWide.stanfit <- function(object) {
  samples <-
    do.call(rbind, llply(object@sim[["samples"]],
                         function(y) do.call(cbind, y)))

  chains <-
    ldply(object@sim[["samples"]],
          function(DF) {
            y <- as.data.frame(Filter(Negate(is.null), attr(DF, "args")))
            y[["adaptation_info"]] <- attr(DF, "adaptation_info")
            y
          })

  iters <- 
    mdply(chains[ , c("chain_id", "iter_save", "warmup")],
          function(chain_id, iter_save, warmup) {
            data.frame(chain_id = chain_id,
                       iter = seq_len(iter_save),
                       warmup = (seq_len(iter_save) <= warmup))
          })

  sampler_params <- 
    ldply(object@sim[["samples"]],
          function(object) {
            as.data.frame(attr(object, "sampler_params"))
          })
  iters <- cbind(iters, sampler_params)

  flatpar_chains <-
    expand.grid(flatpar = as.character(colnames(samples)),
                chain_id = chains[["chain_id"]])
  inits <- ldply(seq_along(object@inits),
                 function(i) {
                   init <- mcmcdb_flatten(object@inits[[i]],
                                          FUN = mcmc_parnames_bugs)
                   data.frame(chain_id = chains[["chain_id"]][i],
                              flatpar = names(init),
                              init = unname(init))
                 })
  flatpar_chains <- merge(flatpar_chains, inits, all.object=TRUE)
  
  metadata <- list()
  metadata[["model_name"]] <- object@model_name
  metadata[["date"]] <- object@date
  metadata[["stanmodel"]] <- object@stanmodel

  # stanfit objects use BUGS-style names for some reason
  McmcdbWide(samples,
             parameters = mcmc_parparser_bugs,
             flatpar_chains = McmcdbFlatparChains(flatpar_chains),
             chains = McmcdbChains(chains),
             iters = McmcdbIters(iters),
             metadata = metadata)
}

#' @rdname McmcdbWide-methods
#' @aliases McmcdbWide,stanfit-method
setMethod("McmcdbWide", "stanfit", McmcdbWide.stanfit)
