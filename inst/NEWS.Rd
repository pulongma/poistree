\name{NEWS}
\title{News for Package \pkg{poistree}}

\section{Changes in poistree version 0.2.0 (2026-07-28)}{
  \itemize{
    \item Renamed the package from \pkg{PPTree} to \pkg{poistree}.
    \item Added the unified \link{ppt_fit} interface, separating gating, intensity scales, and sampler.
    \item Added hard-gated leaf PPT inference through both SMC and RJ MCMC.
    \item Restricted the current implementation to hard-gated terminal-leaf
      PPT while retaining the component-based \link{ppt_fit} wrapper for
      future models.
    \item Exposed the hard-leaf Particle Gibbs with ancestor sampling backend
      through \code{ppt_fit(..., sampler = pgas)}. Non-PPT configurations are
      rejected before sampler dispatch.
    \item Added class \link{ppt} methods and the companion functions: \link{ppt_predict}, \link{ppt_summary}, \link{ppt_logLik}, \link{ppt_evidence}, \link{ppt_lppd}, \link{ppt_diagnostics}, and \link{ppt_sim}.
    \item Retained the legacy hard-PPT functions during the API transition.
  }
}

\section{Changes in PPTree version 0.1.1 (2026-06-29): correctness fixes}{
	\itemize{
		\item \link{PPT.SMC} no longer returns \code{NULL} when \code{XX} is missing; \code{XX} now defaults to the observed points.
		\item \code{lppd} is now summed over the data points (a valid point-process predictive density) instead of over the prediction grid.
		\item The leaf marginal and the leaf-intensity draw now honour the Gamma(a,b) prior (with the Jeffreys b=0 limit) instead of hardcoding a=0.5,b=0.
		\item Removed the spurious \code{-log(n!)} term from \code{base_mloglik} (point-process likelihood) and the double-counted Beta-branch term from \code{split_mloglik}.
		\item Minor: 1/P weight initialisation and NaN guard in the SMC; boundary-consistent \code{predict_lambda}; particles serialised after leaf intensities are drawn.
		\item See \code{REVIEW.md} for details and verification.
	}
}
