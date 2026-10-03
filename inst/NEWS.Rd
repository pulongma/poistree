\name{NEWS}
\title{News for Package \pkg{poistree}}

\section{Changes in poistree version 0.2.1 (2026-08-25)}{
  \subsection{New models}{
    \itemize{
      \item Added the soft-gated terminal-leaf model S-PPT and the hard and
        soft multiscale models MPPT and S-MPPT, all fitted by RJ-MCMC through
        the unified \link{ppt_fit} interface.
      \item Multiscale models accept \code{scale_prior = "independent"} or
        \code{"markov"}. Independent-scale S-MPPT is now the default model,
        and setting only \code{scales = "leaf"} selects S-PPT. Hard-gated PPT
        and MPPT and Markov-scale variants remain available through explicit
        component choices.
    }
  }
  \subsection{New functions and interface}{
    \itemize{
      \item Added exact informed reversible-jump MCMC for independent-scale
        MPPT and S-MPPT through \code{sampler = "irjmcmc"}. The kernel always
        combines collapsed-score informed cuts with sequential conditional
        allocation proposals. All forward and reverse proposal probabilities
        enter the Metropolis--Hastings ratio.
      \item Added \link{ppt_lambda} to evaluate the posterior intensity
        \eqn{\lambda(s)} of any fit at arbitrary locations or over an
        automatic grid (one and two dimensions), returning
        \pkg{ggplot2}-ready pointwise summaries or the raw
        locations-by-draws matrix.
      \item The default minimum leaf occupancy \code{min_leaf_n} is now 1
        for every backend (previously 5), so the prior alone governs tree
        depth unless the user asks for a stricter cut constraint.
      \item \link{ppt_predict} now accepts \code{newdata} anywhere in the
        fitted region; listing locations in \code{predict_at} before fitting
        is no longer required.
      \item \link{ppt_lppd} gained a \code{test} argument for held-out
        log pointwise predictive density evaluation after fitting.
      \item Extended \link{ppt_marginal} with no-reference posterior-state
        integration for hard and soft leaf and multiscale fits. A
        user-supplied \code{reference} matrix remains available as an explicit
        uniform reference-design approximation.
      \item Removed the hard-coded aspect-ratio limit of 6 from hard-leaf PPT
        SMC and PGAS. The new \code{max_aspect_ratio} control defaults to
        \code{Inf}; finite values apply the same shape check to ordinary,
        midpoint, and median-fallback cut candidates.
      \item Added the \code{soft-tree-surface} demo: a complete two-dimensional
        S-MPPT workflow using post-hoc \link{ppt_lambda} evaluation without
        \code{predict_at} at fitting time.
    }
  }
  \subsection{Corrections}{
    \itemize{
      \item Replaced the previous invalid conditional-SMC ancestor replay with
        an exact retained-reference Particle-Gibbs kernel. The legacy
        \code{sampler = "pgas"} token is preserved for compatibility, but
        ancestor sampling is disabled until a full tree-suffix backward weight
        is implemented.
      \item Corrected nondefault Gamma-rate handling, preservation of supplied
        nonunit regions and SMC resampling thresholds, and serialization of
        split probabilities.
    }
  }
}

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
    \item Added class \link{ppt} methods and the companion functions: \link{ppt_predict}, \link{ppt_summary}, \link{ppt_logLik}, \code{ppt_evidence}, \link{ppt_lppd}, \link{ppt_diagnostics}, and \link{ppt_sim}.
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
