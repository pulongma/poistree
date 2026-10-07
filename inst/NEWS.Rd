\name{NEWS}
\title{News for Package \pkg{poistree}}

\section{Changes in poistree version 0.5.8 (2026-10-07)}{

  \subsection{Samplers and proposals}{
    \itemize{
      \item Added cached informed Metropolis--Hastings tree updates through
        \code{sampler = "irjmcmc"} for hard and soft terminal-leaf models.
        Forward and reverse proposal probabilities enter the acceptance ratio;
        retained draws have ordinary posterior weights.
      \item Added soft Particle Gibbs with ancestor sampling
        (\code{sampler = "pgas"}), jointly updating trees and allocations.
        Small-node proposals sum over allocation counts exactly; larger nodes
        use a Laplace scorer with sequential or auxiliary-rate allocation
        proposals and full importance correction. Node- or level-wise
        resampling and an optional ESS threshold are supported.
      \item Added partially collapsed Gibbs (\code{sampler = "pcg"}) for
        soft trees. Gates are updated jointly with labels integrated out,
        followed by label restoration and collapsed tree updates. Robust
        adaptive Metropolis controls and gate-adaptation diagnostics are
        available; adaptation stops no later than the end of burn-in.
      \item Added optional hard-surrogate proposals for soft PCG
        (\code{informed = TRUE}) and large-node soft PGAS
        (\code{proposal_score = "hard"}). Both use \code{proposal_temperature}
        and \code{proposal_defensive}, with sampler-specific mixture laws and
        exact proposal corrections. Defaults remain uninformed PCG and
        Laplace-scored PGAS.
      \item All samplers default to \code{cut_candidates = 50}.
        The supported intensity scale is \code{scales = "leaf"} for both
        hard PPT and soft S-PPT.
    }
  }

  \subsection{Gates and physical-domain integration}{
    \itemize{
      \item Soft models default to dimension-specific gate parameters;
        shared gates remain available. Soft RJ-MCMC and PGAS update every
        gate coordinate once per iteration.
      \item Added \code{gate_scale = "node"} for logistic gates with
        RJ-MCMC, informed MH, and PCG. Existing logistic calls retain root
        scaling; compact gates retain node scaling. Prediction retains the
        scale choice, and earlier saved fits remain supported. PGAS requires
        root-scaled logistic gates.
      \item Added \code{ppt_fit_quadrature()} for intensity models using
        background covariates and physical exposure weights. Supported
        samplers are soft PCG, including informed proposals, and hard shared
        SMC. Weighted integrals are retained for post-fit scoring.
      \item Improved stability of logistic gate tails, conditional routing,
        and exposure calculations. Node-scaled logistic paths with unequal
        slopes use adaptive numerical integration; weighted background
        memberships use linear-scale products for exposure calculations.
    }
  }

  \subsection{Computation and prediction}{
    \itemize{
      \item Hard SMC uses shared-path storage by default; the dense engine
        remains available. Soft PGAS also shares node computations across
        particles, avoiding repeated expansion and full-tree copies.
      \item Exact soft PGAS reuses Poisson-binomial count distributions
        across cuts by exponential tilting. It shares triangular conditional
        allocation tables, skips deterministic tables, reuses root scores
        with coordinate-specific gate invalidation, and uses the cancelled
        exact importance ratio for ordinary and reference particles.
      \item The optional hard PGAS scorer uses cumulative split counts and
        cached path exposures, evaluating true soft gates and exposures only
        for selected or forced splits.
      \item Soft fitting reuses ancestor memberships, leaf exposures, and
        candidate cuts. PCG reuses event normalizers when restoring labels;
        grow/change proposals reuse parent memberships. Weighted-integration
        fits also reuse memberships for exactly matching background rows.
      \item Post-fit prediction shares ancestor calculations and reuses
        retained draws at matching, reordered, or repeated locations.
        Mean-only predictions skip quantiles, and requested quantiles share
        a sort while retaining the weighted quantile convention.
      \item Held-out scoring adds intensity \code{scale} and
        \code{type = "plugin"}; stored test patterns can be rescored.
        Soft PGAS retains a log-likelihood trace for mixing diagnostics.
    }
  }

  \subsection{Robustness and diagnostics}{
    \itemize{
      \item Added integer depth checks and native allocation guards, and
        removed stale precompiled objects from source installations.
      \item Corrected hard SMC/Particle-Gibbs routing at split ties and
        outer-region boundaries, root-only shared SMC fits, dense SMC cut
        counts, and single-location hard marginal draw retention.
      \item Hard SMC and Particle Gibbs report the conditional point-process
        log likelihood of sampled rates. SMC distinguishes relative and
        target normalizers; absolute Bayesian evidence remains unavailable
        because the existing tree-prior normalization is unknown.
    }
  }

  \subsection{Interface and documentation}{
    \itemize{
      \item Reduced the exported interface to \code{ppt_fit},
        \code{ppt_fit_quadrature}, \code{ppt_predict}, \code{ppt_sim},
        \code{ppt_lambda}, \code{ppt_lppd}, \code{ppt_marginal}, and
        \code{ppt_diagnostics}.
        Marginal prediction remains public for existing analysis notebooks.
        Sampler diagnostics and diagnostic plotting remain public.
      \item Made \code{ppt_logLik}, \code{ppt_integral},
        \code{ppt_quadrature_geometry}, and \code{ppt_summary} internal.
        Removed the \code{logLik.ppt}
        forwarding method and S3 registration; \code{summary(fit)} remains
        available.
      \item Removed redundant tests and consolidated overlapping cases,
        reducing declared test cases from 148 to 142 while retaining
        numerical, posterior, boundary, and cache regression coverage.
      \item Expanded argument documentation with types, defaults,
        restrictions, and organized option lists. Added \code{ppt_controls},
        a grouped reference for the 48 fitting controls passed through
        \code{ppt_fit(...)}.
      \item Clarified backend-specific controls, weighted-domain integration,
        and prediction quantiles. Removed implementation comments while
        retaining brief function-purpose notes and required metadata.
    }
  }
}

\section{Changes in poistree version 0.2.1 (2026-08-25)}{
  \subsection{New models}{
    \itemize{
      \item Added the soft-gated terminal-leaf model S-PPT, fitted by RJ-MCMC through
        the unified \link{ppt_fit} interface.
    }
  }
  \subsection{New functions and interface}{
    \itemize{
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
        integration for hard and soft terminal-leaf fits. A
        user-supplied \code{reference} matrix remains available as an explicit
        uniform reference-design approximation.
      \item Removed the hard-coded aspect-ratio limit of 6 from hard-leaf PPT
        SMC and PGAS. The new \code{max_aspect_ratio} control defaults to
        \code{Inf}; finite values apply the same shape check to ordinary,
        midpoint, and median-fallback cut candidates.
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
