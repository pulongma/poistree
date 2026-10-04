\name{NEWS}
\title{News for Package \pkg{poistree}}

\section{Changes in poistree version 0.3.0 (2026-10-04)}{
  \subsection{Informed MH and candidate controls}{
    \itemize{
      \item Added \code{sampler = "irjmcmc"} for the hard and soft
        terminal-leaf models. Square-root balanced proposals include the
        reverse base-proposal probability and the full Metropolis correction.
        The sampler targets the corresponding RJ-MCMC posterior and returns
        ordinary, equally weighted posterior draws.
      \item Hard informed MH caches candidate cuts, local collapsed scores,
        and the current neighborhood across rejections. Soft informed MH
        exactly aggregates latent-allocation proposals by Poisson-binomial
        count recursion, retaining separate grow/prune and change kernels,
        allocation Gibbs sweeps, and gate updates.
      \item Every sampler now defaults to \code{cut_candidates = 50}.
        Hard Particle Gibbs accepts and forwards this control instead of
        hardcoding 30. Explicit user-specified counts and soft PGAS cut grids
        retain their existing meaning; duplicate and invalid cuts are removed.
    }
  }
  \subsection{Audit fixes}{
    \itemize{
      \item Added scalar integer depth validation and native allocation guards
        before complete-tree size calculations. Unsafe depths and oversized
        particle-tree allocations now fail with explicit errors; supported
        controls retain their existing posterior and proposal rules.
      \item Removed precompiled object files from the source checkout so direct
        installations compile the current native interfaces and samplers.
      \item Clarified the backend-specific availability of \code{min_leaf_n},
        \code{prediction_draws}, and \code{gate_family}.
      \item Unified hard SMC and Particle-Gibbs split ties with retained-state
        evaluation: ties go to the right child on every input. Split routing
        preserves observations on unrelated outer-region boundaries.
      \item Shared SMC now honors \code{max_depth = 0}; dense SMC honors
        \code{cut_candidates}, including a single-candidate grid.
      \item Single-location hard marginal queries retain every posterior draw.
      \item Hard SMC and Particle Gibbs report the conditional point-process
        log likelihood of the sampled rates, consistently with RJ-MCMC.
      \item SMC reports \code{log_relative_normalizer} and, for \code{b > 0},
        \code{log_target_normalizer} including the root marginal factor.
        Absolute \code{log_evidence} is unavailable because the existing
        tree-prior weights have an uncomputed normalizing constant. The tree
        weights and backend-specific prior choices are preserved.
      \item The internal soft-tree diagnostic updates all dimension-specific
        gate parameters, matching the public fitting backend.
    }
  }
  \subsection{Model interface}{
    \itemize{
      \item The package implements two terminal-leaf intensity
        models, PPT (\code{gating = "hard"}) and S-PPT
        (\code{gating = "soft"}, the default), with samplers
        \code{rjmcmc}, the new terminal-leaf \code{irjmcmc},
        \code{smc} (PPT only), and \code{pgas}.
        The \code{scales} argument accepts only \code{"leaf"}.
    }
  }
  \subsection{New sampler}{
    \itemize{
      \item Added Particle Gibbs with ancestor sampling for the soft
        terminal-leaf model through \code{ppt_fit(gating = "soft",
        scales = "leaf", sampler = "pgas")}. The conditional SMC runs over
        the tree and the latent allocation variables of the observations jointly on an extended
        complete-tree space with a fixed per-input cut grid, which makes the
        full suffix target ratio of the ancestor-sampling weights available.
        The forward proposal is the exact one-step-ahead posterior of the
        node decision (a Poisson-binomial expectation over the allocations of the node's observations)
        for nodes with at most \code{exact_max} points (log-space
        Poisson-binomial recursion) and a Laplace approximation of its
        Beta-integral representation otherwise; the gate is updated
        by Metropolis--Hastings and the labels by Gibbs sweeps. Resampling
        (with ancestor sampling) takes place after every heap node at which
        a particle advanced (\code{resampling = "node"}) or after every
        level (\code{"level"}), subject to an optional ESS threshold
        (\code{ess_threshold}). Above \code{exact_max} the allocation of a
        node's points to its children is proposed by sequential imputation
        (\code{allocation = "sequential"}) or through auxiliary child rates
        (\code{allocation = "rates"}); both are exact importance proposals.
      \item The hard PPT SMC sampler runs on a shared-path store
        (\code{src/hard_shared_smc.cpp}, \code{engine = "shared"}, the
        default): nodes are boxes shared by every particle that reaches
        them, candidate cuts are scored once per distinct node with sorted
        coordinates and binary-search counts, observation sets are released
        once a node is decided, and resampling copies record maps instead
        of trees. Same law as the dense sampler (\code{engine = "dense"});
        at d = 10, n = 11,000, 60 particles and depth 5 the sweep takes
        0.2 s instead of 25 s. New diagnostics \code{expanded_nodes} and
        \code{resampling_events} (NA where not applicable) in every backend.
      \item Soft models default to one gating parameter per input
        (\code{gate_structure = "dimension"}); the Metropolis update of the
        gating parameters now proposes every coordinate once per iteration
        (systematic scan) instead of one coordinate chosen at random, in
        both the RJ-MCMC and the PGAS backends. Implemented in \code{SMCtree.cpp}
        (class \code{SoftSMCtree}) and \code{soft_smc.h}.
      \item The soft particle system is stored as a shared-path trie of
        coloured path nodes (\code{SoftPathNode} in \code{soft_smc.h}):
        the expansion of a frontier node is computed once per distinct node
        and shared by all particles holding it, gate values are tabulated
        once per gate vector, and resampling copies handle maps instead of
        trees. Typically 4--10 times faster than the dense version.
      \item Regenerated \code{RcppExports}.
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
