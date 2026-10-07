test_that("hard informed transitions preserve the RJ posterior exactly", {
  designs <- list(
    matrix(c(.02,.04,.07,.1,.2,.7,.86,.95), ncol=1),
    cbind(c(.04,.12,.23,.32,.63,.71,.84,.96),
          c(.09,.19,.25,.37,.69,.61,.79,.91)),
    cbind(c(0,.3,.6,1), c(0,0,0,0))
  )
  for (x in designs) for (b in c(0, .1)) {
    exact <- informed_exact(x, b=b)
    keys <- vapply(exact$states, function(s) exact$key(s$tree), "")
    P <- Q <- matrix(0, length(keys), length(keys))
    for (i in seq_along(keys)) {
      law <- poistree:::PPT_IMCMC_transition(
        x, exact$region, exact$splits(exact$states[[i]]$tree),
        max_depth=2L, min_leaf_n=1L, cut_grid_n=3L,
        a=.5, b=b, alpha=.65, eta=1)
      expect_equal(law$log_target, exact$logtarget[i], tolerance=1e-11)
      expect_equal(law$log_normalizer, log(exact$hard$Z[i]), tolerance=1e-11)
      for (j in seq_along(law$trees)) {
        mat <- as.matrix(law$trees[[j]]); tree <- list()
        for (k in seq_len(nrow(mat))) tree[[as.character(mat[k,1])]] <- c(mat[k,2]+1,mat[k,3])
        at <- match(exact$key(tree), keys)
        expect_false(is.na(at))
        P[i,at] <- P[i,at] + law$neighbors$transition_probability[j]
        Q[i,at] <- Q[i,at] + exp(law$neighbors$log_q0[j])
      }
      P[i,i] <- P[i,i] + law$self_probability
    }
    expect_equal(Q, exact$hard$Q, tolerance=1e-12)
    expect_equal(P, exact$hard$P, tolerance=1e-11)
    expect_equal(rowSums(P), rep(1,nrow(P)), tolerance=1e-12)
    flow <- exact$pi*P
    expect_equal(flow, t(flow), tolerance=1e-12)
  }
})

test_that("soft informed count mixtures match enumerated labeled transitions", {
  designs <- list(matrix(c(.1,.3,.6,.9), ncol=1),
                  cbind(c(.1,.3,.6), c(.5,.5,.5)))
  for (x in designs) for (family in c("logistic", "compact")) {
    exact <- informed_exact(x, soft=TRUE, family=family)
    keys <- vapply(exact$states, function(s) exact$key(s$tree,s$labels), "")
    for (kind in 0:1) {
      ref <- if (kind == 0L) exact$gp else exact$change
      K <- Q <- matrix(0,length(keys),length(keys)); Z <- numeric(length(keys))
      for (i in seq_along(keys)) {
        s <- exact$states[[i]]
        law <- poistree:::ppstree_informed_transition(
          x, exact$region, exact$splits(s$tree), as.integer(s$labels),
          rep(5,ncol(x)), a=.5,b=.1,alpha=.65,eta=1,Dmax=2L,
          nmin=1L,cut_mode=1L,ncand=3L,
          gate_family=as.integer(family == "compact"),kind=kind)
        Z[i] <- exp(law$log_normalizer)
        for (neighbor in law$neighbors) {
          mat <- as.matrix(neighbor$splits); tree <- list()
          for (j in seq_len(nrow(mat))) tree[[as.character(mat[j,1])]] <- c(mat[j,2]+1,mat[j,3])
          at <- match(exact$key(tree,as.integer(neighbor$labels)),keys)
          expect_false(is.na(at))
          K[i,at] <- K[i,at]+neighbor$probability
          Q[i,at] <- Q[i,at]+exp(neighbor$log_q0)
        }
      }
      expect_equal(Q,ref$Q,tolerance=1e-11)
      expect_equal(Z,ref$Z,tolerance=1e-8)
      expect_equal(K,ref$K,tolerance=1e-8)
      P <- K*pmin(1,outer(Z,Z,"/")); P[!is.finite(P)] <- 0
      diag(P) <- diag(P)+1-rowSums(P)
      expect_equal(P,ref$P,tolerance=1e-8)
      flow <- exact$pi*P
      expect_equal(flow,t(flow),tolerance=1e-9)
    }
  }
})

test_that("informed samplers return unweighted draws from enumerated targets", {
  x <- matrix(c(.1,.3,.6,.9), ncol=1)
  region <- matrix(c(0,1),1)
  for (gating in c("hard","soft")) {
    exact <- informed_exact(x,soft=gating == "soft")
    want <- sum(exact$pi*vapply(exact$states,function(s) length(s$info$leaves),0L))
    args <- list(x=x,region=region,gating=gating,sampler="irjmcmc",
                 max_depth=2L,cut_candidates=3L,a=.5,b=.1,alpha=.65,eta=1,
                 chains=1L,iter=13000L,burn=1000L,seed=77L,verbose=FALSE)
    if (gating == "soft") {
      args <- c(args,list(gate=5,update_gate=FALSE,thin=1L,
                          tree_moves=1L,change_moves=1L))
    } else args$prediction_draws <- 12000L
    fit <- do.call(ppt_fit,args)
    trace <- fit$diagnostics$leaf_count_trace
    batches <- colMeans(matrix(trace,nrow=200L))
    mcse <- stats::sd(batches)/sqrt(length(batches))
    expect_lt(abs(mean(trace)-want),6*mcse+.005)
    expect_length(fit$posterior$particle_weights,0L)
    expect_equal(fit$prediction$mean,rowMeans(fit$prediction$draws),tolerance=1e-12)
    if (gating == "soft") {

      want_count <- sum(exact$pi*vapply(exact$states,function(s)
        sum(tabulate(match(s$labels,names(s$info$leaves)),length(s$info$leaves))^2),0.0))
      count_trace <- vapply(fit$posterior$state$nodes,function(st)
        sum(st[st[,2]<0,6]^2),0.0)
      batches <- colMeans(matrix(count_trace,nrow=200L))
      mcse <- stats::sd(batches)/sqrt(length(batches))
      expect_lt(abs(mean(count_trace)-want_count),6*mcse+.01)
    }
  }
})

test_that("soft reverse probabilities remain correct with multiple terminal splits", {
  x <- matrix(c(.1,.25,.4,.6,.75,.9),ncol=1)
  region <- matrix(c(0,1),1)
  splits <- rbind(c(1,0,.6),c(2,0,.25),c(3,0,.75))
  labels <- c(4L,5L,6L,7L,5L,6L)
  reference <- informed_exact(x,soft=TRUE,ncand=4L,geometry_only=TRUE)
  decode <- function(mat) {
    t <- list()
    for (k in seq_len(nrow(mat))) t[[as.character(mat[k,1])]] <- c(mat[k,2]+1,mat[k,3])
    t
  }
  original_key <- reference$key(decode(splits),labels)
  original_target <- reference$logtarget(decode(splits),labels)
  law <- function(s,z,kind) poistree:::ppstree_informed_transition(
    x,region,s,as.integer(z),5,a=.5,b=.1,alpha=.65,eta=1,Dmax=2L,
    nmin=1L,cut_mode=1L,ncand=4L,gate_family=0L,kind=kind)
  for (kind in 0:1) {
    current <- law(splits,labels,kind)

    self <- vapply(current$neighbors,function(z)
      reference$key(decode(z$splits),as.integer(z$labels)) == original_key,TRUE)
    if (kind == 1L) expect_equal(sum(self),2L)
    expect_equal(sum(vapply(current$neighbors,`[[`,0.0,"probability")),1,tolerance=1e-12)
    nonself <- current$neighbors[!self]
    for (next_state in nonself) {
      reverse <- law(next_state$splits,next_state$labels,kind)
      back <- Filter(function(z) reference$key(decode(z$splits),as.integer(z$labels)) == original_key,
                     reverse$neighbors)
      expect_length(back,1L)
      dtarget <- reference$logtarget(decode(next_state$splits),next_state$labels)-original_target
      expect_equal(next_state$log_ratio,dtarget+back[[1]]$log_q0-next_state$log_q0,tolerance=1e-10)
      p <- next_state$probability*min(1,exp(current$log_normalizer-reverse$log_normalizer))
      pr <- back[[1]]$probability*min(1,exp(reverse$log_normalizer-current$log_normalizer))
      expect_equal(p,exp(dtarget)*pr,tolerance=1e-10)
    }
  }
})
