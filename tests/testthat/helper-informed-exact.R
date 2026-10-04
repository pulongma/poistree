# Independent, tiny-state reference for the hard and augmented soft RJ targets.
# It deliberately uses neither native transition ratios nor native cut builders.
informed_exact <- function(x, max_depth = 2L, ncand = 3L, soft = FALSE,
                           gate = rep(5, ncol(x)), family = "logistic",
                           a = .5, b = .1, alpha = .65, eta = 1, geometry_only = FALSE) {
  d <- ncol(x); n <- nrow(x); region <- cbind(rep(0, d), rep(1, d))
  key_tree <- function(t) {
    if (!length(t)) return("root")
    ids <- sort(as.integer(names(t)))
    paste(vapply(ids, function(id) sprintf("%d:%d:%.17g", id,
      t[[as.character(id)]][1], t[[as.character(id)]][2]), ""), collapse = ";")
  }
  key <- function(t, z = integer()) paste(key_tree(t), paste(z, collapse = ","))
  cuts <- function(idx, box) lapply(seq_len(d), function(j) {
    s <- sort(x[idx, j]); m <- length(s)
    if (m < 2L) return(numeric())
    c <- unique(s[pmax(1L, floor((m - 1) * seq(.05, .95, length.out = ncand) + 1))])
    c[vapply(c, function(v) v > box[j, 1] + .001 &&
      v < box[j, 2] - .001 && sum(s < v) >= 1L && sum(s >= v) >= 1L, TRUE)]
  })
  rho <- function(depth) alpha * (1 + depth)^(-eta)
  logQ <- function(m, v) {
    if (b == 0) return(lgamma(a + m) - (a + m) * log(v))
    lgamma(a + m) - lgamma(a) + a * log(b) - (a + m) * log(b + v)
  }
  gate_value <- function(z, axis, cut, side, width) {
    p <- if (family == "logistic") plogis(gate[axis] * (z - cut)) else {
      # Compact gates scale their transition band by the parent width.
      u <- pmax(-1, pmin(1, gate[axis] * (z - cut) / width))
      .5 + .75 * u - .25 * u^3
    }
    if (side < 0) 1 - p else p
  }
  detail <- function(t) {
    nodes <- list(); lp <- 0
    visit <- function(id, idx, box, depth, path) {
      tag <- as.character(id); vc <- cuts(idx, box)
      legal <- depth < max_depth && any(lengths(vc) > 0L)
      split <- t[[tag]]
      nodes[[tag]] <<- list(id = id, idx = idx, box = box, depth = depth,
                           vc = vc, legal = legal, leaf = is.null(split), path = path)
      if (is.null(split)) lp <<- lp + if (legal) log1p(-rho(depth)) else 0
      else {
        j <- as.integer(split[1]); c <- split[2]
        stopifnot(legal, c %in% vc[[j]])
        lp <<- lp + log(rho(depth)) - log(d) - log(length(vc[[j]]))
        bl <- br <- box; bl[j, 2] <- c; br[j, 1] <- c
        width <- box[j, 2] - box[j, 1]
        visit(2L * id, idx[x[idx, j] < c], bl, depth + 1L,
              rbind(path, c(j, c, -1, width)))
        visit(2L * id + 1L, idx[x[idx, j] >= c], br, depth + 1L,
              rbind(path, c(j, c, 1, width)))
      }
    }
    visit(1L, seq_len(n), region, 0L, matrix(numeric(), 0, 4))
    leaves <- Filter(function(v) v$leaf, nodes)
    G <- Filter(function(v) v$legal, leaves)
    P <- Filter(function(v) !v$leaf && nodes[[as.character(2L*v$id)]]$leaf &&
                  nodes[[as.character(2L*v$id+1L)]]$leaf, nodes)
    H <- vapply(leaves, function(v) {
      if (!soft) return(prod(v$box[, 2] - v$box[, 1]))
      prod(vapply(seq_len(d), function(j) {
        path <- v$path[v$path[, 1] == j, , drop = FALSE]
        if (!nrow(path)) return(1)
        integrate(function(xx) {
          ans <- rep(1, length(xx))
          for (k in seq_len(nrow(path))) ans <- ans * gate_value(xx, j,
              path[k, 2], path[k, 3], path[k, 4])
          ans
        }, 0, 1, rel.tol = 1e-10, subdivisions = 300L)$value
      }, 0.0))
    }, 0.0)
    phi <- vapply(leaves, function(v) {
      out <- rep(1, n)
      for (k in seq_len(nrow(v$path))) {
        g <- v$path[k, ]; out <- out * gate_value(x[, g[1]], g[1], g[2], g[3], g[4])
      }
      out
    }, numeric(n))
    list(nodes = nodes, leaves = leaves, G = G, P = P, lp = lp, H = H,
         phi = matrix(phi, n, length(leaves), dimnames = list(NULL, names(leaves))))
  }
  if (geometry_only) return(list(detail = detail, key = key,
    logtarget = function(t, z) {
      info <- detail(t); ids <- names(info$leaves)
      m <- tabulate(match(z, ids), length(ids))
      info$lp + sum(logQ(m, info$H)) +
        sum(log(info$phi[cbind(seq_len(n), match(z, ids))]))
    }))
  trees <- list(list()); keys <- key_tree(trees[[1]]); details <- list(); k <- 1L
  while (k <= length(trees)) {
    info <- detail(trees[[k]]); details[[k]] <- info
    for (v in info$G) for (j in seq_len(d)) for (c in v$vc[[j]]) {
      t <- trees[[k]]; t[[as.character(v$id)]] <- c(j, c); ky <- key_tree(t)
      if (!ky %in% keys) { trees[[length(trees)+1L]] <- t; keys <- c(keys, ky) }
    }
    k <- k + 1L
  }
  states <- list(); lw <- numeric()
  for (k in seq_along(trees)) {
    info <- details[[k]]; ids <- as.integer(names(info$leaves))
    labels <- if (soft) as.matrix(expand.grid(rep(list(ids), n))) else matrix(integer(), 1, 0)
    for (r in seq_len(nrow(labels))) {
      z <- as.integer(labels[r, ])
      if (soft) {
        m <- tabulate(match(z, ids), length(ids))
        spatial <- sum(log(info$phi[cbind(seq_len(n), match(z, ids))]))
      } else {
        m <- lengths(lapply(info$leaves, `[[`, "idx")); spatial <- 0
      }
      logp <- info$lp + sum(logQ(m, info$H)) + spatial
      if (!is.finite(logp)) next
      states[[length(states)+1L]] <- list(tree = trees[[k]], labels = z, info = info)
      lw <- c(lw, logp)
    }
  }
  state_keys <- vapply(states, function(s) key(s$tree, s$labels), "")
  size <- length(states); Qgp <- Qc <- matrix(0, size, size)
  for (i in seq_along(states)) {
    s <- states[[i]]; info <- s$info
    add <- function(t, z, p, change) {
      to <- match(key(t, z), state_keys)
      if (!is.na(to)) {
        if (change) Qc[i, to] <<- Qc[i, to] + p else Qgp[i, to] <<- Qgp[i, to] + p
      }
    }
    split <- function(v, j, c, prob, change) {
      t <- s$tree; t[[as.character(v$id)]] <- c(j, c)
      if (!soft) { add(t, integer(), prob, change); return() }
      ii <- which(s$labels %in% if (change) c(2L*v$id, 2L*v$id+1L) else v$id)
      allocs <- as.matrix(expand.grid(rep(list(c(0L, 1L)), length(ii))))
      if (!length(ii)) allocs <- matrix(integer(), 1, 0)
      pr <- gate_value(x[ii, j], j, c, 1, v$box[j, 2]-v$box[j, 1])
      for (r in seq_len(nrow(allocs))) {
        zz <- allocs[r, ]; p <- prod(ifelse(zz == 1L, pr, 1-pr))
        if (p == 0) next
        z <- s$labels; z[ii] <- 2L*v$id + zz; add(t, z, prob*p, change)
      }
    }
    for (v in info$G) for (j in seq_len(d)) for (c in v$vc[[j]])
      split(v, j, c, 1/(2*length(info$G)*d*length(v$vc[[j]])), FALSE)
    for (v in info$P) {
      t <- s$tree; t[[as.character(v$id)]] <- NULL
      z <- s$labels; z[z %in% c(2L*v$id, 2L*v$id+1L)] <- v$id
      add(t, z, 1/(2*length(info$P)), FALSE)
      for (j in seq_len(d)) for (c in v$vc[[j]])
        split(v, j, c, 1/(length(info$P)*d*length(v$vc[[j]])), TRUE)
    }
  }
  pi <- exp(lw - max(lw)); pi <- pi/sum(pi)
  # Hard kernel combines three equally weighted families and excludes self.
  Q <- (2*Qgp + Qc)/3
  diag(Q) <- 0
  balanced <- function(q) {
    E <- sqrt(q * t(q) * exp(outer(lw, lw, function(x, y) y-x)))
    Z <- rowSums(E); K <- matrix(0, size, size); good <- Z > 0
    K[good, ] <- E[good, , drop = FALSE]/Z[good]
    P <- K * pmin(1, outer(Z, Z, "/")); P[!is.finite(P)] <- 0
    diag(P) <- diag(P) + 1-rowSums(P)
    list(Q = q, K = K, Z = Z, P = P)
  }
  list(states = states, pi = pi, logtarget = lw, hard = balanced(Q),
       gp = balanced(Qgp), change = balanced(Qc), x = x, region = region,
       control = list(a=a,b=b,alpha=alpha,eta=eta,max_depth=max_depth,ncand=ncand,gate=gate),
       splits = function(t) {
         if (!length(t)) return(matrix(numeric(), 0L, 3L))
         ids <- sort(as.integer(names(t)))
         t(vapply(ids, function(id) c(id,t[[as.character(id)]][1]-1,t[[as.character(id)]][2]), numeric(3)))
       }, key = key)
}
