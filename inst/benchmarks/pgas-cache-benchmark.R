args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
if(args[1]!='default') .libPaths(c(normalizePath(args[1]),.libPaths()))
library(poistree)
out <- args[2]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(127)
x <- matrix(runif(240*3),240,3)
u <- x[,1]; left <- u < .5
x[left,1] <- qbeta(2*u[left],3,11)
x[!left,1] <- qbeta(2*u[!left]-1,11,3)
colnames(x) <- c('x','y','z')
cases <- list(
 fixed_sparse=list(n=240,d=3,particles=12,iter=120,burn=40,thin=2,max_depth=10,cut_candidates=8,gate=8,update_gate=FALSE,exact_max=40),
 dimension_update=list(n=150,d=3,particles=8,iter=70,burn=20,thin=2,max_depth=6,cut_candidates=7,gate=c(5,8,12),update_gate=TRUE,exact_max=25),
 shared_update=list(n=100,d=2,particles=8,iter=60,burn=20,thin=2,max_depth=5,cut_candidates=6,gate=8,gate_structure='shared',update_gate=TRUE,exact_max=150),
 level_exact=list(n=100,d=2,particles=8,iter=60,burn=20,thin=2,max_depth=8,cut_candidates=6,gate=8,update_gate=FALSE,exact_max=150,resampling='level'),
 adaptive_approx=list(n=100,d=2,particles=8,iter=60,burn=20,thin=2,max_depth=6,cut_candidates=6,gate=c(5,10),update_gate=TRUE,exact_max=0,ess_threshold=.7),
 pg_rates=list(n=90,d=2,particles=6,iter=50,burn=10,thin=2,max_depth=5,cut_candidates=6,gate=8,update_gate=FALSE,exact_max=0,ancestor_sampling=FALSE,allocation='rates'),
 one_level=list(n=40,d=2,particles=6,iter=30,burn=10,thin=2,max_depth=1,cut_candidates=4,gate=8,update_gate=TRUE,exact_max=150,resampling='level'))
timings <- list()
for (i in seq_along(cases)) {
 name <- names(cases)[i]; g <- cases[[i]]; n <- g$n; d <- g$d; g$n <- g$d <- NULL
 X <- x[seq_len(n),seq_len(d),drop=FALSE]
 region <- cbind(lower=rep(0,d),upper=rep(1,d))
 common <- list(x=X,region=region,predict_at=X[seq(1,n,length.out=25),,drop=FALSE],
                test=X[seq(2,n,length.out=12),,drop=FALSE],gating='soft',scales='leaf',sampler='pgas',
                a=.5,b=.5/n,a_gate=4,b_gate=4/g$gate,sd_gate=.15,
                alpha=.85,eta=.7,label_sweeps=1,seed=8700+i,verbose=FALSE)
 started <- proc.time()
 fit <- do.call(ppt_fit,c(common,g))
 elapsed <- proc.time()-started
 saveRDS(list(posterior=fit$posterior,prediction=fit$prediction,diagnostics=fit$diagnostics,
              rng=.Random.seed,control=fit$control),file.path(out,paste0(name,'.rds')))
 timings[[i]] <- data.frame(case=name,n=n,d=d,iter=g$iter,seconds=elapsed[['elapsed']],
                            cpu_seconds=elapsed[['user.self']]+elapsed[['sys.self']],
                            leaves=fit$posterior$mean_leaves,expanded=fit$diagnostics$expanded_nodes,
                            resampled=fit$diagnostics$resampling_events)
 write.csv(do.call(rbind,timings),file.path(out,'timing.csv'),row.names=FALSE)
 cat(name,elapsed[['elapsed']],'seconds\n'); flush.console()
}
saveRDS(list(x=x,cases=cases),file.path(out,'design.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
