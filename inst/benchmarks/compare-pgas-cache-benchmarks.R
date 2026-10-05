args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L)
old <- read.csv(file.path(args[1],'timing.csv'))
new <- read.csv(file.path(args[2],'timing.csv'))
rows <- lapply(old$case,function(id) {
 a <- readRDS(file.path(args[1],paste0(id,'.rds')))
 b <- readRDS(file.path(args[2],paste0(id,'.rds')))
 data.frame(case=id,posterior_identical=identical(a$posterior,b$posterior),
    predictions_identical=identical(a$prediction,b$prediction),
    diagnostics_identical=identical(a$diagnostics,b$diagnostics),rng_identical=identical(a$rng,b$rng),
    max_prediction_error=max(abs(a$prediction$draws-b$prediction$draws)),
    baseline_seconds=old$seconds[old$case==id],optimized_seconds=new$seconds[new$case==id])
})
result <- do.call(rbind,rows)
result$speedup <- result$baseline_seconds/result$optimized_seconds
write.csv(result,args[3],row.names=FALSE)
print(result,row.names=FALSE)
stopifnot(all(result$posterior_identical),all(result$predictions_identical),
          all(result$diagnostics_identical),all(result$rng_identical))
cat('PASS: all posterior states, predictions, diagnostics and RNG streams identical.\n')
