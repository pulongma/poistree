#ifndef POISTREE_MCMC_PROGRESS_H
#define POISTREE_MCMC_PROGRESS_H

#include <Rcpp.h>
#include <progress.hpp>
#include <exception>

class PPTMCMCProgress {
 public:
  PPTMCMCProgress(unsigned long iterations, bool verbose)
      : progress_(iterations, verbose),
        exceptions_(std::uncaught_exceptions()), verbose_(verbose) {}

  ~PPTMCMCProgress() {
    if (std::uncaught_exceptions() > exceptions_) {

      Progress::monitor().abort();
      if (verbose_) REprintf("\n");
    }
  }

  void increment() {
    Rcpp::checkUserInterrupt();
    progress_.increment();
  }

  PPTMCMCProgress(const PPTMCMCProgress&) = delete;
  PPTMCMCProgress& operator=(const PPTMCMCProgress&) = delete;

 private:
  Progress progress_;
  const int exceptions_;
  const bool verbose_;
};

#endif
