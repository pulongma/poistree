## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2 and above, April 2020.

## Copyright (C) 2025-present by Pulong Ma

# Display the package startup and citation message.
.onAttach <- function(...) {

	date <- date()
	x <- regexpr("[0-9]{4}", date)
	this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)

	packageStartupMessage("\n#########################################################")
	packageStartupMessage("## poistree: Bayesian Poisson Point-Process Trees")
	packageStartupMessage("## Copyright (C) 2025-", this.year,
			" by Pulong Ma <plma@iastate.edu>", sep="")

	packageStartupMessage("## Please cite poistree including its version number.")
	packageStartupMessage("##########################################################")
}

# Unload the package shared library.
.onUnload <- function(libpath) {
	library.dynam.unload("poistree", libpath)
}
