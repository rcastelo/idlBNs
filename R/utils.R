## from https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## function: isPackageLoaded
## purpose: to check whether the package specified by the name given in
##          the input argument is loaded. this function is borrowed from
##          the discussion on the R-help list found in this url:
##          https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## parameters: name - package name
## return: TRUE if the package is loaded, FALSE otherwise
          
.isPackageLoaded <- function(name) {
    (paste("package:", name, sep="") %in% search()) ||
    (name %in% loadedNamespaces())
}

#' @importFrom utils installed.packages
.load_suggested_package <- function(pkgname) {
    instpkgs <- installed.packages(noCache=TRUE)[, "Package"]
    installed <- pkgname %in% instpkgs
    loaded <- .isPackageLoaded(pkgname)
    if (!loaded)
        loaded <- suppressPackageStartupMessages(requireNamespace(pkgname,
                                                                  quietly=TRUE))
    return(installed & loaded)
}
