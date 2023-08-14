tmp.install.packages <- function(pack, dependencies=TRUE, ...) {
  path <- tempdir()
  ## Add 'path' to .libPaths, and be sure that it is not
  ## at the first position, otherwise any other package during
  ## this session would be installed into 'path'
  firstpath <- .libPaths()[1]
  .libPaths(c(firstpath, path))
  install.packages(pack, dependencies=dependencies, lib=path, ...)
}

tmp.install.packages("lasso2_1.2-22.tar.gz", type="source", repos=NULL)