#' generate a countmat matrix given gene number and cell number
#' @param ngenes gene number
#' @param ncells cell number
#' @param loc location parameter for the library size log-normal distribution, together with \code{libscale} controls library size. Determine different sequencing depth (number of reads) between cells.
#' @param scale scale parameter for the library size log-normal distribution, together with \code{libloc} controls library size. Determine different sequencing depth (number of reads) between ceels.
#' @return a countmat matrix
#' @export

generate_countmat = function(ngenes, ncells, loc, scale) {
  params = newSplatParams()
  params = setParam(params, "nGenes", ngenes)
  params = setParam(params, "batchCells", ncells)
  params = setParam(params, "lib.loc", loc)
  params = setParam(params, "lib.scale", scale)
  params = setParam(params, "dropout.present", TRUE)
  readmat = counts(splatSimulate(params))

  return (readmat)
}
