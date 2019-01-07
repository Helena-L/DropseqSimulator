#' generate a countmat matrix given gene number and cell number
#' @param ngenes gene number
#' @param ncells cell number
#' @return a countmat matrix
#' @export

generate_countmat = function(ngenes, ncells) {
  params = newSplatParams()
  params = setParam(params, "nGenes", ngenes)
  params = setParam(params, "batchCells", ncells)
  params = setParam(params, "lib.loc", 1)
  params = setParam(params, "dropout.present", TRUE)
  readmat = counts(splatSimulate(params))

  return (readmat)
}
