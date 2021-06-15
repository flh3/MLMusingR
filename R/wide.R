#' Wide dataset to be used for growth modeling
#'
#' A dataset containing 30 observations with reading scores
#' taken in the fall kindergarten, spring kindergarten, and
#' spring first grade
#'
#' @format A wide data frame of 30 observations:
#' \describe{
#'   \item{studentid}{Factor indicating student identification}
#'   \item{int}{treatment or control}
#'   \item{female}{1 = female, 0 = male}
#'   \item{fall_k}{Reading scores in fall kindergarten}
#'   \item{spring_k}{Reading scores in spring kindergarten}
#'   \item{spring_g1}{Reading scores in spring first grade}
#'
#' }
"wide"
