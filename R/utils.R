#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#' @keywords internal
#' @importFrom rlang %||%
NULL


#' Rank inverse normal transformation
#'
#' @keywords internal
rint <- function(x) {
  stats::qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}
