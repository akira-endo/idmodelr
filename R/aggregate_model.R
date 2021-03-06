
#' A Function to Apply Multiple Aggregations to Model Output
#'
#' @inherit aggregate_model_internal
#' @param aggregate_to A character vector or list specifying the aggregation operations to perform on the
#' model output. Operations are carried out in the order specified. Implemented options are; disease, demographic,
#' and incidence.
#' @param compartments A character vector or list specifying the unique compartments to aggregate. May either be
#' specified once for all aggregation functions or for each function separately.
#' @param hold_out_var A character vector or list specifying the unique compartments that will not be aggregated. May either be
#' specified once for all aggregation functions or for each function separately. If compartments is set then this argument does
#' not need to be used.
#' @param test, Logical defaults to \code{FALSE}. For testing, returns the processed inputs rather than
#' performing the aggregation.
#' @export
#'
#' @examples
#'
#' df <- data.frame(A1 = 1, B1 = 1, A2 = 1, B2 = 1, A3 = 1, B3 = 1)
#' aggregate_model(df, aggregate_to = "incidence",
#'                 compartments = c("A", "B"), strat = 3,
#'                 summary_var = TRUE)
#'
#'
aggregate_model <- function(df, aggregate_to = NULL, compartments = NULL, strat = NULL,
                            hold_out_var = NULL, id_col = NULL, groups = NULL,
                            new_var = "incidence", total_pop = TRUE, summary_var = FALSE,
                            test = FALSE) {
  if (length(aggregate_to) > 1) {
    if (length(compartments) == 1) {
      compartments <- rep(list(compartments), length(aggregate_to))
    }else if (!(length(compartments) == length(aggregate_to))) {
      stop("The aggregate_to contains ", length(aggregate_to),
           " actions, whilst compartments contains ", length(compartments),
           ". Compartments must either be of a vector of compartments or a list of equal length to aggregate_to.")
    }

    if (length(hold_out_var) == 1) {
      hold_out_var <- rep(list(hold_out_var), length(aggregate_to))
    }else if (!(length(hold_out_var) == length(aggregate_to))) {
      stop("The aggregate_to contains ", length(aggregate_to),
           " actions, whilst hold_out_var contains ", length(hold_out_var),
           ". hold_out_var must either be of a vector of compartments or a list of equal length to aggregate_to.")
    }

    if (length(strat) == 1) {
      strat <- rep(strat, length(aggregate_to))
    }else if (!(length(strat) == length(aggregate_to))) {
      stop("The aggregate_to contains ", length(aggregate_to),
           " actions, whilst strat contains ", length(strat),
           ". strat must either be of a single number or a numeric vector of equal length to aggregate_to.")
    }
  }else{
      compartments <- list(compartments)
      hold_out_var <- list(hold_out_var)
  }

  if (!test) {
    for (i in 1:length(aggregate_to)) {
      df <- aggregate_model_internal(df, aggregate_to = aggregate_to[i],
                                     compartments = compartments[[i]],
                                     strat = strat[i], hold_out_var = hold_out_var[[i]],
                                     new_var = new_var,
                                     id_col = id_col, groups = groups, total_pop = total_pop,
                                     summary_var = summary_var)
    }
  }else{
    df <- c(aggregate_to, compartments, strat,
               hold_out_var,  new_var, id_col, groups)
  }


  return(df)
}
