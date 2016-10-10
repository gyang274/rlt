#------------------------------------------------------------------------------#
#------------------------------- rlt::skm_gdp.r -------------------------------#
#------------------------- author: gyang274@gmail.com -------------------------#
#------------------------------------------------------------------------------#

#--------+---------+---------+---------+---------+---------+---------+---------#
#234567890123456789012345678901234567890123456789012345678901234567890123456789#

#------------------------------------------------------------------------------#
#------------------------------------ main ------------------------------------#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#------------------------------------ test ------------------------------------#
#------------------------------------------------------------------------------#

#' skm_gdp_val0_r: selective skm with greedy propgation on val0 cells.
#' skm_gdp_val0_r: with an input m x n matrix - matrix is init as all rows indentical
#' and sum to 1 followed by mask a selective columns for each row as 0 - m are stores
#' sources, n are zipcode desinations and v the cell values are percentage population
#' lives in the destination, and cells are masked as 0 when s - t within a given dist
#' X miles, want to select m sequentially so that always covers as much as popluation
#' objective is sequential max(sum( min(dist(s, t) - all s) < Xmile ? p : 0 - all t))
#'
#' @param x: m x n matrix: s - t - p
#' @return s: index 1 - m gives highest overall p:
#'   max(sum( min(dist(s, t) - all s) < Xmile ? p : 0 - all t))
skm_gdp_val0_r <- function(x) {

  message("skm_gdp_val0_r: selective skm with greedy propgation on val0 cells ...\n")

  .ptc <- proc.time()

  rownames(x) <- 1:nrow(x)

  s = numeric(nrow(x))

  for ( i in 1L:nrow(x) ) {

    message("skm_gdp_val0_r: optimizing at it <", i, "> ...\n")

    # .ik_ptc <- proc.time()

    if ( ncol(x) > 0 ) {

      #- rowSum as the %_pop outside Xmile for each location
      rs <- rowSums(x)

      #- s[i] wrt. the %_pop outside Xmile for each location
      rk <- order(rs)[1L]

      s[i] <- as.integer(names(rs[rk]))

      # update x by removing column t with 0 (dist within Xmile) to s selected
      col_val0_idx = which(x[rk, , drop=TRUE] == 0)

      if ( length(col_val0_idx) > 0 ) {

        x <- x[-rk, -col_val0_idx, drop = FALSE]

      } else {

        x <- x[-rk,              , drop = FALSE]

      }

      # cat("skm_gdp_val0_r: optimizing at it <", i, "> ... update s to\n", print(s), "\n")

      # cat("skm_gdp_val0_r: optimizing at it <", i, "> ... update x to\n", print(x), "\n")

    } else {

      # cat("skm_gdp_val0_r: optimizing at it <", i, "> ... break ... \n")

      break;

    }

    .ik_ptd <- proc.time() - .ik_ptc

    message("skm_gdp_val0_r: optimizing at it <", i, "> consumes ", .ik_ptd[3], " seconds.\n")

    message("skm_gdp_val0_r: optimizing at it <", i, "> ... done.\n")
  }

  if ( nrow(x) > 0 ) {

    s[(length(s) - nrow(x) + 1):length(s)] = rownames(x)

  }

  # .ptd <- proc.time() - .ptc

  # message("skm_gdp_val0_r: cosumes ", .ptd[3], " seconds.\n")

  # message("skm_gdp_val0_r: selective skm with greedy propgation on val0 cells ... done.\n")

  return(s)

}

#------------------------------------------------------------------------------#
