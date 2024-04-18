bootstats <- function(b,
                      conf,
                      type) {
  # Get the number of estimated quantities
  p <- length(b$t0)

  # Define the labels of the output
  lab <- c("original", "bias", "std.error", "low.ci", "high.ci")

  # Initialize the output structure
  tmp <- as.data.frame(matrix(nrow = p,
                              ncol = length(lab),
                              dimnames = list(NULL, lab)))

  # For each estimated quantity evaluate the
  for (i in 1:p) {
    # Original estimation, bias, standard deviation
    tmp[i, "original"] <- b$t0[i]
    tmp[i, "bias"] <- mean(b$t[, i]) - b$t0[i]
    tmp[i, "std.error"] <- stats::sd(b$t[, i])

    # Confidence interval
    if (type == "norm") {
      ci <- boot::boot.ci(b, index = i, type = "norm", conf = conf)

      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$norm[2]
        tmp[i, "high.ci"] <- ci$norm[3]
      }

    } else if (type == "basic") {
      ci <- boot::boot.ci(b, index = i, type = "basic", conf = conf)
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$basic[4]
        tmp[i, "high.ci"] <- ci$basic[5]
      }

    } else if (type == "percent") {
      ci <- boot::boot.ci(b, index = i, type = "perc", conf = conf)

      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$percent[4]
        tmp[i, "high.ci"] <- ci$percent[5]
      }

    } else if (type == "bca") {
      ci <- boot::boot.ci(b, index = i, conf = conf)

      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$bca[4]
        tmp[i, "high.ci"] <- ci$bca[5]
      }
    }
  }

  return(tmp)
}
