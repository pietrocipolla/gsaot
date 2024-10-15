# Function to extract bootstrap estimates from boot object
# The function comes from https://github.com/arnaldpuy/sensobol/blob/master

bootstats <- function(b,
                      conf,
                      type) {
  # Get the number of estimated quantities
  p <- length(b$t0)

  # Define the labels of the output
  lab <- c("original", "low.ci", "high.ci")

  # Initialize the output structure
  tmp <- as.data.frame(matrix(nrow = p,
                              ncol = length(lab),
                              dimnames = list(NULL, lab)))

  # For each estimated quantity evaluate the
  for (i in 1:p) {
    # Central estimate
    bias <- mean(b$t[, i]) - b$t0[i]
    tmp[i, "original"] <- b$t0[i] - bias

    # Confidence interval
    if (type == "norm") {
      ci <- boot::boot.ci(b, index = i, type = "norm", conf = conf)

      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$normal[2]
        tmp[i, "high.ci"] <- ci$normal[3]
      }

    } else if (type == "basic") {
      ci <- boot::boot.ci(b, index = i, type = "basic", conf = conf)
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$basic[4]
        tmp[i, "high.ci"] <- ci$basic[5]
      }

    } else if (type == "perc") {
      tmp[i, "original"] <- b$t0[i]
      ci <- boot::boot.ci(b, index = i, type = "perc", conf = conf)

      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$percent[4]
        tmp[i, "high.ci"] <- ci$percent[5]
      }

    } else if (type == "bca") {
      ci <- boot::boot.ci(b, index = i, type = "bca", conf = conf)

      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$bca[4]
        tmp[i, "high.ci"] <- ci$bca[5]
      }
    }
  }

  return(tmp)
}
