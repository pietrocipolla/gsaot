#' Plot Optimal Transport inner statistics
#'
#' Plot Optimal Transport based inner statistics for each partition using
#' `ggplot2` package. If provided, it plots also the uncertainty estimates.
#'
#' @inheritParams plot.gsaot_indices
#'
#' @returns A \code{patchwork} object that, if called, will print.
#' @export
#'
#' @examples
#' N <- 1000
#'
#' mx <- c(1, 1, 1)
#' Sigmax <- matrix(data = c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#'
#' x1 <- rnorm(N)
#' x2 <- rnorm(N)
#' x3 <- rnorm(N)
#'
#' x <- cbind(x1, x2, x3)
#' x <- mx + x %*% chol(Sigmax)
#'
#' A <- matrix(data = c(4, -2, 1, 2, 5, -1), nrow = 2, byrow = TRUE)
#' y <- t(A %*% t(x))
#'
#' x <- data.frame(x)
#'
#' M <- 25
#'
#' # Get sensitivity indices
#' sensitivity_indices <- ot_indices(x, y, M)
#' plot_inner_stats(sensitivity_indices)
#'
plot_inner_stats <- function(x,
                             ranking = NULL,
                             wb_all = FALSE, ...) {
  # Get the number of inputs
  K <- length(x$indices)

  # Select only the inputs requested by `ranking`
  # ----------------------------------------------------------------------------
  if (!is.null(ranking)) {
    # Check if ranking is an integer and is less or equal than the total number
    # of inputs
    if (ranking %% 1 == 0 & abs(ranking) <= K) {
      inputs_to_plot <-
        ifelse(rep(sign(ranking), each = abs(ranking)) > 0,
               seq(ranking),
               seq(from = K + ranking + 1, to = K))

      # Get the indices of the inputs to be plotted
      inputs_to_plot <- which(rank(-x$indices) %in% inputs_to_plot)
    } else
      stop("`ranking` should be an integer with absolute value less than the number of inputs")
  } else
    inputs_to_plot <- seq(K)

  # Order the plots from the most important to the least one
  inputs_to_plot <- inputs_to_plot[order(x$indices[inputs_to_plot], decreasing = TRUE)]

  # Get the names of the inputs
  input_names <- names(x$indices)

  # Generate the list of plots
  gg_plots <- list()

  for (k in inputs_to_plot) {
    # Get the partition for the k-th input and its dimension
    partition <- x$partitions[, k]
    M <- max(partition)

    # If the input is numeric, get the mean of each partition as a x-point to be plotted
    if (is.numeric(x$x[, k])) {
      x_means <- array(dim = M)

      # Get the mean for each partition as reference for x
      for (m in seq(M))
        x_means[m] <- mean(x$x[which(partition == m), k])

      # Build the plotting structure
      data_to_plot <- data.frame(
        x = x_means,
        y = x$inner_statistics[[k]][1,],
        Component = x$method
      )

      # Build the plot of the inner statistic
      gg_plots[[k]] <- ggplot2::ggplot(data = data_to_plot,
                                       ggplot2::aes(x = .data[["x"]],
                                                    y = .data[["y"]])) +
        ggplot2::geom_point() +
        ggplot2::geom_line(linetype = 2, linewidth = 0.1) +
        ggplot2::labs(x = input_names[k],
                      y = "Inner statistics")  +
        ggplot2::guides(color = "none")

      # If the bootstrap estimates are present, plot them
      if (x$boot) {
        data_to_plot <- cbind(data_to_plot, x$inner_statistics_ci[[k]][1:M, ])
        colnames(data_to_plot)[4:5] <- c("low_ci", "high_ci")

        gg_plots[[k]] <- gg_plots[[k]] +
          ggplot2::geom_errorbar(data = data_to_plot, ggplot2::aes(ymin = .data[["low_ci"]],
                                                                   ymax = .data[["high_ci"]]),
                                 width = 0.5 / M * (max(x$x[, k]) - min(x$x[, k])))
      }

      # If the advective and diffusive components are present, plot them
      if (exists("adv", where = x) & wb_all) {
        data_to_plot <- data.frame(
          x = x_means,
          y = x$inner_statistics[[k]][2,],
          Component = x$method
        )

        gg_plots[[K + k]] <- ggplot2::ggplot(data = data_to_plot,
                                             ggplot2::aes(x = .data[["x"]],
                                                          y = .data[["y"]])) +
          ggplot2::geom_point() +
          ggplot2::geom_line(linetype = 2, linewidth = 0.1) +
          ggplot2::labs(x = input_names[k],
                        y = "Inner statistics (Advective)")  +
          ggplot2::guides(color = "none")

        if (x$boot) {
          data_to_plot <- cbind(data_to_plot, x$inner_statistics_ci[[k]][(M + 1):(2 * M), ])
          colnames(data_to_plot)[4:5] <- c("low_ci", "high_ci")

          gg_plots[[K + k]] <- gg_plots[[K + k]] +
            ggplot2::geom_errorbar(data = data_to_plot, ggplot2::aes(ymin = .data[["low_ci"]],
                                                                     ymax = .data[["high_ci"]]),
                                   width = 0.5 / M * (max(x$x[, k]) - min(x$x[, k])))
        }

        data_to_plot <- data.frame(
          x = x_means,
          y = x$inner_statistics[[k]][3,],
          Component = x$method
        )

        gg_plots[[2 * K + k]] <-
          ggplot2::ggplot(data = data_to_plot,
                          ggplot2::aes(x = .data[["x"]],
                                       y = .data[["y"]])) +
          ggplot2::geom_point() +
          ggplot2::geom_line(linetype = 2, linewidth = 0.1) +
          ggplot2::labs(x = input_names[k],
                        y = "Inner statistics (Diffusive)")  +
          ggplot2::guides(color = "none")

        if (x$boot) {
          data_to_plot <- cbind(data_to_plot, x$inner_statistics_ci[[k]][(2 * M + 1):(3 * M), ])
          colnames(data_to_plot)[4:5] <- c("low_ci", "high_ci")

          gg_plots[[2 * K + k]] <- gg_plots[[2 * K + k]] +
            ggplot2::geom_errorbar(data = data_to_plot, ggplot2::aes(ymin = .data[["low_ci"]],
                                                                     ymax = .data[["high_ci"]]),
                                   width = 0.5 / M * (max(x$x[, k]) - min(x$x[, k])))
        }
      }
    } else {
      # If the input is non-numeric, use barplots
      if (is.factor(x$x[, k])) {
        x_means <- levels(x$x[, k])
      } else {
        x_means <- unique(x$x[, k])
      }

      data_to_plot <- data.frame(
        x = x_means,
        y = x$inner_statistics[[k]][1,],
        Component = x$method
      )

      gg_plots[[k]] <- ggplot2::ggplot(data = data_to_plot,
                                       ggplot2::aes(x = .data[["x"]],
                                                    y = .data[["y"]])) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(x = input_names[k],
                      y = "Inner statistics")  +
        ggplot2::guides(fill = "none")

      # If the bootstrap estimates are present, plot them
      if (x$boot) {
        data_to_plot <- cbind(data_to_plot, x$inner_statistics_ci[[k]][1:M, ])
        colnames(data_to_plot)[4:5] <- c("low_ci", "high_ci")

        gg_plots[[k]] <- gg_plots[[k]] +
          ggplot2::geom_errorbar(data = data_to_plot, ggplot2::aes(ymin = .data[["low_ci"]],
                                                                   ymax = .data[["high_ci"]]),
                                 width = 0.5 / M * (max(x$x[, k]) - min(x$x[, k])))
      }

      if (exists("adv", where = x) & wb_all) {
        data_to_plot <- data.frame(
          x = x_means,
          y = x$inner_statistics[[k]][2,],
          Component = x$method
        )

        gg_plots[[2 * K + k]] <-
          ggplot2::ggplot(data = data_to_plot,
                          ggplot2::aes(x = .data[["x"]],
                                       y = .data[["y"]])) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::labs(x = input_names[k],
                        y = "Inner statistics (Advective)")  +
          ggplot2::guides(fill = "none")

        if (x$boot) {
          data_to_plot <- cbind(data_to_plot, x$inner_statistics_ci[[k]][(M + 1):(2 * M), ])
          colnames(data_to_plot)[4:5] <- c("low_ci", "high_ci")

          gg_plots[[K + k]] <- gg_plots[[K + k]] +
            ggplot2::geom_errorbar(data = data_to_plot, ggplot2::aes(ymin = .data[["low_ci"]],
                                                                     ymax = .data[["high_ci"]]),
                                   width = 0.5 / M * (max(x$x[, k]) - min(x$x[, k])))
        }

        data_to_plot <- data.frame(
          x = x_means,
          y = x$inner_statistics[[k]][3,],
          Component = x$method
        )

        gg_plots[[2 * K + k]] <-
          ggplot2::ggplot(data = data_to_plot,
                          ggplot2::aes(x = .data[["x"]],
                                       y = .data[["y"]])) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::labs(x = input_names[k],
                        y = "Inner statistics (Diffusive)")  +
          ggplot2::guides(fill = "none")

        if (x$boot) {
          data_to_plot <- cbind(data_to_plot, x$inner_statistics_ci[[k]][(2 * M + 1):(3 * M), ])
          colnames(data_to_plot)[4:5] <- c("low_ci", "high_ci")

          gg_plots[[2 * K + k]] <- gg_plots[[2 * K + k]] +
            ggplot2::geom_errorbar(data = data_to_plot, ggplot2::aes(ymin = .data[["low_ci"]],
                                                                     ymax = .data[["high_ci"]]),
                                   width = 0.5 / M * (max(x$x[, k]) - min(x$x[, k])))
        }
      }
    }
  }

  # Order the elements in decreasing importance
  gg_plots <- gg_plots[inputs_to_plot]

  # Remove the NULL plots before the patchwork
  if (any(sapply(gg_plots, is.null))) {
    null_indices <- which(sapply(gg_plots, is.null))
    gg_plots <- gg_plots[-null_indices]
  }

  # Plot the inner statistics according to the grid
  patchwork::wrap_plots(
    gg_plots,
    nrow = length(inputs_to_plot),
    ncol = ifelse(wb_all, 3, 1),
    byrow = FALSE
  ) +
    patchwork::plot_layout(axis_titles = "collect_y")
}
