#' Plot Optimal Transport separation measures
#'
#' Plot Optimal Transport based separation measures for each partition using
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
#' plot_separations(sensitivity_indices)
#'
plot_separations <- function(x,
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

  # Reassign the number of plots to do
  K <- length(inputs_to_plot)

  # Order the plots from the most important to the least one
  inputs_to_plot <- inputs_to_plot[order(x$indices[inputs_to_plot], decreasing = TRUE)]

  # Get the names of the inputs
  input_names <- names(x$indices)

  # Generate the list of plots
  gg_plots <- list()

  for (input in seq_along(inputs_to_plot)) {
    k <- inputs_to_plot[input]

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
        y = x$separation_measures[[k]][1,],
        Component = x$method
      )

      # Build the plot of the inner statistic
      gg_plots[[input]] <- ggplot2::ggplot(data = data_to_plot,
                                       ggplot2::aes(x = .data[["x"]],
                                                    y = .data[["y"]])) +
        ggplot2::geom_point() +
        ggplot2::geom_line(linetype = 2, linewidth = 0.1) +
        ggplot2::labs(x = input_names[k],
                      y = "Separation Measures")  +
        ggplot2::guides(color = "none")

      # If the bootstrap estimates are present, plot them
      if (x$boot) {
        data_to_plot <- cbind(data_to_plot, x$separation_measures_ci[[k]][1:M, ])
        colnames(data_to_plot)[4:5] <- c("low_ci", "high_ci")

        gg_plots[[input]] <- gg_plots[[input]] +
          ggplot2::geom_errorbar(data = data_to_plot, ggplot2::aes(ymin = .data[["low_ci"]],
                                                                   ymax = .data[["high_ci"]]),
                                 width = 0.2 / M * (max(x$x[, k]) - min(x$x[, k])))
      }

      # If the advective and diffusive components are present, plot them
      if (exists("adv", where = x) & wb_all) {
        data_to_plot <- data.frame(
          x = x_means,
          y = x$separation_measures[[k]][2,],
          Component = x$method
        )

        gg_plots[[K + input]] <- ggplot2::ggplot(data = data_to_plot,
                                             ggplot2::aes(x = .data[["x"]],
                                                          y = .data[["y"]])) +
          ggplot2::geom_point() +
          ggplot2::geom_line(linetype = 2, linewidth = 0.1) +
          ggplot2::labs(x = input_names[k],
                        y = "Separation Measures (Advective)")  +
          ggplot2::guides(color = "none")

        if (x$boot) {
          data_to_plot <- cbind(data_to_plot, x$separation_measures_ci[[k]][(M + 1):(2 * M), ])
          colnames(data_to_plot)[4:5] <- c("low_ci", "high_ci")

          gg_plots[[K + input]] <- gg_plots[[K + input]] +
            ggplot2::geom_errorbar(data = data_to_plot, ggplot2::aes(ymin = .data[["low_ci"]],
                                                                     ymax = .data[["high_ci"]]),
                                   width = 0.2 / M * (max(x$x[, k]) - min(x$x[, k])))
        }

        data_to_plot <- data.frame(
          x = x_means,
          y = x$separation_measures[[k]][3,],
          Component = x$method
        )

        gg_plots[[2 * K + input]] <-
          ggplot2::ggplot(data = data_to_plot,
                          ggplot2::aes(x = .data[["x"]],
                                       y = .data[["y"]])) +
          ggplot2::geom_point() +
          ggplot2::geom_line(linetype = 2, linewidth = 0.1) +
          ggplot2::labs(x = input_names[k],
                        y = "Separation Measures (Diffusive)")  +
          ggplot2::guides(color = "none")

        if (x$boot) {
          data_to_plot <- cbind(data_to_plot, x$separation_measures_ci[[k]][(2 * M + 1):(3 * M), ])
          colnames(data_to_plot)[4:5] <- c("low_ci", "high_ci")

          gg_plots[[2 * K + input]] <- gg_plots[[2 * K + input]] +
            ggplot2::geom_errorbar(data = data_to_plot, ggplot2::aes(ymin = .data[["low_ci"]],
                                                                     ymax = .data[["high_ci"]]),
                                   width = 0.2 / M * (max(x$x[, k]) - min(x$x[, k])))
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
        y = x$separation_measures[[k]][1,],
        Component = x$method
      )

      gg_plots[[input]] <- ggplot2::ggplot(data = data_to_plot,
                                       ggplot2::aes(x = .data[["x"]],
                                                    y = .data[["y"]])) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(x = input_names[k],
                      y = "Separation Measures")  +
        ggplot2::guides(fill = "none")

      # If the bootstrap estimates are present, plot them
      if (x$boot) {
        data_to_plot <- cbind(data_to_plot, x$separation_measures_ci[[k]][1:M, ])
        colnames(data_to_plot)[4:5] <- c("low_ci", "high_ci")

        gg_plots[[input]] <- gg_plots[[input]] +
          ggplot2::geom_errorbar(data = data_to_plot, ggplot2::aes(ymin = .data[["low_ci"]],
                                                                   ymax = .data[["high_ci"]]),
                                 width = 0.2 / M * (length(unique(x$x[, k]))))
      }

      if (exists("adv", where = x) & wb_all) {
        data_to_plot <- data.frame(
          x = x_means,
          y = x$separation_measures[[k]][2,],
          Component = x$method
        )

        gg_plots[[K + input]] <-
          ggplot2::ggplot(data = data_to_plot,
                          ggplot2::aes(x = .data[["x"]],
                                       y = .data[["y"]])) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::labs(x = input_names[k],
                        y = "Separation Measures (Advective)")  +
          ggplot2::guides(fill = "none")

        if (x$boot) {
          data_to_plot <- cbind(data_to_plot, x$separation_measures_ci[[k]][(M + 1):(2 * M), ])
          colnames(data_to_plot)[4:5] <- c("low_ci", "high_ci")

          gg_plots[[K + input]] <- gg_plots[[K + input]] +
            ggplot2::geom_errorbar(data = data_to_plot, ggplot2::aes(ymin = .data[["low_ci"]],
                                                                     ymax = .data[["high_ci"]]),
                                   width = 0.2 / M * (length(unique(x$x[, k]))))
        }

        data_to_plot <- data.frame(
          x = x_means,
          y = x$separation_measures[[k]][3,],
          Component = x$method
        )

        gg_plots[[2 * K + input]] <-
          ggplot2::ggplot(data = data_to_plot,
                          ggplot2::aes(x = .data[["x"]],
                                       y = .data[["y"]])) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::labs(x = input_names[k],
                        y = "Separation Measures (Diffusive)")  +
          ggplot2::guides(fill = "none")

        if (x$boot) {
          data_to_plot <- cbind(data_to_plot, x$separation_measures_ci[[k]][(2 * M + 1):(3 * M), ])
          colnames(data_to_plot)[4:5] <- c("low_ci", "high_ci")

          gg_plots[[2 * K + input]] <- gg_plots[[2 * K + input]] +
            ggplot2::geom_errorbar(data = data_to_plot, ggplot2::aes(ymin = .data[["low_ci"]],
                                                                     ymax = .data[["high_ci"]]),
                                   width = 0.2 / M * (length(unique(x$x[, k]))))
        }
      }
    }
  }

  # # Order the elements in decreasing importance
  # gg_plots <- gg_plots[inputs_to_plot]

  # Remove the NULL plots before the patchwork
  if (any(sapply(gg_plots, is.null))) {
    null_indices <- which(sapply(gg_plots, is.null))
    gg_plots <- gg_plots[-null_indices]
  }

  # Plot the Separation Measures according to the grid
  patchwork::wrap_plots(
    gg_plots,
    nrow = length(inputs_to_plot),
    ncol = ifelse(wb_all, 3, 1),
    byrow = FALSE
  ) +
    patchwork::plot_layout(axis_titles = "collect_y")
}

#' Compare sensitivity indices across methods
#'
#' This function takes a list of `gsaot_indices` objects
#' and generates a bar plot comparing the sensitivity indices across different methods.
#'
#' @param x_list A list of S3 objects of class `"gsaot_indices"`, each representing
#'   sensitivity analysis results for a different solver.
#' @param wb_all (default `FALSE`) Logical that defines whether or not to plot
#'   the Advective and Diffusive components of the Wasserstein-Bures indices.
#'
#' @return A `ggplot` object representing the bar plot of sensitivity indices
#'   grouped by input and colored by method.
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
#' # Calculate sensitivity indices
#' ind_wb <- ot_indices_wb(x, y, M)
#' ind_sink <- ot_indices(x, y, M)
#'
#' plot_comparison(list(ind_wb, ind_sink))
#'
plot_comparison <- function(x_list,
                            wb_all = FALSE) {
  all_data <- list()

  # Plot all the inputs
  K <- nrow(x_list[[1]]$indices)
  inputs_to_plot <- seq(K)

  for (i in seq_along(x_list)) {
    obj <- x_list[[i]]

    # Read the method name for plotting
    method_name <- obj$method

    # Format method name with epsilon if method is "sinkhorn" or "sinkhorn_log"
    if (method_name %in% c("sinkhorn", "sinkhorn_log")) {
      epsilon_val <- format(obj$solver_optns$epsilon, scientific = FALSE)
      method_name <- paste0(method_name, " (eps=", epsilon_val, ")")
    }

    if (method_name == "transport") {
      method_name <- if (!is.null(obj$solver_optns$method)) {
        obj$solver_optns$method
      } else {
        "networksimplex"
      }
    }

    if (!exists("adv", where = obj) | !wb_all) {
      inputs <- names(obj$indices[order(obj$indices, decreasing = TRUE)])[inputs_to_plot]
      indices <- unname(obj$indices[order(obj$indices, decreasing = TRUE)])[inputs_to_plot]

      df <- data.frame(
        Method = method_name,
        Inputs = factor(inputs, levels = unique(inputs)),
        Indices = indices
      )
    } else {
      # Add advective and diffusive parts if required
      ordered <- order(obj$indices, decreasing = TRUE)
      inputs <- names(obj$indices[ordered])[inputs_to_plot]

      df <- data.frame(
        Method = rep(c("wass-bures", "advective", "diffusive"), each = length(inputs_to_plot)),
        Inputs = factor(rep(inputs, times = 3), levels = unique(inputs)),
        Indices = c(
          obj$indices[ordered][inputs_to_plot],
          obj$adv[ordered][inputs_to_plot],
          obj$diff[ordered][inputs_to_plot]
        )
      )
    }

    all_data[[i]] <- df
  }

  # Aggregate all the data in a single data.frame
  final_data <- do.call(rbind, all_data)

  # Define shared dodge object
  dodge <- ggplot2::position_dodge2(width = 0.6, padding = 0.2)

  # Plot the comparison
  p <- ggplot2::ggplot(final_data,
                       ggplot2::aes(x = .data[["Inputs"]],
                                    y = .data[["Indices"]],
                                    fill = .data[["Method"]],
                                    group = .data[["Method"]])) +
    ggplot2::geom_bar(stat = "identity",
                      position = dodge,
                      width = 0.6) +
    ggplot2::labs(
      x = "Inputs",
      y = "Indices"
    )

  return(p)
}


