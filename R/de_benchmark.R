#' Merge results from differential expression tools
#'
#' Merge results from differential expression tools to benchmark calls. Plots and stuff can be made from this object as well
#'
#' @param de_list a list of de results
#' @param de_label a character vector of the same length as de_list with labels for each method
#' @param oracle a data.frame with the columns target_id, is_de and other optional columns (TBD)
#' @param de_colors a named list of what colors to assign to what tool. The name corresponds to the method. If not named, will assign in alphabetical order. If \code{NULL}, colors will be chosen automatically.
#' @export
new_de_benchmark <- function(de_list, de_labels, oracle, de_colors = NULL,
  join_mode = 'intersect') {
  stopifnot( is(de_list, "list") )
  stopifnot( length(de_list) == length(de_labels) )
  stopifnot( is(de_labels, "character") )

  names(de_list) <- de_labels

  original_oracle <- oracle
  original_de_list <- de_list

  if (!is.null(de_colors)) {
    # user supplied colors
    stopifnot( length(de_list) == length(de_colors) )
    if (!is.null(names(de_colors))) {
      # ensure that the user supplied names match the actual names
      de_colors <- de_colors[de_labels]
      stopifnot(sort(names(de_colors)) == sort(de_labels))
    } else {
      names(de_colors) <- sort(de_labels)
    }
  } else {
    # since the user didn't supply any colors, give them some colorblind ones

    # thank you: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
    colorblind_colors <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7",
      "#0072B2", "#D55E00", "#FF0000")

    # if we don't have enough colorblind colors, add some random ones
    if (length(de_list) > length(colorblind_colors)) {
      set.seed(1)
      remaining_colors <- sample(grDevices::colors(),
        length(de_list) - length(colorblind_colors))
      colorblind_colors <- c(colorblind_colors, remaining_colors)
    }

    de_colors <- colorblind_colors[1:length(de_labels)]
    names(de_colors) <- sort(de_labels)
  }

  # extract relevant columns and rename
  de_list <- lapply(seq_along(de_list),
    function(i) {
      res <- de_list[[i]] %>%
        dplyr::select(target_id, pval, qval) %>%
        data.table::data.table()

      c('pval', 'qval') %>%
        data.table::setnames(res, ., paste0(., '_', de_labels[i]))

      res
    })

  # TODO: consider s/inner_join/full_join because of filtering
  all_res <- NULL
  if (join_mode == 'intersect') {
    all_res <- Reduce(function(x, y) dplyr::inner_join(x, y, by = c("target_id")),
      de_list)
  } else {
    all_res <- Reduce(function(x, y) dplyr::full_join(x, y, by = c("target_id")),
      de_list)
  }

  oracle <- dplyr::select(oracle, target_id, is_de, log_fc, matches('qval'))
  if ('qval' %in% colnames(oracle)) {
    oracle <- dplyr::rename(oracle, oracle_fdr = qval)
  }

  melt_by <- function(data, unit_by) {
    m_unit <- data %>%
      dplyr::select(target_id, starts_with(unit_by)) %>%
      reshape2::melt(id.vars = "target_id", variable.name = "method")
    ret <- data.table::data.table(oracle) %>%
      dplyr::inner_join(data.table::data.table(m_unit), by = "target_id") %>%
      dplyr::rename(estimate = value)
    # data.table::setnames(ret, paste0(unit_by, "_oracle"), "oracle")
    ret
  }
  m_pval <- melt_by(all_res, "pval")
  m_qval <- melt_by(all_res, "qval")

  if (join_mode == 'intersect') {
    all_res <- all_res %>%
      dplyr::inner_join(data.table::data.table(oracle), by = "target_id")
  } else {
    all_res <- dplyr::left_join(data.table::data.table(oracle),
      data.table::data.table(all_res), by = "target_id")
  }


  names(de_list) <- de_labels

  ret <- list(
    all_data = as_df(all_res),
    m_pval = as_df(m_pval),
    m_qval = as_df(m_qval),
    labels = de_labels,
    color_mapping = de_colors,
    original_data = original_de_list,
    oracle = original_oracle
    )
  class(ret) <- "de_benchmark"
  ret
}

#' Filter a benchmark
#'
#' Filter a benchmark
#' @param de_bench a 'de_benchmark' object
#' @param to_remove a character string of labels indicating which results to remove
#' @param keep_colors if \code{TRUE}, keep the original colors
#' @return a \code{de_benchmark} object with \code{to_remove} removed
#' @export
filter_benchmark <- function(de_bench, to_remove, join_mode,
  keep_colors = FALSE) {
  stopifnot( is(de_bench, 'de_benchmark') )

  keep_index <- !names(de_bench$original_data) %in% to_remove
  original_data <- de_bench$original_data[keep_index]
  color_mapping <- NULL
  if (keep_colors) {
    keep_index <- !(names(de_bench$color_mapping) %in% to_remove)
    color_mapping <- de_bench$color_mapping[keep_index]
  }

  new_de_benchmark(original_data, names(original_data), oracle = de_bench$oracle,
    de_colors = color_mapping, join_mode = join_mode)
}

#' rename a benchmark
#'
#' rename a benchmark object
#' @param de_bench a \code{de_benchmark} object
#' @param original_names the names of the original method to be renamed
#' @param new_names the new names to replace the original names
#' @return a new \code{de_benchmark} object
#' @export
rename_benchmark <- function(de_bench, original_names, new_names, join_mode) {
  stopifnot( is(de_bench, 'de_benchmark') )

  original_data <- de_bench$original_data

  original_labels <- de_bench$labels
  new_labels <- original_labels
  for (i in seq_along(original_labels)) {
    j <- which(original_names[i] == original_labels)
    new_labels[j] <- new_names[i]
  }

  color_mapping <- de_bench$color_mapping[original_labels]
  names(color_mapping) <- new_labels
  new_de_benchmark(original_data, new_labels, oracle = de_bench$oracle,
    de_colors = color_mapping, join_mode = join_mode)
}

#' Plot the (estimated) FDR against the TPR
#'
#' Plot the (estimated) FDR against the TPR.
#' Currently can only do the estimated FDR.
#' Will eventually support the true FDR
#'
#' @param de_bench a \code{de_benchmark} object
#' @return a \code{ggplot} object
#' @export
fdr_tpr_plot <- function(de_bench) {
  stopifnot( is(de_bench, "de_benchmark") )

  tmp <- dplyr::mutate(de_bench$m_qval, method = sub("qval_", "", method))

  tmp <- dplyr::group_by(tmp, method)
  tmp <- dplyr::arrange(tmp, estimate)
  tmp <- dplyr::mutate(tmp, tpr = cummean(is_de))

  p <- ggplot(tmp, aes(estimate, tpr, group = method))
  p <- p + geom_line(aes(colour = method))
  p <- p + xlab("estimated FDR")
  p <- p + ylab("TPR")
  p <- p + xlim(0, 1)
  p <- p + ylim(0, 1)
  p <- p + scale_color_manual(values = de_bench$color_mapping)

  p
}

average_bench_fdr <- function(list_bench,
  sim_filter = FALSE,
  jagged_summary = FALSE) {
  stopifnot(all(sapply(list_bench, class) == 'de_benchmark'))

  all_fdr <- lapply(seq_along(list_bench),
    function(i) {
      bench <- list_bench[[i]]
      res <- calculate_fdr(bench, sim_filter = sim_filter)
      res$pvals <- dplyr::mutate(res$pvals, sample = i)
      res$pvals
    })

  all_fdr <- dplyr::bind_rows(all_fdr)
  all_fdr <- dplyr::group_by(all_fdr, method, nde)
  all_fdr <- dplyr::mutate(all_fdr, n_samples = length(nde))
  if (!jagged_summary) {
    all_fdr <- dplyr::filter(all_fdr, n_samples == length(list_bench))
  }

  # TODO: figure out why getting NAs in the sd columns
  all_fdr <- dplyr::summarize(all_fdr,
    sd_qval = sd(qval),
    qval = mean(qval),
    sd_tFDR = sd(tFDR),
    tFDR = mean(tFDR),
    sd_sensitivity = sd(sensitivity),
    sensitivity = mean(sensitivity),
    specificity = mean(specificity),
    true_fdr = mean(true_fdr),
    n_samples = length(p),
    p = mean(p),
    n = mean(n)
    )

  list(pvals = all_fdr)
}

#' Get the sensitivity and specificity at a level
#'
#' Given some levels, calculate the (estimated) sensitivity and specificity.
#' Note: this requires that the oracle has a column 'qval'
#' @param de_bench a \code{de_benchmark} object
#' @param fdr_level a vector of fdr levels (e.g. \code{c(0.01, 0.05, 0.10)})
#' @return a data frame with sensitivity, specificity, precision, accuracy
#' @export
get_sensitivity_specificity <- function(de_bench,
  fdr_level = c(0.01, 0.05, 0.10)) {
  p <- sum(de_bench$oracle$is_de)
  n <- sum(!de_bench$oracle$is_de)

  m_qval <- de_bench$m_qval
  m_qval <- dplyr::group_by(m_qval, method)
  m_qval <- dplyr::filter(m_qval, !is.na(estimate) & !is.na(oracle_fdr))
  res <- lapply(fdr_level,
    function(f_level) {
      m_qval <- dplyr::mutate(m_qval, is_de = oracle_fdr <= f_level)
      dplyr::summarize(m_qval,
        tp = sum(estimate <= f_level & is_de),
        fp = sum(estimate <= f_level & !is_de),
        fn = sum(estimate > f_level & is_de),
        tn = sum(estimate > f_level & !is_de),
        # p = sum(is_de),
        # n = sum(!is_de),

        sensitivity = tp / p,
        specificity = tn / n,
        precision = tp / (tp + fp),
        accuracy = (tp + tn) / (tp + fp + tn + fn),
        true_fdr = fp / (tp + fp),
        fdr_level = f_level
      )
    })

  dplyr::bind_rows(res)
}

#' Get the sensitivity and specificity at a level
#'
#' Given some levels, calculate the (estimated) sensitivity and specificity.
#' Note: this requires that the oracle does NOT have a column 'qval'
#' @param de_bench a \code{de_benchmark} object
#' @param fdr_level a vector of fdr levels (e.g. \code{c(0.01, 0.05, 0.10)})
#' @return a data frame with sensitivity, specificity, precision, accuracy
#' @export
get_sensitivity_specificity_oracle <- function(de_bench,
  fdr_level = c(0.01, 0.05, 0.10)) {
  m_qval <- de_bench$m_qval
  m_qval <- dplyr::mutate(m_qval, method = sub('qval_', '', method))
  m_qval <- dplyr::group_by(m_qval, method)
  m_qval <- dplyr::filter(m_qval, !is.na(estimate))

  p <- sum(de_bench$oracle$is_de)
  n <- sum(!de_bench$oracle$is_de)

  res <- lapply(fdr_level,
    function(f_level) {
      dplyr::summarize(m_qval,
        tp = sum(estimate <= f_level & is_de),
        fp = sum(estimate <= f_level & !is_de),
        fn = sum(estimate > f_level & is_de),
        tn = sum(estimate > f_level & !is_de),
        # p = sum(is_de),
        # n = sum(!is_de),

        sensitivity = tp / p,
        specificity = tn / n,
        precision = tp / (tp + fp),
        accuracy = (tp + tn) / (tp + fp + tn + fn),
        true_fdr = fp / (tp + fp),
        fdr_level = f_level
      )
    })

  dplyr::bind_rows(res)
}

#' @export
average_sensitivity_specificity <- function(de_bench_list,
  fdr_level = c(0.01, 0.05, 0.10), use_oracle = FALSE) {
  stopifnot(is(de_bench_list, 'list'))
  stopifnot(all(sapply(de_bench_list, is, 'de_benchmark')))

  res <- lapply(seq_along(de_bench_list),
    function(i) {
      bench <- de_bench_list[[i]]
      res <- if (use_oracle) {
        get_sensitivity_specificity_oracle(bench, fdr_level)
      } else {
        get_sensitivity_specificity(bench, fdr_level)
      }
      dplyr::mutate(res, sample = i)
    })

  res <- dplyr::bind_rows(res)
  # res <- dplyr::group_by(res, method, fdr_level)
  #
  # dplyr::summarize(res,
  #   sensitivity = mean(sensitivity, na.rm = TRUE),
  #   specificity = mean(specificity, na.rm = TRUE),
  #   precision = mean(precision, na.rm = TRUE),
  #   accuracy = mean(accuracy, na.rm = TRUE)
  #   )
  res
}

#' @export
sensitivity_specificity_plot <- function(de_bench,
  fdr_level = c(0.01, 0.05, 0.10)) {

  s_summary <- NULL
  if (is(de_bench, 'list')) {
    s_summary <- average_sensitivity_specificity(de_bench, fdr_level)
  } else {
    s_summary <- get_sensitivity_specificity(de_bench, fdr_level)
  }

  sensitivity_plot <- ggplot(s_summary,
    aes(factor(fdr_level), sensitivity, fill = method)) +
    geom_boxplot(aes(fill = method)) +
    ylim(0, 1)

  precision_plot <- ggplot(s_summary,
    aes(factor(fdr_level), precision, fill = method)) +
    geom_boxplot(aes(fill = method)) +
    ylim(0, 1)

  specificity_plot <- ggplot(s_summary,
    aes(factor(fdr_level), specificity, fill = method)) +
    geom_boxplot(aes(fill = method)) +
    ylim(0, 1)

  fdr_plot <- ggplot(s_summary,
    aes(factor(fdr_level), true_fdr, fill = method)) +
    geom_boxplot(aes(fill = method)) +
    ylim(0, 1) +
    geom_hline(yintercept = 0.01, linetype = 2) +
    geom_hline(yintercept = 0.10, linetype = 3) +
    geom_hline(yintercept = 0.05, linetype = 4)

  together <- cowplot::plot_grid(sensitivity_plot,
    precision_plot,
    fdr_plot,
    ncol = 1)

  list(
    sensitivity = sensitivity_plot,
    precision = precision_plot,
    specificity = specificity_plot,
    fdr = fdr_plot,
    together = together
    )
}

calculate_fdr <- function(de_bench, sim_filter = FALSE) {
  stopifnot( is(de_bench, "de_benchmark") )

  p <- NULL
  n <- NULL
  if (sim_filter) {
    de_low_expression <- dplyr::filter(de_bench$oracle, !sim_filt & is_de)
    # de_low_expression <- dplyr::filter(de_bench$oracle, allZero & is_de)
    message(paste0('filtering out: ', nrow(de_low_expression)))
    # de_bench$all_data <- dplyr::anti_join(de_bench$all_data, de_low_expression,
    #   by = 'target_id')
    # de_bench$m_pval <- dplyr::anti_join(de_bench$m_pval, de_low_expression,
    #   by = 'target_id')
    # de_bench$m_qval <- dplyr::anti_join(de_bench$m_qval, de_low_expression,
    #   by = 'target_id')

    de_bench$all_data <- dplyr::mutate(de_bench$all_data,
      is_de = ifelse(target_id %in% de_low_expression$target_id, FALSE, is_de))
    de_bench$m_pval <- dplyr::mutate(de_bench$m_pval,
      is_de = ifelse(target_id %in% de_low_expression$target_id, FALSE, is_de))
    de_bench$m_qval <- dplyr::mutate(de_bench$m_qval,
      is_de = ifelse(target_id %in% de_low_expression$target_id, FALSE, is_de))

    p <- sum(de_bench$all_data$is_de)
    n <- sum(!de_bench$all_data$is_de)

  } else {
    p <- sum(de_bench$oracle$is_de)
    n <- sum(!de_bench$oracle$is_de)
  }

  n_true_de <- sum(de_bench$all_data$is_de)
  message('Intersection of targets: ', nrow(de_bench$all_data))
  message('Number of truly DE: ', n_true_de)

  pvals <- dplyr::mutate(de_bench$m_pval, method = sub("pval_", "", method))

  qvals <- dplyr::select(de_bench$m_qval, target_id, method, estimate)
  qvals <- dplyr::mutate(qvals, method = sub("qval_", "", method))
  qvals <- dplyr::rename(qvals, qval = estimate)

  pvals <- dplyr::inner_join(pvals, qvals, by = c("target_id", "method"))

  pvals <- dplyr::mutate(pvals, method = sub("pval_", "", method))
  pvals <- dplyr::group_by(pvals, method)
  # sort them by the qval rather than the pval.
  # some programs (DESeq2) oddly report NA when the pval is small
  # pvals <- dplyr::arrange(pvals, estimate)
  pvals <- dplyr::arrange(pvals, qval)
  pvals <- dplyr::mutate(pvals, nde = 1:n(), tFDR = cummean(!is_de))
  # get the number of qval <=
  # pvals <- dplyr::mutate(pvals, tTPR = cumsum(is_de) / nde)
  pvals <- dplyr::mutate(pvals, tp = cumsum(is_de), fp = cumsum(!is_de))
    # p = sum(is_de), n = sum(!is_de))
  pvals <- dplyr::mutate(pvals, fn = p - tp, tn = n - fp)
  pvals <- dplyr::mutate(pvals,
    sensitivity = tp / p,
    specificity = tn / n,
    precision = tp / (tp + fp),
    accuracy = (tp + tn) / (tp + fp + tn + fn),
    true_fdr = fp / (tp + fp),
    p = p,
    n = n
    )
  message('positives: ', p, ' negatives: ', n)

  list(n_true_de = n_true_de, pvals = pvals)
}

#' @export
get_fdr <- function(de_bench,
  sim_filter = FALSE) {
  fdr_obj <- NULL
  if (is(de_bench, 'list')) {
    stopifnot( all(sapply(de_bench, class) == 'de_benchmark') )
    fdr_obj <- average_bench_fdr(de_bench, sim_filter = sim_filter)
  } else {
    stopifnot( is(de_bench, "de_benchmark") )
    fdr_obj <- calculate_fdr(de_bench)
  }

  fdr_obj
}

#' @export
tpr_fdr_plot <- function(de_bench) {
  fdr_obj <- get_fdr(de_bench)
  pvals <- fdr_obj$pvals

  p <- ggplot(pvals, aes(qval, tTPR, group = method))

  # TODO: eventually at the estimated
  # if (estimate) {
  #   p <- p + geom_line(aes(nde, qval, colour = method), linetype = 3)
  # }

  p <- p + geom_line(aes(colour = method, linetype = method),
    size = 0.8, alpha = 0.8)

  if (is(de_bench, 'list')) {
    p <- p + ylab("TPR")
    p <- p + xlab(paste0("FDR (average, n = ", length(de_bench), ")"))
    p <- p + scale_color_manual(values = de_bench[[1]]$color_mapping)
  } else {
    p <- p + geom_vline(xintercept = fdr_obj$n_true_de, linetype = 3)
    p <- p + xlab("FDR")
    p <- p + ylab("TPR")
    p <- p + scale_color_manual(values = de_bench$color_mapping)
  }
  p <- p + geom_hline(yintercept = 0.10, linetype = 3)
  p <- p + xlim(0, 1)
  p <- p + ylim(0, 1)

  p
}

#' true fdr versus power
#'
#' plot the true fdr versus power
#' @param list_bench a list of \code{de_benchmark}s
#' @return a ggplot2 object
#' @export
fdr_power_plot <- function(list_bench, sim_filter = FALSE) {
  suppressMessages(mean_fdr <- get_fdr(list_bench, sim_filter = sim_filter)$pvals)

  p <- ggplot(mean_fdr, aes(true_fdr, sensitivity))
  p <- p + geom_path(aes(color = method, linetype = method),
    size = 0.8, alpha = 0.8)
  p <- p + xlab(paste0('FDR (average, n = ', length(list_bench), ')'))
  p <- p + ylab('power (sensitivity)')
  p <- p + scale_color_manual(values = list_bench[[1]]$color_mapping)

  p
}

#' fdr versus power
#'
#' true fdr versus power. also point out where the estimated fdr is.
#' also draws isolines specifying the number of things in the ranking
#'
#' @param mean_fdr an object from \code{get_fdr}
#' @param fdr_levels a vector of levels
#' @param start where to start the isolines
#' @param jump how big each step is between the isolines
#' @param rank_fdr if \code{TRUE}, legend based off the fdr
#' @param sim_filter if \code{TRUE}, remove things that didn't pass the simulation filter from the truth
#' @param ignore_estimated_fdr if not \code{NULL}, a character vector of method names to ignore when plotting estimated fdr point
#' @return a \code{ggplot2} object
#' @export
fdr_efdr_power_plot <- function(
  mean_fdr,
  fdr_levels = c(0.01, 0.05, 0.10),
  start = 50,
  jump = 50,
  rank_fdr = NULL,
  method_colors = NULL,
  fdr_level_position = -0.001,
  ignore_estimated_fdr = NULL,
  isolines = TRUE
  ) {
  # suppressMessages(mean_fdr <- get_fdr(list_bench)$pvals)
  mean_fdr <- dplyr::group_by(mean_fdr, method)
  which_levels <- lapply(fdr_levels,
    function(current_level) {
      res <- dplyr::filter(mean_fdr, qval <= current_level)
      res <- dplyr::filter(res, nde == max(nde))
      res <- dplyr::mutate(res, fdr_level = current_level)
      res
    })
  which_levels <- dplyr::bind_rows(which_levels)
  which_levels <- dplyr::mutate(which_levels,
    sensitivity_low = sensitivity - 2 * sd_sensitivity,
    sensitivity_high = sensitivity + 2 * sd_sensitivity)
  if (!is.null(ignore_estimated_fdr)) {
    which_levels <- dplyr::anti_join(which_levels,
      data.frame(method = ignore_estimated_fdr, stringsAsFactors = FALSE))
  }
  max_nde <- max(mean_fdr$nde)
  jump_fdr <- dplyr::inner_join(mean_fdr,
    data.frame(nde = seq(start, max_nde, by = jump)),
    by = "nde")

  # get the ranking so that we can rename the legend
  if (!is.null(rank_fdr)) {
    rank_fdr_df <- dplyr::filter(mean_fdr, true_fdr <= rank_fdr)
    rank_fdr_df <- dplyr::filter(rank_fdr_df, max(nde) == nde)
    rank_fdr_df <- dplyr::ungroup(rank_fdr_df)
    rank_fdr_df <- dplyr::arrange(rank_fdr_df, desc(sensitivity))

    if (nrow(rank_fdr_df) != length(unique(mean_fdr$method))) {
      remaining_methods <- setdiff(unique(mean_fdr$method),
        unique(rank_fdr_df$method))
      rank_fdr_df <- dplyr::bind_rows(rank_fdr_df,
        data.frame(method = sort(remaining_methods)))
      warning(paste0('Not all the methods reach ', rank_fdr,
        ". The ones that don't reach are simply getting sorted alphabetically."))
    }

    mean_fdr <- dplyr::ungroup(mean_fdr)
    mean_fdr <- dplyr::mutate(mean_fdr,
      method = factor(method, levels = rank_fdr_df$method, ordered = TRUE))
    mean_fdr <- dplyr::group_by(mean_fdr, method)
  }

  p <- ggplot(mean_fdr, aes(true_fdr, sensitivity))

  # draw isolines method 2
  if (isolines) {
    total_max <- mean_fdr$p[1] + mean_fdr$n[1]
    fdr_lines <- data.frame(nde = seq(start, total_max, by = jump), p = mean_fdr$p[1],
      n = mean_fdr$n[1])
    fdr_lines <- dplyr::mutate(fdr_lines,
      top_x = ifelse( (nde / p) <= 1, 0, 1 - (p / nde)),
      top_y = ifelse( (nde / p) <= 1, (nde / p), 1),
      bottom_x = 1,
      bottom_y = 0)
    fdr_lines <- dplyr::distinct(
      dplyr::select(fdr_lines, nde, top_x, top_y, bottom_x, bottom_y))

    p <- p + geom_segment(
      aes(x = top_x, y = top_y, xend = bottom_x, yend = bottom_y),
      data = fdr_lines, color = 'gray', alpha = 0.8)
    # p + xlim(0, 1.2) + ylim(0, 1.2)

    # draw the isolines
    fdr_lines <- dplyr::mutate(jump_fdr, n_intercept = nde / p)
    fdr_lines <- dplyr::ungroup(fdr_lines)
    fdr_lines <- dplyr::distinct(dplyr::select(fdr_lines, nde, n_intercept))
    # p <- p + geom_abline(aes(intercept = n_intercept, slope = -n_intercept),
    #   data = fdr_lines, color = 'gray', alpha = 0.8)
    # label the isolines
    fdr_lines <- dplyr::mutate(fdr_lines, x = 0)
    p <- p + geom_text(aes(x, n_intercept, label = nde), size = 8, color = 'black',
      alpha = 0.6, data = fdr_lines)

  }

  # draw the fdr path
  p <- p + geom_path(aes(color = method),
    size = 1.2, alpha = 0.8)

  # draw the dashed line to the end
  max_nde <- dplyr::filter(mean_fdr, nde == max(nde))
  max_nde <- dplyr::mutate(max_nde, max_fdr = n / (n + p))
  p <- p + geom_segment(
    aes(x = true_fdr, y = sensitivity, xend = max_fdr, yend = sensitivity,
      color = method), data = max_nde, linetype = 2)
  p <- p + geom_segment(
    aes(x = max_fdr, y = sensitivity, xend = max_fdr, yend = 1),
    data = max_nde, linetype = 2
    )

  # put the rank labels on the line
  # p <- p + geom_text(aes(true_fdr, sensitivity,
  #   color = method, label = nde), data = jump_fdr, size = 8)

  # put points of the estimated fdr versus
  p <- p + geom_point(aes(true_fdr, sensitivity,
    group = method,
    color = method,
    shape = factor(fdr_level)),
    data = which_levels,
    size = 7,
    show.legend = FALSE
    )

  p <- p + geom_errorbar(aes(true_fdr,
    ymax = sensitivity_high,
    ymin = sensitivity_low,
    color = method),
    data = which_levels,
    width = 0.005,
    show.legend = FALSE)

  # put vertical lines on useful fdrs
  fdr_levels_lines <- data.frame(x = c(), y = c(), xend = c(), yend = c())
  for (current_level in fdr_levels) {
    fdr_levels_lines <- dplyr::bind_rows(fdr_levels_lines,
      data.frame(x = current_level, y = 0, xend = current_level,
      yend = 1))
  }
  p <- p + geom_segment(aes(x = x, y = y, yend = yend, xend = xend),
    data = fdr_levels_lines, linetype = 3)

  p <- p + xlab('false discovery rate')
  p <- p + ylab('sensitivity')
  p <- p + geom_point(aes(x, y, shape = factor(x)),
    data = data.frame(x = fdr_levels, y = rep(fdr_level_position,
      length(fdr_levels))),
    # shape = (16, 17, 15):(15 + length(fdr_levels)),
    size = 7, show.legend = FALSE)

  if (!is.null(method_colors)) {
    p <- p + scale_color_manual(values = method_colors)
  }

  p
}

#' @export
fdr_specificity_plot <- function(list_bench) {
  suppressMessages(mean_fdr <- get_fdr(list_bench)$pvals)

  p <- ggplot(mean_fdr, aes(true_fdr, specificity))
  p <- p + geom_path(aes(color = method, linetype = method),
    size = 0.8, alpha = 0.8)
  p <- p + xlab(paste0('FDR (average, n = ', length(list_bench), ')'))
  p <- p + ylab('specificity')
  p <- p + scale_color_manual(values = list_bench[[1]]$color_mapping)

  p
}

#' FDR versus estimated FDR
#'
#' Plot the true fdr versus the estimated fdr
#'
#' @param de_bench a \code{de_benchmark} object resulting from \code{merge_results}
#' @return a \code{ggplot} object
#' @export
fdr_efdr_plot <- function(de_bench) {
  fdr_obj <- get_fdr(de_bench)

  p <- ggplot(fdr_obj$pvals, aes(qval, tFDR))
  p <- p + geom_path(aes(color = method, linetype = method))
  if (is(de_bench, 'list')) {
    p <- p + xlab(paste0('estimated FDR (average, n = ', length(de_bench), ')'))
    p <- p + ylab('FDR')
    p <- p + scale_color_manual(values = de_bench[[1]]$color_mapping)
  } else {
    p <- p + xlab("estimated FDR")
    p <- p + ylab("FDR")
    p <- p + scale_color_manual(values = de_bench$color_mapping)
  }

  p
}

#' @export
boxplot_prep <- function(list_bench, fdr_level = c(0.01, 0.05, 0.10)) {
  stopifnot( is(list_bench, "list") )

  all_summaries <- lapply(seq_along(list_bench),
    function(i) {
      bench <- list_bench[[i]]
      ss <- get_sensitivity_specificity_oracle(bench, fdr_level)
      dplyr::mutate(ss, sample = i)
    })
  all_summaries <- dplyr::bind_rows(all_summaries)

  all_summaries
}

#' True fdr versus estimated fdr
#'
#' Plot the true fdr versus the estimated fdr.
#' Note, this requires that \code{boxplot_prep} be run ahead of time
#'
#' @param bp a \code{boxplot_prep} object
#' @param facet_fdr if \code{TRUE}, use \code{facet_wrap} across the fdr levels
#' @return a \code{ggplot2} object
#' @export
fdr_efdr_boxplot <- function(bp, facet_fdr = TRUE) {
  stopifnot( is(bp, 'data.frame') )

  p <- ggplot(bp, aes(method, true_fdr))
  p <- p + geom_boxplot(aes(color = method))
  p <- p + geom_hline(aes(yintercept = fdr_level), linetype = 2,
    color = 'black')
  p <- p + ylab('true FDR')
  if (facet_fdr) {
    p <- p + facet_wrap(~fdr_level)
  }

  p
}

#' power at a specific estimated fdr
#'
#' Plot the estimated fdr versus power
#' Note, this requires that \code{boxplot_prep} be run ahead of time
#'
#' @param bp a \code{boxplot_prep} object
#' @param facet_fdr if \code{TRUE}, use \code{facet_wrap} across the fdr levels
#' @return a \code{ggplot2} object
#' @export
power_efdr_boxplot <- function(bp, facet_fdr = TRUE) {
  stopifnot( is(bp, 'data.frame') )

  p <- ggplot(bp, aes(method, sensitivity))
  p <- p + geom_boxplot(aes(color = method))
  p <- p + ylab('power (sensitivity)')
  if (facet_fdr) {
    p <- p + facet_wrap(~fdr_level)
  }

  p
}

#' FDR versus ranked list
#'
#' Plot the true FDR versus the ranking induced by p-value
#'
#' @param de_bench a \code{de_benchmark} object resulting from \code{merge_results}
#' @param estimate if \code{TRUE}, then plot the estimated FDR with dotted lines
#' @return a \code{ggplot} object
#' @export
fdr_nde_plot <- function(de_bench, estimate = TRUE) {

  fdr_obj <- get_fdr(de_bench)
  pvals <- fdr_obj$pvals

  p <- ggplot(pvals, aes(nde, tFDR, group = method))

  # TODO: eventually at the estimated
  # if (estimate) {
  #   p <- p + geom_line(aes(nde, qval, colour = method), linetype = 3)
  # }

  p <- p + geom_line(aes(colour = method, linetype = method),
    size = 0.8, alpha = 0.8)

  if (is(de_bench, 'list')) {
    p <- p + xlab("Number of features called DE")
    p <- p + ylab(paste0("FDR (average, n = ", length(de_bench), ")"))
    p <- p + scale_color_manual(values = de_bench[[1]]$color_mapping)
  } else {
    p <- p + geom_vline(xintercept = fdr_obj$n_true_de, linetype = 3)
    p <- p + xlab("Number of features called DE")
    p <- p + ylab("FDR")
    p <- p + scale_color_manual(values = de_bench$color_mapping)
  }
  p <- p + geom_hline(yintercept = 0.10, linetype = 3)
  p <- p + ylim(0, 1)

  p
}

#' @export
pval_distribution <- function(de_bench) {
  stopifnot( is(de_bench, "de_benchmark") )

  de_bench$m_pval %>%
    ggplot(aes(estimate, y = ..density..)) + # nolint
    geom_histogram(binwidth = 0.02) +
    facet_wrap( ~ method )
}

#' @export
de_rank_scatter <- function(de_bench, cutoff = 5000) {
  stopifnot( is(de_bench, "de_benchmark") )
  rank_diff <- function(x) {
    x$true_rank <- rank(-abs(x$log_fc))
    x$est_rank <- rank(x$estimate)
    tmp <- head(arrange(x, est_rank), cutoff)
    mutate(tmp, relative_rank = rank(true_rank))
  }

  grped <- as.data.frame(de_bench$m_qval) %>%
    group_by(method)

  tmp <- grped %>%
    filter(method == 'qval_DESeq2') %>%
    ungroup()

  do(grped, rank_diff(.))
}
