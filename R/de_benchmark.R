#' Merge results from differential expression tools
#'
#' Merge results from differential expression tools to benchmark calls. Plots and stuff can be made from this object as well
#'
#' @param de_list a list of de results
#' @param de_label a character vector of the same length as de_list with labels for each method
#' @param oracle a data.frame with the columns target_id, is_de and other optional columns (TBD)
#' @param de_colors a named list of what colors to assign to what tool. The name corresponds to the method. If not named, will assign in alphabetical order.
#' @export
new_de_benchmark <- function(de_list, de_labels, oracle, de_colors = NULL) {
  stopifnot( is(de_list, "list") )
  stopifnot( length(de_list) == length(de_labels) )
  stopifnot( is(de_labels, "character") )

  if (!is.null(de_colors)) {
    # user supplied colors
    stopifnot( length(de_list) == length(de_colors) )
    if (!is.null(names(de_colors))) {
      # ensure that the user supplied names match the actual names
      stopifnot(sort(names(de_colors)) == sort(de_labels))
    } else {
      names(de_colors) <- sort(de_labels)
    }
  } else {
    # since the user didn't supply any colors, give them some colorblind ones

    # thank you: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
    colorblind_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7")

    # if we don't have enough colorblind colors, add some random ones
    if (length(de_list) > length(colorblind_colors)) {
      set.seed(1)
      remaining_colors <- sample(colors(),
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

  all_res <- Reduce(function(x, y) dplyr::inner_join(x, y, by = c("target_id")),
    de_list)

  oracle <- oracle %>%
    dplyr::select(target_id, is_de, log_fc)

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

  all_res <- all_res %>%
    dplyr::inner_join(data.table::data.table(oracle), by = "target_id")

  ret <- list(
    all_data = as_df(all_res),
    m_pval = as_df(m_pval),
    m_qval = as_df(m_qval),
    labels = de_labels,
    color_mapping = de_colors
    )
  class(ret) <- "de_benchmark"
  ret
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
  p <- p + xlab("eFDR")
  p <- p + ylab("TPR")
  p <- p + xlim(0, 1)
  p <- p + ylim(0, 1)
  p <- p + scale_color_manual(values = de_bench$color_mapping)


}

calculate_fdr <- function(de_bench) {
  stopifnot( is(de_bench, "de_benchmark") )

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
  pvals <- dplyr::arrange(pvals, estimate)
  pvals <- dplyr::mutate(pvals, nde = 1:n(), tFDR = cummean(!is_de))

  list(n_true_de = n_true_de, pvals = pvals)
}

#' FDR versus estimated FDR
#'
#' Plot the true fdr versus the estimated fdr
#'
#' @param de_bench a \code{de_benchmark} object resulting from \code{merge_results}
#' @return a \code{ggplot} object
#' @export
fdr_efdr_plot <- function(de_bench) {
  stopifnot( is(de_bench, "de_benchmark") )

  fdr_obj <- calculate_fdr(de_bench)
  pvals <- fdr_obj$pvals

  p <- ggplot(pvals, aes(qval, tFDR))
  p <- p + geom_line(aes(color = method, linetype = method))
  p <- p + xlab("eFDR")
  p <- p + ylab("FDR")
  p <- p + scale_color_manual(values = de_bench$color_mapping)

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
  stopifnot( is(de_bench, "de_benchmark") )

  fdr_obj <- calculate_fdr(de_bench)
  pvals <- fdr_obj$pvals

  plt <- ggplot(pvals, aes(nde, estimate, group = method))

  if (estimate) {
    plt <- plt + geom_line(aes(nde, qval, colour = method), linetype = 3)
  }

  plt <- plt +
    geom_line(aes(nde, tFDR, colour = method, linetype = method),
      size = 0.8, alpha = 0.8) +
    geom_vline(xintercept = fdr_obj$n_true_de, linetype = 3) +
    geom_hline(yintercept = 0.10, linetype = 3) +
    #theme_bw() +
    xlab("Number of features called DE") +
    ylab("FDR") +
    ylim(0, 1)
  plt <- plt + scale_color_manual(values = de_bench$color_mapping)

  plt
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
