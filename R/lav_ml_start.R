# lav_ml_start.R: Level-specific starting values for multilevel CFA
#
# Addresses the problem where lavaan's default start = "default" uses the
# pooled (total) covariance for all level-specific parameters, causing
# between-level uniquenesses and loadings to be on the wrong numerical scale
# when the intraclass correlation coefficient (ICC) is low.
#
# Method:
#   - Computes within-level (sigma_w) and between-level (sigma_b) covariances
#     directly from the clustered data using pooled level-specific moments.
#   - Runs fabin (single-factor approximation) on each level's covariance
#     separately, using the correct standardization for each factor.
#   - When sigma_b is not usable (indefinite, negative diagonals), falls back
#     to within-informed loadings matched by indicator set.
#   - Derives factor variance starts via conditional least-squares.
#   - Sets uniquenesses to the residual variance (total - implied factor),
#     floored at 1% of total variance.
#   - Uses level-specific means for intercepts.
#
# Returns a numeric vector suitable for lavaan's start argument.

#' Level-specific starting values for multilevel CFA
#'
#' Computes level-specific starting values for two-level CFA models,
#' addressing the problem where the default starting values use the pooled
#' (total) covariance for all level-specific parameters. This causes
#' between-level uniquenesses and loadings to be on the wrong numerical scale
#' when the intraclass correlation coefficient (ICC) is low.
#'
#' @param model A character string specifying the lavaan model syntax.
#'   The model must be a two-level CFA with fully explicit syntax (all
#'   loadings, variances, and uniquenesses declared).
#' @param data A data frame containing the observed variables.
#' @param cluster A character string giving the name of the column in
#'   \code{data} that identifies the clustering variable (e.g., school).
#' @param group An optional character string giving the name of the column
#'   in \code{data} that identifies the grouping variable.
#' @param group.label An optional character vector of group labels, in the
#'   order the groups should be processed. Must match the unique values of
#'   \code{data[[group]]}.
#' @param meanstructure Logical. If \code{TRUE} (default), intercepts are
#'   included and their starting values are set from level-specific means.
#' @param auto_fix_first Logical. If \code{TRUE} (default), the first
#'   loading of each factor is automatically fixed to 1. Must match the
#'   value used in the \code{cfa()} call.
#' @param verbose Logical. If \code{TRUE}, progress notes are printed.
#'
#' @return A numeric vector of starting values for all free parameters, in the
#'   order of the free rows of the parameter table, suitable for passing to the
#'   \code{start} argument of \code{\link{cfa}}.
#'
#' @details
#' The function computes within-level and between-level covariance matrices
#' from the clustered data, then runs a single-factor approximation (FABIN)
#' separately at each level to obtain loading starting values on the correct
#' scale. Factor variance starting values are derived via conditional
#' least-squares. Uniqueness starting values are set to the residual variance
#' (total minus factor-implied), floored at 1\% of the total variance.
#'
#' For between-only variables (appearing only at the between level),
#' the between covariance is estimated from the unweighted sample
#' covariance of cluster means. For unequal cluster sizes this is a
#' simplification of lavaan's size-weighted between covariance, but
#' is adequate for starting-value purposes.
#'
#' The \code{cfa()} call that uses the returned starting values must use
#' the same \code{meanstructure}, \code{auto.fix.first}, \code{group}, and
#' \code{group.label} settings, otherwise the start vector length will not
#' match the parameter table.
#'
#' @section Requirements:
#' \itemize{
#'   \item Two-level CFA only; no regressions or exogenous covariates.
#'   \item All modeled observed variables must exist as numeric columns in
#'     \code{data}.
#'   \item Cross-loadings are supported.
#'   \item Free latent covariances and latent means start at zero.
#'   \item Free observed residual covariances are not supported.
#'   \item The model must be fully explicit (do not rely on lavaan's
#'     auto-completion).
#'   \item Continuous indicators only. No categorical outcomes.
#'   \item At least 2 clusters per group.
#'   \item At least 10 valid observations per level and group (after
#'     removing missing values).
#' }
#'
#' @examples
#' library(lavaan)
#'
#' # Simulate two-level data
#' set.seed(123)
#' n_schools <- 50
#' n_students <- 20
#' school <- rep(1:n_schools, each = n_students)
#' n <- n_schools * n_students
#' f_b <- rnorm(n_schools, sd = sqrt(0.25))
#' f_w <- rnorm(n)
#' y1 <- 0.8 * (f_b[school] + f_w) + rnorm(n, sd = sqrt(0.6)) + 2
#' y2 <- 0.9 * (f_b[school] + f_w) + rnorm(n, sd = sqrt(0.5)) + 2.5
#' y3 <- 0.7 * (f_b[school] + f_w) + rnorm(n, sd = sqrt(0.7)) + 3
#' y4 <- 1.0 * (f_b[school] + f_w) + rnorm(n, sd = sqrt(0.4)) + 2
#' dat <- data.frame(y1, y2, y3, y4, school)
#'
#' model <- "
#'   level: 1
#'     f =~ y1 + y2 + y3 + y4
#'     y1 ~~ y1
#'     y2 ~~ y2
#'     y3 ~~ y3
#'     y4 ~~ y4
#'   level: 2
#'     f =~ y1 + y2 + y3 + y4
#'     y1 ~~ y1
#'     y2 ~~ y2
#'     y3 ~~ y3
#'     y4 ~~ y4
#'     y1 ~ 1
#'     y2 ~ 1
#'     y3 ~ 1
#'     y4 ~ 1
#'     f ~ 0*1
#' "
#'
#' starts <- lavMLstart(model, dat, cluster = "school")
#'
#' fit <- cfa(model, data = dat, cluster = "school",
#'            meanstructure = TRUE,
#'            start = starts, test = "none", baseline = FALSE)
#' summary(fit)
#'
#' @export
lavMLstart <- function(model, data, cluster, group = NULL,
                       group.label = NULL, meanstructure = TRUE,
                       auto_fix_first = TRUE, verbose = FALSE) {

  # --- Input validation ---
  if (!cluster %in% names(data)) {
    lav_msg_stop(gettextf("Cluster variable '%s' not found in data.", cluster))
  }
  if (any(is.na(data[[cluster]]))) {
    lav_msg_stop(gettext("Cluster variable contains missing values."))
  }
  n_clusters <- length(unique(data[[cluster]]))
  if (n_clusters < 2L) {
    lav_msg_stop(gettextf("Need at least 2 clusters, got %d.", n_clusters))
  }

  if (!is.null(group)) {
    if (!group %in% names(data)) {
      lav_msg_stop(gettextf("Group variable '%s' not found in data.", group))
    }
    if (any(is.na(data[[group]]))) {
      lav_msg_stop(gettext("Group variable contains missing values."))
    }
  }
  ngroups <- if (!is.null(group)) length(unique(data[[group]])) else 1L

  # --- Step 1: Build the target parameter table ---
  # Use cfa(do.fit = FALSE) to get the parameter table exactly as cfa() will
  # use it internally. This avoids discrepancies between lavaanify() and the
  # actual fit (e.g., auto-completed factor variances, first-loading handling).
  # do.fit = FALSE only parses the model and constructs the parameter table;
  # no optimization or moment computation is performed, so the full data is
  # used without cost.
  dummy_fit <- tryCatch(
    cfa(model, data = data,
        cluster = if (is.null(cluster)) NULL else cluster,
        group = if (is.null(group)) NULL else group,
        group.label = if (is.null(group)) NULL else group.label,
        estimator = "ML", meanstructure = meanstructure,
        auto.fix.first = auto_fix_first,
        test = "none", baseline = FALSE,
        fit.by.level = FALSE, do.fit = FALSE),
    error = function(e) e
  )
  if (inherits(dummy_fit, "error")) {
    lav_msg_stop(gettextf("Failed to construct the model: %s.",
                         conditionMessage(dummy_fit)))
  }
  target <- parTable(dummy_fit)
  if (!"group" %in% names(target) || is.null(target$group)) {
    target$group <- rep(1L, nrow(target))
  }
  free <- target$free > 0L
  n_free <- sum(free)
  nlevels <- if (length(target$level)) max(target$level) else 1L
  if (nlevels > 2L) {
    lav_msg_stop(gettext("Only two-level models are supported."))
  }
  if (nlevels < 2L) {
    lav_msg_stop(gettext("The model must have at least 2 levels."))
  }

  latent <- unique(target$lhs[target$op == "=~"])
  observed <- setdiff(
    unique(c(
      target$rhs[target$op == "=~"],
      target$lhs[target$op == "~~" & target$lhs == target$rhs],
      target$lhs[target$op == "~1"]
    )),
    latent
  )
  missing_items <- setdiff(observed, names(data))
  if (length(missing_items)) {
    lav_msg_stop(gettextf("Modeled indicators missing from data: %s.",
                         lav_msg_view(missing_items)))
  }
  nonnumeric <- observed[!vapply(data[observed], is.numeric, logical(1))]
  if (length(nonnumeric)) {
    lav_msg_stop(gettextf("Modeled indicators must be numeric: %s.",
                         lav_msg_view(nonnumeric)))
  }

  if (is.null(target$group)) {
    target$group <- rep(1L, nrow(target))
  }

  # Check for unsupported model features
  if (any(target$op == "~" & !target$rhs %in% latent)) {
    lav_msg_stop(gettext("Models with exogenous covariates are not supported."))
  }
  if (any(target$op == "~")) {
    lav_msg_stop(gettext("Regression coefficients are not supported."))
  }
  if (any(target$op == "|")) {
    lav_msg_stop(gettext("Categorical/ordinal outcomes are not supported."))
  }
  if (any(free & target$op == "~~" & target$lhs != target$rhs &
          target$lhs %in% observed & target$rhs %in% observed)) {
    lav_msg_stop(gettextf(
      "Free observed residual covariances (e.g., %s) are not yet supported by %s. Please remove them or use the default starting values.",
      "y1 ~~ y2", "lavMLstart()"))
  }

  if (verbose) {
    lav_msg_note(gettextf("Free parameters: %d, group(s): %d, level(s): %d.",
                         n_free, ngroups, nlevels))
  }

  # --- Step 2: Split data by group and compute level-specific moments ---
  if (!is.null(group)) {
    if (is.null(group.label)) {
      group.label <- as.character(unique(data[[group]]))
    } else {
      group.label <- as.character(group.label)
      if (length(unique(group.label)) != length(group.label)) {
        lav_msg_stop(gettextf("The %s argument must have unique values.",
                             "group.label="))
      }
      if (length(group.label) != ngroups) {
        lav_msg_stop(gettextf("The %s argument length (%d) does not match number
                              of groups in data (%d).",
                              "group.label=", length(group.label), ngroups))
      }
      if (!setequal(group.label, as.character(unique(data[[group]])))) {
        lav_msg_stop(gettextf("The %s argument values do not match data group values.",
                             "group.label="))
      }
    }
    group_data <- lapply(group.label, function(g) data[data[[group]] == g, ])
    names(group_data) <- group.label
    for (i in seq_along(group_data)) {
      if (nrow(group_data[[i]]) == 0L) {
        lav_msg_stop(gettextf("Group '%s' has no data.",
                             names(group_data)[i]))
      }
    }
  } else {
    group_data <- list(data)
    names(group_data) <- "1"
  }

  level_items <- lapply(seq_len(ngroups), function(g) {
    lapply(1:2, function(l) {
      rows <- target$group == g & target$level == l
      intersect(observed, unique(c(
        target$rhs[rows & target$op == "=~"],
        target$lhs[rows & target$op %in% c("~~", "~1")],
        target$rhs[rows & target$op == "~~"]
      )))
    })
  })
  group_moments <- lapply(seq_len(ngroups), function(g) {
    sets <- level_items[[g]]
    if (any(lengths(sets) == 0L)) {
      lav_msg_stop(gettext("Each level must contain observed indicators."))
    }
    w <- lav_ml_compute_level_moments(group_data[[g]], cluster,
                                       items = sets[[1]])
    b <- lav_ml_compute_level_moments(group_data[[g]], cluster,
                                       items = sets[[2]],
                                       between_only = setdiff(sets[[2]],
                                                              sets[[1]]))
    items <- unique(unlist(sets))
    expand <- function(sigma, selected) {
      out <- matrix(NA_real_, length(items), length(items))
      pos <- match(selected, items)
      out[pos, pos] <- sigma
      out
    }
    mu_w <- mu_b <- rep(NA_real_, length(items))
    mu_w[match(w$items, items)] <- w$overall_mean
    mu_b[match(b$items, items)] <- b$mu_b[[1L]]
    list(items = items, sigma_w = expand(w$sigma_w, w$items),
         sigma_b = list(expand(b$sigma_b[[1L]], b$items)),
         overall_mean = mu_w, mu_b = list(mu_b),
         J = min(w$J, b$J), df_within = w$df_within)
  })
  names(group_moments) <- names(group_data)

  for (i in seq_along(group_moments)) {
    m <- group_moments[[i]]
    if (m$J < 2L) {
      lav_msg_stop(gettextf("Group '%s' has fewer than 2 clusters after
                            removing NAs.", names(group_moments)[i]))
    }
    if (m$df_within < 1L) {
      lav_msg_stop(gettextf("Group '%s' has no within-cluster degrees of
                            freedom (all clusters size 1).",
                            names(group_moments)[i]))
    }
  }

  # --- Step 3: Fill starts per group ---
  fill_starts <- function(g) {
    values <- rep(NA_real_, nrow(target))
    fixed <- target$group == g & target$free == 0L & is.finite(target$ustart)
    values[fixed] <- target$ustart[fixed]
    moments <- group_moments[[g]]
    items_all <- moments$items

    # Within-level fabin for all factors (fallback for between-level)
    within_fabin <- list()
    if (nlevels >= 1L) {
      sigma_w <- moments$sigma_w
      g_1_idx <- which(target$group == g & target$level == 1L)
      if (length(g_1_idx)) {
        lv_names_w <- unique(target$lhs[g_1_idx][target$op[g_1_idx] == "=~"])
        for (f in lv_names_w) {
          load_idx_w <- which(target$group == g & target$level == 1L &
                               target$lhs == f & target$op == "=~")
          if (!length(load_idx_w)) next
          ind_w <- target$rhs[load_idx_w]
          pos_w <- match(ind_w, items_all)
          if (any(is.na(pos_w))) next
          cov_w <- sigma_w[pos_w, pos_w, drop = FALSE]
          fab_w <- tryCatch(
            lav_cfa_1fac_fabin(cov_w, std_lv = FALSE,
                                lambda_only = TRUE, method = "fabin3"),
            error = function(e) NULL
          )
          if (!is.null(fab_w$lambda) && all(is.finite(fab_w$lambda))) {
            within_fabin[[f]] <- list(lambda = fab_w$lambda, ind = ind_w)
          }
        }
      }
    }

    for (l in seq_len(nlevels)) {
      sigma_l <- if (l == 1L) moments$sigma_w else moments$sigma_b[[l - 1L]]
      mu_l <- if (l == 1L) NULL else moments$mu_b[[l - 1L]]
      if (is.null(sigma_l)) next

      g_l_idx <- which(target$group == g & target$level == l)
      if (!length(g_l_idx)) next

      lv_names <- unique(target$lhs[g_l_idx][target$op[g_l_idx] == "=~"])

      for (f in lv_names) {
        load_idx <- which(target$group == g & target$level == l &
                           target$lhs == f & target$op == "=~")
        if (!length(load_idx)) next

        ind_names <- target$rhs[load_idx]
        ind_pos <- match(ind_names, items_all)
        if (any(is.na(ind_pos))) next

        cov_f <- sigma_l[ind_pos, ind_pos, drop = FALSE]

        # Determine if covariance is usable for fabin
        eig_raw <- eigen(cov_f, symmetric = TRUE, only.values = TRUE)$values
        min_eig <- min(eig_raw)
        max_diag <- max(diag(cov_f))
        cov_usable <- (min_eig > -0.05 * max(max_diag, 1e-6)) &&
                      all(diag(cov_f) > 0.01 * max(max_diag, 1e-6))

        # PD repair
        if (min_eig < 0) {
          eps <- 0.01 * max(max_diag, 1e-6)
          ridge <- abs(min_eig) + eps
          cov_f <- cov_f + diag(length(eig_raw)) * ridge
        }
        d <- diag(cov_f)
        d[d <= 0] <- 0.01 * max(d[d > 0], 1e-4)
        diag(cov_f) <- d

        # Read actual fixed variance value
        var_idx <- which(target$group == g & target$level == l &
                         target$lhs == f & target$op == "~~" & target$rhs == f)
        std_lv <- FALSE
        fixed_var <- NULL
        if (length(var_idx) == 1L) {
          if (target$free[var_idx] == 0L) {
            fixed_var <- lav_ml_extract_fixed_value(target, var_idx)
            if (!is.null(fixed_var) && abs(fixed_var - 1) < 1e-10) {
              std_lv <- TRUE
            }
          }
        }

        # Run fabin on the level-specific covariance
        lambda_start <- NULL
        if (cov_usable) {
          fabin_result <- tryCatch(
            lav_cfa_1fac_fabin(cov_f, std_lv = FALSE,
                                lambda_only = TRUE, method = "fabin3"),
            error = function(e) NULL
          )
          if (!is.null(fabin_result$lambda) &&
              all(is.finite(fabin_result$lambda)) &&
              all(abs(fabin_result$lambda) < 10)) {
            lambda_start <- fabin_result$lambda
          }
        }

        # Fallback: within-informed loadings
        if (is.null(lambda_start) && l > 1L) {
          for (wf_name in names(within_fabin)) {
            wf <- within_fabin[[wf_name]]
            if (setequal(wf$ind, ind_names)) {
              lambda_start <- wf$lambda[match(ind_names, wf$ind)]
              if (verbose) {
                lav_msg_note(gettextf("Using within-informed loadings (from
                                      '%s') for factor '%s' g%d l%d.",
                                      wf_name, f, g, l))
              }
              break
            }
          }
        }

        if (is.null(lambda_start)) {
          lambda_start <- c(1, rep(0.8, length(ind_names) - 1))
          if (verbose) {
            lav_msg_note(gettextf("No fabin available for factor '%s' group
                                  %d level %d; using generic fallback.",
                                  f, g, l))
          }
        }

        # Rescale for non-unit fixed variance
        if (!is.null(fixed_var) && fixed_var > 0) {
          psi_implied <- lav_ml_compute_var_start(lambda_start, cov_f,
                                                   std_lv = FALSE)
          if (psi_implied > 1e-10) {
            lambda_start <- lambda_start * sqrt(psi_implied / fixed_var)
          }
        }

        # Override fixed loadings
        fixed_load_idx <- load_idx[target$free[load_idx] == 0L]
        if (is.null(fixed_var) && length(fixed_load_idx)) {
          for (fi in fixed_load_idx) {
            fv <- target$ustart[fi]
            pos_in_lambda <- match(fi, load_idx)
            if (is.finite(fv) && abs(fv) > 1e-10 &&
                abs(lambda_start[pos_in_lambda]) > 1e-10) {
              lambda_start <- lambda_start * fv / lambda_start[pos_in_lambda]
              break
            }
          }
        }
        if (length(fixed_load_idx) > 0L) {
          for (fi in fixed_load_idx) {
            fv <- target$ustart[fi]
            if (!is.na(fv) && is.finite(fv)) {
              pos_in_lambda <- match(fi, load_idx)
              lambda_start[pos_in_lambda] <- fv
            }
          }
        }

        values[load_idx] <- lambda_start

        # Factor variance start
        free_var_idx <- which(target$group == g & target$level == l &
                              target$lhs == f & target$op == "~~" &
                              target$rhs == f & target$free > 0L)
        if (length(free_var_idx) == 1L) {
          values[free_var_idx] <- lav_ml_compute_var_start(
            lambda_start, cov_f, std_lv, verbose = verbose
          )
        }
      }

      # Latent covariances
      latent_cov <- g_l_idx[target$op[g_l_idx] == "~~" &
                            target$lhs[g_l_idx] %in% lv_names &
                            target$rhs[g_l_idx] %in% lv_names &
                            target$free[g_l_idx] > 0L]
      for (ri in latent_cov[!is.finite(values[latent_cov])]) {
        values[ri] <- if (target$lhs[ri] == target$rhs[ri]) {
          max(0.5 * min(diag(sigma_l)[match(level_items[[g]][[l]],
                                            items_all)]), 1e-6)
        } else 0
      }

      # Intercept starts
      int_idx <- which(target$group == g & target$level == l &
                       target$op == "~1" & target$free > 0L &
                       !target$lhs %in% lv_names)
      values[g_l_idx[target$op[g_l_idx] == "~1" &
                       target$lhs[g_l_idx] %in% lv_names &
                       target$free[g_l_idx] > 0L]] <- 0
      if (length(int_idx) > 0L) {
        int_names <- target$lhs[int_idx]
        if (!is.null(mu_l)) {
          int_pos <- match(int_names, items_all)
          values[int_idx] <- mu_l[int_pos]
        } else {
          values[int_idx] <- moments$overall_mean[match(int_names, items_all)]
        }
      }

      # Standalone variances
      all_uniq_idx <- which(target$group == g & target$level == l &
                            target$op == "~~" & target$lhs == target$rhs &
                            target$free > 0L)
      factor_inds <- unlist(lapply(lv_names, function(f) {
        target$rhs[which(target$group == g & target$level == l &
                         target$lhs == f & target$op == "=~")]
      }))
      standalone <- all_uniq_idx[!target$lhs[all_uniq_idx] %in% factor_inds]
      if (length(standalone) > 0L) {
        for (si in standalone) {
          pos <- match(target$lhs[si], items_all)
          if (!is.na(pos)) {
            v <- sigma_l[pos, pos]
            values[si] <- max(0.5 * abs(v), 0.01 * max(abs(v), 1e-4))
          }
        }
      }
    }
    values
  }

  # --- Step 4: Fill starts for each group ---
  all_values <- rep(NA_real_, nrow(target))
  for (g in seq_len(ngroups)) {
    g_values <- fill_starts(g)
    g_rows <- which(target$group == g)
    all_values[g_rows] <- g_values[g_rows]
  }

  # --- Step 5: Handle cross-group label sharing ---
  all_labels <- target$label[free]
  valid_labels <- all_labels[nzchar(all_labels) & !is.na(all_labels)]
  shared_labels <- valid_labels[duplicated(valid_labels) |
                                  valid_labels %in%
                                  rev(valid_labels[duplicated(valid_labels)])]
  label_changed <- FALSE
  for (lbl in unique(shared_labels)) {
    lbl_rows <- which(free & target$label == lbl)
    if (length(lbl_rows) < 2L) next
    finite_rows <- lbl_rows[is.finite(all_values[lbl_rows])]
    if (length(finite_rows) > 0L) {
      first_finite <- finite_rows[1L]
      if (any(abs(all_values[lbl_rows] - all_values[first_finite]) > 1e-12,
              na.rm = TRUE)) {
        label_changed <- TRUE
      }
      all_values[lbl_rows] <- all_values[first_finite]
    }
  }

  # --- Step 5b: Recompute dependent starts if labels changed loadings ---
  if (label_changed) {
    for (g in seq_len(ngroups)) {
      moments <- group_moments[[g]]
      items_all <- moments$items
      for (l in seq_len(nlevels)) {
        sigma_l <- if (l == 1L) moments$sigma_w else moments$sigma_b[[l - 1L]]
        if (is.null(sigma_l)) next
        g_l_idx <- which(target$group == g & target$level == l)
        if (!length(g_l_idx)) next
        lv_names <- unique(target$lhs[g_l_idx][target$op[g_l_idx] == "=~"])
        for (f in lv_names) {
          load_idx <- which(target$group == g & target$level == l &
                            target$lhs == f & target$op == "=~")
          if (!length(load_idx)) next
          ind_names <- target$rhs[load_idx]
          ind_pos <- match(ind_names, items_all)
          if (any(is.na(ind_pos))) next
          cov_f <- sigma_l[ind_pos, ind_pos, drop = FALSE]

          lambda_final <- all_values[load_idx]
          if (any(!is.finite(lambda_final))) next

          var_idx <- which(target$group == g & target$level == l &
                           target$lhs == f & target$op == "~~" &
                           target$rhs == f)
          fixed_var <- NULL
          std_lv <- FALSE
          if (length(var_idx) == 1L && target$free[var_idx] == 0L) {
            fixed_var <- lav_ml_extract_fixed_value(target, var_idx)
            if (!is.null(fixed_var) && abs(fixed_var - 1) < 1e-10) {
              std_lv <- TRUE
            }
          }
          free_var_idx <- var_idx[target$free[var_idx] > 0L]
          if (length(free_var_idx) == 1L) {
            all_values[free_var_idx] <- lav_ml_compute_var_start(
              lambda_final, cov_f, std_lv, verbose = FALSE
            )
          }
        }
      }
    }
    for (lbl in unique(shared_labels)) {
      lbl_rows <- which(free & target$label == lbl)
      if (length(lbl_rows) < 2L) next
      finite_rows <- lbl_rows[is.finite(all_values[lbl_rows])]
      if (length(finite_rows) > 0L) {
        all_values[lbl_rows] <- all_values[finite_rows[1L]]
      }
    }
  }

  # Compute residuals with final constrained loadings
  for (g in seq_len(ngroups)) {
    moments <- group_moments[[g]]
    for (l in 1:2) {
      rows <- which(target$group == g & target$level == l)
      loads <- rows[target$op[rows] == "=~"]
      factors <- unique(target$lhs[loads])
      if (!length(factors)) next
      items <- level_items[[g]][[l]]
      lambda_mat <- matrix(0, length(items), length(factors),
                       dimnames = list(items, factors))
      lambda_mat[cbind(match(target$rhs[loads], items),
                   match(target$lhs[loads], factors))] <- all_values[loads]
      psi_mat <- matrix(0, length(factors), length(factors),
                    dimnames = list(factors, factors))
      covs <- rows[target$op[rows] == "~~" &
                     target$lhs[rows] %in% factors &
                     target$rhs[rows] %in% factors]
      for (ri in covs) {
        v <- all_values[ri]
        if (!is.finite(v)) {
          lav_msg_stop(gettextf("Latent covariance has no finite start: %s ~~ %s.",
                                target$lhs[ri], target$rhs[ri]))
        }
        psi_mat[target$lhs[ri], target$rhs[ri]] <- v
        psi_mat[target$rhs[ri], target$lhs[ri]] <- v
      }
      absent <- setdiff(factors,
                        target$lhs[covs[target$lhs[covs] == target$rhs[covs]]])
      if (length(absent)) {
        lav_msg_stop(gettextf("The model must explicitly declare latent variances: %s.",
                             lav_msg_view(absent)))
      }
      implied <- diag(lambda_mat %*% psi_mat %*% t(lambda_mat))
      sigma <- if (l == 1L) moments$sigma_w else moments$sigma_b[[1L]]
      uniq <- rows[target$op[rows] == "~~" &
                    target$lhs[rows] == target$rhs[rows] &
                    target$lhs[rows] %in% target$rhs[loads] &
                    target$free[rows] > 0L]
      for (ri in uniq) {
        item <- target$lhs[ri]
        total <- sigma[match(item, moments$items), match(item, moments$items)]
        all_values[ri] <- max(total - implied[match(item, items)],
                              0.01 * max(abs(total), 1e-6))
      }
    }
  }
  # Residual equality constraints
  for (lbl in unique(shared_labels)) {
    rows <- which(free & target$label == lbl)
    all_values[rows] <- all_values[rows[1L]]
  }

  # --- Step 6: Final checks ---
  free_values <- all_values[free]
  if (any(!is.finite(free_values))) {
    rows <- which(free)[!is.finite(free_values)]
    lav_msg_stop(gettextf("Unsupported or uninitialized parameters: %s.",
                         lav_msg_view(sprintf("%s %s %s [group %s, level %s]",
                                       target$lhs[rows], target$op[rows],
                                       target$rhs[rows], target$group[rows],
                                        target$level[rows]))))
  }

  var_free_idx <- which(free & target$op == "~~" & target$lhs == target$rhs)
  if (length(var_free_idx) && any(free_values[
    match(var_free_idx, which(free))] <= 0)) {
    if (verbose) {
      lav_msg_note(gettext("Some variance starts are non-positive."))
    }
  }

  free_values
}


# --- Internal: compute within and between level moments from clustered data ---
lav_ml_compute_level_moments <- function(data, cluster, items,
                                          between_only = character()) {
  cluster_var <- data[[cluster]]
  if (any(is.na(cluster_var))) {
    lav_msg_stop(gettext("Cluster variable has NA values."))
  }

  if (!length(items)) {
    lav_msg_stop(gettext("The model has no observed indicators."))
  }

  y <- as.matrix(data[, items, drop = FALSE])
  cl <- cluster_var

  ok <- complete.cases(y) & apply(is.finite(y), 1, all)
  if (sum(ok) < 10L) {
    lav_msg_stop(gettext("Fewer than 10 valid observations after removing NAs."))
  }
  y <- y[ok, , drop = FALSE]
  cl <- cl[ok]

  clusters <- sort(unique(cl))
  J <- length(clusters)
  if (J < 2L) {
    lav_msg_stop(gettextf("Need at least 2 clusters, got %d.", J))
  }
  n_vec <- vapply(clusters, function(c) sum(cl == c), integer(1))
  N <- sum(n_vec)

  overall_mean <- colMeans(y)

  # Within-level covariance (pooled)
  within_ss <- matrix(0, ncol(y), ncol(y))
  df_within <- 0L
  for (c in clusters) {
    rows <- which(cl == c)
    nc <- length(rows)
    if (nc <= 1L) next
    yc <- y[rows, , drop = FALSE]
    yc_c <- sweep(yc, 2L, colMeans(yc))
    within_ss <- within_ss + crossprod(yc_c)
    df_within <- df_within + (nc - 1L)
  }
  if (df_within < 1L && !all(items %in% between_only)) {
    lav_msg_stop(gettext("No within-cluster df (all clusters have size 1)."))
  }
  sigma_w <- if (df_within > 0L) within_ss / df_within else within_ss

  # Between-level covariance (matching lavaan's lav_samp_cl_patterns)
  n_items <- ncol(y)
  cluster_means <- matrix(NA, nrow = J, ncol = n_items)
  for (j in seq_len(J)) {
    cluster_means[j, ] <- colMeans(y[cl == clusters[j], , drop = FALSE])
  }

  cluster_means_c <- sweep(cluster_means, 2L, overall_mean)

  s_corr <- (N - sum(n_vec^2) / N) / (J - 1L)

  between_ss <- matrix(0, n_items, n_items)
  for (j in seq_len(J)) {
    between_ss <- between_ss + n_vec[j] *
      outer(cluster_means_c[j, ], cluster_means_c[j, ])
  }
  s_b <- between_ss / (J - 1L)

  sigma_b <- (s_b - sigma_w) / s_corr

  mu_b <- as.numeric(n_vec %*% cluster_means / sum(n_vec))

  only <- match(between_only, items)
  if (length(only)) {
    raw_between <- stats::cov(cluster_means)
    sigma_b[only, ] <- raw_between[only, , drop = FALSE]
    sigma_b[, only] <- raw_between[, only, drop = FALSE]
    mu_b[only] <- colMeans(cluster_means)[only]
  }

  list(
    items = items,
    sigma_w = sigma_w,
    sigma_b = list(sigma_b),
    mu_b = list(mu_b),
    overall_mean = overall_mean,
    J = J,
    df_within = df_within
  )
}


# --- Internal: extract a fixed parameter value from the lavaanify table ---
lav_ml_extract_fixed_value <- function(target, idx) {
  if (length(idx) != 1L) return(NULL)
  u <- target$ustart[idx]
  if (!is.na(u) && is.finite(u)) {
    return(as.numeric(u))
  }
  1.0
}


# --- Internal: compute a data-driven factor variance start ---
# Conditional least-squares: psi = sum_{i<j} lambda_i*lambda_j*sigma_ij
#                               / sum_{i<j} (lambda_i*lambda_j)^2
# Scale-consistent: rescaling loadings by 'a' rescales psi by 1/a^2.
lav_ml_compute_var_start <- function(lambda, cov_f, std_lv, verbose = FALSE) {
  n <- length(lambda)
  if (n < 2L) {
    if (abs(lambda[1L]) > 1e-10) {
      return(max(cov_f[1L, 1L] / lambda[1L]^2, 1e-6))
    }
    return(0.05)
  }

  if (std_lv) return(1.0)

  if (sum(lambda^2) < 1e-20) {
    return(0.05)
  }

  i_idx <- combn(n, 2, simplify = FALSE)
  numer <- 0
  denom <- 0
  for (pair in i_idx) {
    i <- pair[1L]
    j <- pair[2L]
    lambda_prod <- lambda[i] * lambda[j]
    if (abs(lambda_prod) < 1e-14) next
    numer <- numer + lambda_prod * cov_f[i, j]
    denom <- denom + lambda_prod^2
  }

  if (denom < 1e-20) {
    psi <- mean(diag(cov_f)) / mean(lambda^2)
  } else {
    psi <- numer / denom
  }

  if (!is.finite(psi) || psi <= 0) {
    psi <- max(mean(diag(cov_f)) / max(sum(lambda^2), 1e-10), 1e-6)
  }

  scale_ref <- max(mean(diag(cov_f)) / max(mean(lambda^2), 1e-10), 1e-6)
  psi_lower <- 1e-4 * scale_ref
  psi_upper <- 100 * scale_ref
  psi_clipped <- min(max(psi, psi_lower), psi_upper)

  if (verbose && abs(psi_clipped - psi) > 1e-10) {
    lav_msg_note(gettextf("Factor variance start clipped from %g to %g.",
                         psi, psi_clipped))
  }

  psi_clipped
}
