box::use(
    tibble[tibble],
    rlang[
        quo, expr, eval_tidy, syms
    ],
    dplyr[group_modify, ungroup, pull, select, group_vars],
    purrr[map, map_dbl, map_dfr, pmap, iwalk, map_chr],
    stats[na.omit, t.test, sd],
    tidyselect[where, eval_select, everything],
    tidyr[expand_grid, pivot_wider, unnest],
    tabhelpers[table_default, table_summary],
    cli,
    ../utils[all_pairs]
)

#' @export
t_test <- function(data, ..., formula = NULL, alternative = "two.sided",
                   mu = 0, paired = FALSE, var.equal = FALSE, combine = FALSE,
                   format = FALSE, na_action = "na.omit", permuted_pairs = "lt") {

    cols <- eval_select(quo(c(...)), data = data)
    n_cols <- length(cols)

    is_divisor <- function(x) n_cols %% x == 0
    expand_vec <- function(x) {
        if (length(x) == 1) return(rep(x, n_cols))
        if (length(x) == n_cols) return(x |> unlist() |> as.vector())
        if (is_divisor(length(x))) {
            repeat_times <- n_cols / length(x)
            x <- x |> unlist() |> as.vector()
            return(rep(x, each = repeat_times))
        }
        stop(sprintf("Length of input (%d) must be either 1, equal to number of selected variables (%d), or a divisor of it",
                     length(x), n_cols))
    }

    mu_vec <- tryCatch(
        expand_vec(mu),
        error = function(e) stop("Invalid mu length: ", e$message)
    )

    alt_vec <- tryCatch(
        expand_vec(alternative),
        error = function(e) stop("Invalid alternative length: ", e$message)
    )

    na_action <- if (is.character(na_action)) get(na_action) else as_function(na_action)

    names(mu_vec) <- names(cols)
    names(alt_vec) <- names(cols)

    if (!is.null(formula)) {
        lhs <- f_lhs(formula)
        rhs <- f_rhs(formula)

        new_data <- data |>
            select(!!lhs, !!rhs) |>
            pivot_wider(
                names_from = !!rhs,
                values_from = !!lhs,
                values_fn = list
            ) |>
            unnest(cols = everything())

        return(
            t_test(
                new_data, where(is.numeric),
                alternative = alternative, mu = mu, paired = paired,
                var.equal = var.equal, combine = combine,
                format = format, na_action = na_action
            )
        )
    }

    if (combine) {
        vars <- map(all_pairs(names(cols), simplify = F, direction = permuted_pairs), syms)
        diff_mean <- function(a, b) mean(a) - mean(b)

        out <-
            data |>
            group_modify( ~ {
                map_dfr(
                    vars, function(.var) {
                        mu_diff <- mu_vec[as.character(.var[[1]])] - mu_vec[as.character(.var[[2]])]

                        alt <- alt_vec[as.character(.var[[1]])]

                        test <- expr(
                            t.test(
                                !!!.var,
                                alternative = alt,
                                mu = mu_diff,
                                paired = {{paired}},
                                var.equal = {{var.equal}},
                                na.action = {{na_action}}
                            )
                        )
                        test <- eval_tidy(test, data = .x)

                        m <- expr(diff_mean(!!!.var)) |> eval_tidy(data = data)

                        cohen_d <- if (paired) {
                            .var <- syms(.var)
                            diff_vals <- expr(!!.var[[1]] - !!.var[[2]]) |>
                                eval_tidy(data = .x)
                            mean(diff_vals) / sd(diff_vals)
                        } else {
                            x1 <- eval_tidy(expr(pull(.x, !!.var[[1]])))
                            x2 <- eval_tidy(expr(pull(.x, !!.var[[2]])))
                            m1 <- mean(x1)
                            m2 <- mean(x2)
                            s1 <- sd(x1)
                            s2 <- sd(x2)
                            n1 <- length(x1)
                            n2 <- length(x2)

                            if (var.equal) {
                                pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
                                (m1 - m2) / pooled_sd
                            } else {
                                (m1 - m2) / sqrt((s1^2 + s2^2) / 2)
                            }
                        }


                        tibble(
                            x = as.character(.var[[1]]),
                            y = as.character(.var[[2]]),
                            diff_mean = m,
                            true_mu_diff = mu_diff,
                            cohen = cohen_d,
                            alternative = test$alternative,
                            df = unname(test$parameter),
                            statistics = unname(test$statistic),
                            pvals = if (format) as.factor(format_pval(test$p.value)) else round(test$p.value, digits = 3)
                        )
                    }
                )
            }) |>
            ungroup()

        test_type <- "Two-Sample"

    } else {
        vars <- syms(names(cols))

        out <-
            data |>
            group_modify( ~ {
                map_dfr(
                    vars, function(.var) {
                        var_name <- deparse(.var)
                        x <- eval_tidy(expr(pull(.x, !!.var)))

                        test <- t.test(
                            x,
                            alternative = alt_vec[var_name],
                            mu = mu_vec[var_name],
                            paired = paired,
                            var.equal = var.equal
                        )

                        m <- mean(x)
                        sd_x <- sd(x)

                        tibble(
                            x = var_name,
                            mean = m,
                            sd = sd_x,
                            true_mu = mu_vec[var_name],
                            cohen = (mean - true_mu) / sd,
                            alternative = test$alternative,
                            df = unname(test$parameter),
                            statistics = unname(test$statistic),
                            pvals = if (format) as.factor(format_pval(test$p.value)) else round(test$p.value, digits = 3)
                        )
                    }
                )
            }) |>
            ungroup()

        test_type <- "One-Sample"
    }

    res <- list(
        out = out,
        test_type = test_type,
        combine = combine,
        paired = paired,
        group_varnames = group_vars(data)
    )
    class(res) <- c("basestats", "ttest")
    res
}

print.ttest <- function(x, ...) {
    data <- x$out
    test_type <- x$test_type
    paired <- x$paired
    combine <- x$combine

    if (test_type == "One-Sample") {
        title <- "One-Sample t-test"
    } else if (test_type == "Two-Sample" && paired) {
        title <- "2 Paired Sample t-test"
    } else {
        title <- "2 Independent Sample t-test"
    }

    cat("\n")
    cli$rule(center = title) |>
        cli$style_bold() |>
        cli$cat_line("\n\n")

    # Check if we have group data - look for non-numeric columns that aren't part of the test variables
    has_groups <- !is.null(x$group_varnames)

    data |>
        pmap(list) |>
        iwalk(function (res, i) {
            # Create base comparison text
            if (test_type == "One-Sample") {
                var1 <- res$x
                true_mu <- res$true_mu
                comparison_text <- glue::glue("{i}. {var1} (μ = {true_mu})")
            } else if (combine) {
                var1 <- res$x
                var2 <- res$y
                comparison_text <- glue::glue("{i}. {var1} vs. {var2}")
            } else {
                var1 <- res$x
                comparison_text <- glue::glue("{i}. {var1}")
            }

            # Add group information if available
            if (has_groups) {
                # Extract group information
                group_info <- setdiff(names(res), c("x", "y", "mean", "sd", "true_mu", "diff_mean",
                                                   "true_mu_diff", "cohen", "alternative", "df",
                                                   "statistics", "pvals"))

                # Format group information
                group_str <- ""
                if (length(group_info) > 0) {
                    group_values <- map_chr(group_info, ~paste(.x, "=", res[[.x]]))
                    group_str <- paste(" when", cli$style_italic(paste(group_values, collapse = ", ")))
                }

                # Add group information to comparison text
                comparison_text <- paste0(comparison_text, group_str)
            }

            alt_text <- switch(res$alternative,
                               "two.sided" = "different from",
                               "less" = "less than",
                               "greater" = "greater than")
            cohen_d <- res$cohen

            cat(comparison_text, "\n")

            result_tbl <- tibble(
                t_statistic = res$statistics,
                `p-value` = res$pvals
            )

            cat("\n\n")
            table_default(
                result_tbl,
                center_table = TRUE,
                n_space = 3,
                style_columns = list(
                    "p-value" = function (ctx) {
                        val <- as.numeric(ctx$formatted_value)
                        if (val < 0.05 & val >= 0.001) {
                            cli$col_red(val)
                        } else if (val < 0.001) {
                            rep_txt <- replace(val, val < 0.001, "<0.001")
                            cli$style_bold(rep_txt)
                        } else {
                            cli$style_italic(val)
                        }
                    }
                ),
                title = "Test Result"
            )
            cat("\n\n")
            table_summary(
                tibble(
                    x = c("Alternative Hypothesis", "True Population", "Estimate", "Cohen's d"),
                    y = c(
                        paste(alt_text, res$true_mu),
                        res$true_mu,
                        round(ifelse(combine, res$diff_mean, res$mean), digits = 3),
                        round(cohen_d, digits = 3)
                    )
                ),
                title = "Other Information",
                center_table = TRUE,
                style = list(sep = "  :  ")
            )
        })
}
# print.ttest <- function(x, ...) {
#     data <- x$out
#     test_type <- x$test_type
#     paired <- x$paired
#     combine <- x$combine
#
#     if (test_type == "One-Sample") {
#         title <- "One-Sample t-test"
#     } else if (test_type == "Two-Sample" && paired) {
#         title <- "2 Paired Sample t-test"
#     } else {
#         title <- "2 Independent Sample t-test"
#     }
#
#     cat("\n")
#     cli$rule(center = title) |>
#         cli$style_bold() |>
#         cli$cat_line("\n\n")
#
#     data |>
#         pmap(list) |>
#         iwalk(function (res, i) {
#             if (test_type == "One-Sample") {
#                 var1 <- res$x
#                 true_mu <- res$true_mu
#                 comparison_text <- glue::glue("{i}. {var1} (μ = {true_mu})")
#             } else if (combine) {
#                 var1 <- res$x
#                 var2 <- res$y
#                 comparison_text <- glue::glue("{i}. {var1} vs. {var2}")
#             } else {
#                 var1 <- res$x
#                 comparison_text <- glue::glue("{i}. {var1}")
#             }
#
#             alt_text <- switch(res$alternative,
#                                "two.sided" = "different from",
#                                "less" = "less than",
#                                "greater" = "greater than")
#             cohen_d <- res$cohen
#
#             cat(comparison_text, "\n")
#
#             result_tbl <- tibble(
#                 t_statistic = res$statistics,
#                 `p-value` = res$pvals
#             )
#
#             cat("\n\n")
#             table_default(
#                 result_tbl,
#                 center_table = TRUE,
#                 n_space = 3,
#                 style_columns = list(
#                     "p-value" = function (ctx) {
#                         val <- as.numeric(ctx$formatted_value)
#                         if (val < 0.05 & val >= 0.001) {
#                             cli$col_red(val)
#                         } else if (val < 0.001) {
#                             rep_txt <- replace(val, val < 0.001, "<0.001")
#                             cli$style_bold(rep_txt)
#                         } else {
#                             cli$style_italic(val)
#                         }
#                     }
#                 ),
#                 title = "Test Result"
#             )
#             cat("\n\n")
#             table_summary(
#                 tibble(
#                     x = c("Alternative Hypothesis", "True Population", "Estimate", "Cohen's d"),
#                     y = c(
#                         paste(alt_text, res$true_mu),
#                         res$true_mu,
#                         round(ifelse(combine, res$diff_mean, res$mean), digits = 3),
#                         round(cohen_d, digits = 3)
#                     )
#                 ),
#                 title = "Other Information",
#                 center_table = TRUE,
#                 style = list(sep = "  :  ")
#             )
#         })
# }

box::register_S3_method("print", "ttest")
