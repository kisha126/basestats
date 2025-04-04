box::use(
    tibble[tibble],
    rlang[
        quo, expr, eval_tidy, syms, quo_is_null, enquo
    ],
    dplyr[group_modify, ungroup, pull, select, group_vars],
    purrr[map, map_dbl, map_dfr, pmap, iwalk, map_chr],
    stats[na.omit, ks.test, pnorm],
    tidyselect[where, eval_select, everything],
    tidyr[expand_grid, pivot_wider, unnest],
    tabhelpers[table_default, table_summary],
    cli,
    ../utils[all_pairs]
)

#' Modified Kolmogorov-Smirnov Test
#'
#' @examples
#' set.seed(123)
#' df_ex <- data.frame(
#'     x1 = rnorm(100, mean = 100, sd = 14),
#'     x2 = rt(100, df = 99, ncp = 100),
#'     y1 = rchisq(100, df = 100),
#'     y2 = rgamma(100, 51, 0.5),
#'     g = sample(rep(letters[1:5], each = 20))
#' )
#'
#' box::use(bs = ./basestats)
#'
#' # One-Sample
#' df_ex |>
#'     bs$npar$ks_test(x1, x2)
#'
#' # Two-Sample
#' df_ex |>
#'     bs$npar$ks_test(x1, x2, y = c(y1, y2))
#'
#'
#' @export
ks_test <- function(data, ..., y = NULL, alternative = "two.sided",
                    exact = NULL, simulate.p_value = FALSE, B = 2000) {

    cols <- eval_select(quo(c(...)), data = data)
    vars <- syms(names(cols))

    y_quo <- enquo(y)

    if (!quo_is_null(y_quo)) {
        y_expr <- rlang::quo_get_expr(y_quo)
        if (is.call(y_expr) || is.function(y_expr) || is.list(y_expr) || rlang::is_formula(y_expr)) {
            if (is.call(y_expr) &&
                (y_expr[[1]] == quote(c) ||
                 y_expr[[1]] == quote(list) ||
                 y_expr[[1]] == quote(cbind))) {
                col_names <- all.vars(y_expr)
                y_cols <- tidyselect::eval_select(rlang::parse_expr(paste0("c(", paste(col_names, collapse = ","), ")")), data)
                y_list <- syms(names(y_cols))
            } else {
                # Handle anonymous functions, formula functions, or list of functions
                if (rlang::is_formula(y_expr)) {

                    y_expr <- rlang::new_formula(
                        lhs = NULL,
                        rhs = y_expr[[2]],
                        env = rlang::quo_get_env(y_quo)
                    )

                    y_list <- list(rlang::as_function(y_expr))  # Convert formula to function
                } else {
                    y_list <- list(rlang::eval_tidy(y_quo))
                }
            }
        } else {
            y_cols <- tidyselect::eval_select(y_quo, data)
            y_list <- syms(names(y_cols))
        }
    } else {
        y_list <- NULL
    }

    out <-
        data |>
        group_modify( ~ {
            purrr::map_dfr(
                seq_along(vars), function(idx) {
                    .var <- vars[[idx]]

                    if (rlang::quo_is_null(y_quo)) {
                        test <- expr(
                            ks.test(
                                !!.var,
                                "pnorm",
                                alternative = {{alternative}},
                                exact = {{exact}},
                                simulate.p.value = {{simulate.p_value}},
                                B = {{B}}
                            )
                        )

                        null_hyp <- "Approximately Normally Distributed"
                    } else if (length(y_list) == 1) {
                        # Single function or column for all variables
                        y_val <- y_list[[1]]
                        if (is.function(y_val) || is.list(y_val)) {
                            x_vals <- rlang::eval_tidy(expr(pull(.x, !!.var)))
                            if (is.list(y_val)) {
                                # Handle list of functions
                                y_val <- y_val[[((idx - 1) %% length(y_val)) + 1]]
                            }
                            test <- ks.test(
                                x_vals, y_val, alternative = alternative,
                                exact = exact,
                                simulate.p.value = simulate.p_value,
                                B = B
                            )
                            null_hyp <- glue::glue("\nThe sample `{deparse(.var)}` follows the theoretical distribution specified by CDF.\n")
                        } else {
                            test <- expr(
                                ks.test(
                                    !!.var,
                                    y = !!y_val,
                                    alternative = alternative,
                                    exact = exact,
                                    simulate.p.value = simulate.p_value,
                                    B = B
                                )
                            )
                            null_hyp <- glue::glue("The sample `{deparse(.var)}` follows the theoretical distribution specified by column `{deparse(y_val)}`.")
                        }
                    } else {
                        # Multiple functions or columns
                        if (length(vars) %% length(y_list) != 0) {
                            stop("The number of columns in `...` must be a multiple of the number of elements in `y`.")
                        }
                        k <- length(vars) / length(y_list)
                        y_val <- y_list[[((idx - 1) %/% k) + 1]]
                        if (is.function(y_val) || is.list(y_val)) {
                            x_vals <- rlang::eval_tidy(expr(pull(.x, !!.var)))
                            if (is.list(y_val)) {
                                # Handle list of functions
                                y_val <- y_val[[((idx - 1) %% length(y_val)) + 1]]
                            }
                            test <- ks.test(
                                x_vals, y_val, alternative = alternative,
                                exact = exact,
                                simulate.p.value = simulate.p_value,
                                B = B
                            )
                            null_hyp <- glue::glue("\nThe sample `{deparse(.var)}` follows the theoretical distribution specified by CDF.\n")
                        } else {
                            test <- expr(
                                ks.test(
                                    !!.var,
                                    y = !!y_val,
                                    alternative = alternative,
                                    exact = exact,
                                    simulate.p.value = simulate.p_value,
                                    B = B
                                )
                            )
                            null_hyp <- glue::glue("The sample `{deparse(.var)}` follows the theoretical distribution specified by column `{deparse(y_val)}`.")
                        }
                    }
                    test <- rlang::eval_tidy(test, data = .x)

                    tibble(
                        null_hyp = null_hyp,
                        columns = deparse(.var),
                        statistics = unname(test$statistic),
                        pval = unname(round(test$p.value, digits = 3))
                    )
                }
            )
        }) |>
        ungroup()

    res <- list(
        out = out,
        test_type = if (is.null(y)) "One-Sample" else "Two-Sample",
        group_varnames = group_vars(data)
    )
    class(res) <- c("basestats", "kstest")
    res
}

print.kstest <- function(x, ...) {
    data <- x$out
    test_type <- x$test_type

    if (test_type == "One-Sample") {
        title <- "One-Sample Kolmogorov-Smirnov Test"
    } else if (test_type == "Two-Sample") {
        title <- "Two-Sample Kolmogorov-Smirnov Test"
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
            var1 <- res$columns

            # null_text <- res$null_hyp
            #
            # cli$rule(center = null_text, line = "") |>
            #     cli$cat_line()

            cat(glue::glue("{i}. {var1}"), "\n")

            result_tbl <- tibble(
                d_statistic = res$statistics,
                `p-value` = res$pval
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
        })
}

box::register_S3_method("print", "kstest")
