box::use(
    tibble[tibble],
    rlang[
        quo, expr, eval_tidy, syms
    ],
    dplyr[group_modify, ungroup, pull, select],
    purrr[map, map_dbl, map_dfr, pmap, iwalk],
    stats[na.omit, wilcox.test],
    tidyselect[where, eval_select, everything],
    tidyr[expand_grid, pivot_wider, unnest],
    tabhelpers[table_default, table_summary],
    cli,
    ../utils[all_pairs, format_pval]
)

wilcoxon_test <- function(data, ..., alternative = "two.sided", mu = 0,
                          paired = FALSE, correct = FALSE, exact = NULL,
                          combine = FALSE, format = FALSE,
                          na_action = "na.omit", permuted_pairs = "lt") {
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

    names(mu_vec) <- names(cols)
    names(alt_vec) <- names(cols)

    if (combine) {
        vars <- map(all_pairs(names(cols), simplify = F, direction = permuted_pairs), syms)

        out <-
            data |>
            group_modify( ~ {
                map_dfr(
                    vars, function(.var) {
                        mu_diff <- mu_vec[as.character(.var[[1]])] - mu_vec[as.character(.var[[2]])]

                        alt <- alt_vec[as.character(.var[[1]])]

                        test <- expr(
                            wilcox.test(
                                !!!.var,
                                alternative = alt,
                                mu = mu_diff,
                                paired = {{paired}},
                                exact = {{exact}},
                                correct = {{correct}},
                                na.action = {{na_action}}
                            )
                        )
                        test <- eval_tidy(test, data = .x)

                        tibble(
                            x = as.character(.var[[1]]),
                            y = as.character(.var[[2]]),
                            alternative = test$alternative,
                            statistics = unname(test$statistic),
                            pvals = if (format) as.factor(format_pval(test$p.value)) else round(test$p.value, digits = 3)
                        )
                    }
                )
            }) |>
            ungroup()

        k_test <- "Two-Sample"

    } else {
        vars <- syms(names(cols))

        out <-
            data |>
            group_modify( ~ {
                map_dfr(
                    vars, function(.var) {
                        var_name <- as.character(.var)
                        x <- eval_tidy(expr(pull(.x, !!.var)))

                        test <- wilcox.test(
                            x,
                            alternative = alt_vec[var_name],
                            mu = mu_vec[var_name],
                            paired = paired,
                            exact = exact,
                            correct = correct,
                            na.action = na_action
                        )

                        m <- mean(x)

                        tibble(
                            x = var_name,
                            true_mu = mu_vec[var_name],
                            alternative = test$alternative,
                            statistics = unname(test$statistic),
                            pvals = if (format) as.factor(format_pval(test$p.value)) else round(test$p.value, digits = 3)
                        )
                    }
                )
            }) |>
            ungroup()

        k_test <- "One-Sample"
    }

    res <- list(out = out, test_type = k_test, combine = combine, paired = paired)
    class(res) <- c("basestats", "wilcoxtest")
    res
}

print.wilcoxtest <- function(x, ...) {
    data <- x$out
    test_type <- x$test_type
    paired <- x$paired
    combine <- x$combine

    if (test_type == "One-Sample") {
        title <- "One-Sample Wilcoxon Signed Rank Test"
        test_stat_name <- "V-stat"
    } else if (test_type == "Two-Sample" && paired) {
        title <- "Paired Wilcoxon Signed Rank Test"
        test_stat_name <- "V-stat"
    } else {
        title <- "Wilcoxon Rank Sum Test (Mann-Whitney U Test)"
        test_stat_name <- "W-stat"
    }

    cat("\n")
    cli::rule(center = title) |>
        cli::style_bold() |>
        cli::cat_line("\n\n")

    data |>
        pmap(list) |>
        iwalk(function (res, i) {
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

            alt_text <- switch(res$alternative,
                               "two.sided" = "different from",
                               "less" = "less than",
                               "greater" = "greater than")

            cat(comparison_text, "\n")

            result_tbl <- tibble::tibble(
                !!test_stat_name := res$statistics,
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
                )
            )
            cat("\n\n")
        })
}

box::register_S3_method("print", "wilcoxtest")
