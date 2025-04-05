box::use(
    tidyselect[eval_select],
    rlang[
        quo, expr, quos, exprs, syms, enquo, enexpr, enquos,
        enexprs, eval_tidy, is_function, quo_is_null, quo_get_expr,
        new_environment, base_env, is_formula,
        f_rhs, f_lhs
    ],
    purrr[pmap, map_dfr, iwalk],
    dplyr[group_modify, pull, ungroup, group_vars],
    glue[glue],
    tibble[tibble],
    cli,
    stats[ks.test, pnorm],
    tabhelpers[table_default]
)

#' Modified Kolmogorov-Smirnov Test
#'
#' @description Provides a pipe-friendly framework to performs Kolmogorov-Smirnov test
#'   Support grouped data and multiple variables.
#'
#' @param data a data frame. Columns are variables.
#' @param ... One or more unquoted expressions (or variable names) separated by
#'   commas. Used to select a variable of interest.
#'
#' @param y If not provided, the test will perform normality with mean = 0 and sd = 1, by default. It can be anonymous
#'   function (\code{function}; \code{\(.)}; or purrr-style syntax), or a column to be compared.
#'
#' @param alternative Indicates the alternative hypothesis.
#' @param exact NULL or a logical indicating whether an exact p-value should be computed. See ‘Details’ for the
#'   meaning of NULL.
#'
#' @param simulate.p_value a logical indicating whether to compute p-values by Monte Carlo simulation. (Ignored for the
#'   one-sample test.)
#'
#' @param B an integer specifying the number of replicates used in the Monte Carlo test.
#'
#' @return  a list for the `kstest` method.
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

    cols <- tidyselect::eval_select(quo(c(...)), data = data)
    vars <- syms(names(cols))
    var_names <- names(cols)

    y_quo <- enquo(y)

    # Special handling for functions in y
    if (!rlang::quo_is_null(y_quo)) {
        y_expr <- rlang::quo_get_expr(y_quo)
        y_env <- rlang::quo_get_env(y_quo)

        # Check if y is a function or list of functions
        y_eval <- try(rlang::eval_tidy(y_quo, data = data), silent = TRUE)
        has_functions <- !inherits(y_eval, "try-error") &&
            (is.function(y_eval) || (is.list(y_eval) && any(sapply(y_eval, is.function))))

        if (has_functions) {
            # We'll handle functions separately

            # For deparsing functions
            deparse_fn <- function(fn) {
                if (is.function(fn)) {
                    fn_text <- deparse(fn)
                    if (length(fn_text) > 3) {
                        # Truncate long function definitions
                        paste0(paste(fn_text[1:3], collapse = " "), "...")
                    } else {
                        paste(fn_text, collapse = " ")
                    }
                } else {
                    "non-function"
                }
            }

            # Create a pairwise grid for all combinations
            if (is.function(y_eval)) {
                # If y is a single function, apply it to all columns
                fn_list <- list(y_eval)
                fn_deparsed <- list(deparse_fn(y_eval))
            } else if (is.list(y_eval)) {
                # If y is a list of functions, use all functions for all columns
                fn_list <- y_eval
                fn_deparsed <- lapply(y_eval, deparse_fn)
            }

            # Create all combinations of columns and functions
            combinations <- tidyr::expand_grid(
                var_idx = seq_along(vars),
                fn_idx = seq_along(fn_list)
            )

            out <- data |>
                group_modify(~ {
                    purrr::map_dfr(
                        seq_len(nrow(combinations)), function(combo_idx) {
                            var_idx <- combinations$var_idx[combo_idx]
                            fn_idx <- combinations$fn_idx[combo_idx]

                            var <- vars[[var_idx]]
                            var_name <- var_names[var_idx]
                            fn <- fn_list[[fn_idx]]
                            fn_text <- fn_deparsed[[fn_idx]]

                            # Extract the column data
                            x_data <- .x[[var_name]]

                            # Create a proper distribution function that ks.test can use
                            dist_fn <- function(q) {
                                return(fn(q))
                            }

                            # Run the test directly
                            test_result <- ks.test(
                                x_data, dist_fn,
                                alternative = alternative,
                                exact = exact,
                                simulate.p.value = simulate.p_value,
                                B = B
                            )

                            # Format the output
                            tibble(
                                test_type = "One-Sample",
                                null_hyp = glue::glue("The sample `{var_name}` follows the theoretical distribution specified by CDF."),
                                columns = var_name,
                                cdf_function = fn_text,  # Store the deparsed function
                                statistics = unname(test_result$statistic),
                                pval = unname(round(test_result$p.value, digits = 3))
                            )
                        }
                    )
                }) |>
                ungroup()

            return(structure(
                list(data = out, group_varnames = group_vars(data)),
                class = c("basestats", "kstest")
            ))
        }

        # Handle string column names
        if (is.character(y_expr)) {
            # User provided string literals for column names
            if (length(y_expr) == 1) {
                y_cols <- y_expr
            } else {
                y_cols <- y_expr
            }
            y_list <- syms(y_cols)
        }
        # Handle calls like c("y1", "y2") or list("y1", "y2")
        else if (is.call(y_expr) &&
                 (y_expr[[1]] == quote(c) ||
                  y_expr[[1]] == quote(list) ||
                  y_expr[[1]] == quote(cbind))) {
            # Extract character arguments if they exist
            args <- as.list(y_expr)[-1]  # Remove the function name

            # Check if all arguments are character strings
            if (all(sapply(args, is.character))) {
                # Direct character strings like c("y1", "y2")
                y_cols <- unlist(args)
                y_list <- syms(y_cols)
            } else {
                # Treat as column names like c(y1, y2)
                col_names <- all.vars(y_expr)

                # Try to select columns from data
                y_expr_str <- paste0("c(", paste(col_names, collapse = ","), ")")
                y_cols <- try(
                    tidyselect::eval_select(rlang::parse_expr(y_expr_str), data),
                    silent = TRUE
                )

                if (inherits(y_cols, "try-error")) {
                    # Handle the case where column selection fails
                    stop("Invalid column names in y: ", paste(col_names, collapse = ", "))
                }

                y_list <- syms(names(y_cols))
            }
        }
        # Handle bare column names like y = y1
        else if (is.name(y_expr)) {
            col_name <- as.character(y_expr)

            # Check if this column exists in the data
            if (col_name %in% names(data)) {
                y_list <- syms(col_name)
            } else {
                # Try to find the column by evaluating the expression in the environment
                y_value <- try(rlang::eval_tidy(y_quo), silent = TRUE)

                if (!inherits(y_value, "try-error") && is.character(y_value)) {
                    # It's a variable containing character strings of column names
                    y_list <- syms(y_value)
                } else {
                    stop("Column '", col_name, "' not found in data")
                }
            }
        }
        # Handle formulas for custom functions
        else if (rlang::is_formula(y_expr)) {
            y_expr <- rlang::new_formula(
                lhs = NULL,
                rhs = y_expr[[2]],
                env = y_env
            )
            y_list <- list(rlang::as_function(y_expr))  # Convert formula to function
        }
        else {
            # Try direct column selection
            y_cols <- try(
                tidyselect::eval_select(y_quo, data),
                silent = TRUE
            )

            if (!inherits(y_cols, "try-error")) {
                y_list <- syms(names(y_cols))
            } else {
                # Last resort: evaluate and check if result is usable
                y_value <- try(rlang::eval_tidy(y_quo), silent = TRUE)

                if (!inherits(y_value, "try-error")) {
                    if (is.character(y_value)) {
                        # Character vector of column names
                        y_list <- syms(y_value)
                    } else {
                        # Some other value - assume it's a function or reference
                        y_list <- list(y_value)
                    }
                } else {
                    stop("Cannot interpret y argument: ", deparse(y_expr))
                }
            }
        }
    } else {
        y_list <- NULL
    }

    # For column tests (two-sample test)
    x_vars <- vars
    if (!is.null(y_list) && !is.function(y_list[[1]]) && !is.list(y_list[[1]])) {
        y_vars <- y_list

        # Create pairing table
        pairs <- tidyr::expand_grid(
            x_idx = seq_along(x_vars),
            y_idx = seq_along(y_vars)
        )

        out <-
            data |>
            group_modify( ~ {
                purrr::map_dfr(
                    seq_len(nrow(pairs)), function(pair_idx) {
                        x_idx <- pairs$x_idx[pair_idx]
                        y_idx <- pairs$y_idx[pair_idx]

                        .var <- x_vars[[x_idx]]
                        y_val <- y_vars[[y_idx]]

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

                        null_hyp <- glue::glue("The sample `{deparse(.var)}` follows the same distribution as sample `{deparse(y_val)}`.")
                        test_type <- "Two-Sample"

                        test <- rlang::eval_tidy(test, data = .x)

                        tibble(
                            test_type = test_type,
                            null_hyp = null_hyp,
                            columns = deparse(.var),
                            y_column = deparse(y_val),
                            statistics = unname(test$statistic),
                            pval = unname(round(test$p.value, digits = 3))
                        )
                    }
                )
            }) |>
            ungroup()
    } else {
        # No y provided or y is a function - default to one-sample test
        out <-
            data |>
            group_modify( ~ {
                purrr::map_dfr(
                    seq_along(vars), function(idx) {
                        .var <- vars[[idx]]

                        if (rlang::quo_is_null(y_quo)) {
                            # Default to normal distribution test
                            test <- expr(
                                ks.test(
                                    !!.var,
                                    pnorm,
                                    alternative = {{alternative}},
                                    exact = {{exact}},
                                    simulate.p.value = {{simulate.p_value}},
                                    B = {{B}}
                                )
                            )

                            null_hyp <- glue::glue("`{deparse(.var)}` ~ N(0,1)")
                            test_type <- "One-Sample"
                        } else {
                            # y is a function
                            y_val <- y_list[[1]]

                            test <- expr(
                                ks.test(
                                    !!.var,
                                    y = !!y_val,
                                    alternative = {{alternative}},
                                    exact = {{exact}},
                                    simulate.p.value = {{simulate.p_value}},
                                    B = {{B}}
                                )
                            )

                            null_hyp <- glue::glue("The sample `{deparse(.var)}` follows the theoretical distribution specified by CDF.")
                            test_type <- "One-Sample"
                        }

                        test <- rlang::eval_tidy(test, data = .x)

                        tibble(
                            test_type = test_type,
                            null_hyp = null_hyp,
                            columns = deparse(.var),
                            statistics = unname(test$statistic),
                            pval = unname(round(test$p.value, digits = 3))
                        )
                    }
                )
            }) |>
            ungroup()
    }

    structure(
        list(data = out, group_varnames = group_vars(data)),
        class = c("basestats", "kstest")
    )
}

print.kstest <- function(x, ...) {
    data <- x$data
    test_type <- unique(data$test_type)

    if (test_type == "One-Sample") {
        title <- "One-Sample Kolmogorov-Smirnov Test"
    } else if (test_type == "Two-Sample") {
        title <- "Two-Sample Kolmogorov-Smirnov Test"
    }

    cat("\n")
    cli$rule(center = title) |>
        cli$style_bold() |>
        cli$cat_line("\n\n")

    data |>
        purrr::pmap(list) |>
        purrr::iwalk(function (res, i) {
            var1 <- res$columns

            cat(glue::glue("{i}. {var1}"), "\n\n")

            null_text <- paste("H0:", res$null_hyp)

            cli$rule(center = null_text, line = " ") |>
                cli$cat_line()

            # Display CDF function if available
            cat("\n")
            if (!is.null(res$cdf_function)) {
                cdf_func <- paste("CDF Function:", res$cdf_function)
                cli$rule(center = cdf_func, line = " ") |>
                    cli$cat_line()
            } else if (!is.null(res$y_column)) {
                y_col <- paste("Compared with:", res$y_column)
                cli$rule(center = y_col, line = " ") |>
                    cli$cat_line()
            }

            result_tbl <- tibble::tibble(
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
