box::use(
    tibble[tibble, rownames_to_column],
    dplyr[select, summarise, group_vars, rename, mutate, filter, pick, group_split, n],
    rlang[expr, quo, enquo],
    tidyselect[eval_select, everything, all_of],
    tidyr[pivot_longer, unnest, expand_grid],
    purrr[pluck, iwalk, pmap, map2_dbl, map_dbl],
    stats[kruskal.test, qt, pnorm],
    cli,
    tabhelpers[table_default]
)

#' @export
kw_test1 <- function(.data, columns, ...) {
    columns_quo <- enquo(columns)
    new_data <- .data[, tidyselect::eval_select(columns_quo, data = .data), drop = FALSE]

    out <- new_data |>
        summarise(
            kw_model = list({
                data_list <- as.list(pick(everything()))
                test <- kruskal.test(data_list)

                tibble(
                    statistic = test$statistic,
                    df = test$parameter,
                    pval = test$p.value
                )
            })
        ) |>
        unnest(kw_model)

    res <- list(
        out = out,
        grouping_var = group_vars(.data)
    )

    class(res) <- c("basestats", "kwtest")
    res
}

#' @export
kw_test2 <- function(.data, x, g, ...) {
    new_data <- select(.data, x, g)

    out <- new_data |>
        summarise(
            kw_model = list({
                data_list <- pick(everything())
                test <- kruskal.test(formula = formula, data = data_list)

                tibble(
                    statistic = test$statistic,
                    df = test$parameter,
                    pval = test$p.value
                )
            })
        ) |>
        unnest(kw_model)

    res <- list(
        out = out,
        grouping_var = group_vars(.data)
    )

    class(res) <- c("basestats", "kwtest")
    res
}

#' @export
kw_test3 <- function(.data, formula, ...) {
    x <- rlang::f_lhs(formula)
    g <- rlang::f_rhs(formula)

    new_data <- select(.data, x, g)

    out <- new_data |>
        summarise(
            kw_model = list({
                data_list <- pick(everything())
                test <- kruskal.test(formula = formula, data = data_list)

                tibble(
                    statistic = test$statistic,
                    df = test$parameter,
                    pval = test$p.value
                )
            })
        ) |>
        unnest(kw_model)

    res <- list(
        out = out,
        grouping_var = group_vars(.data)
    )

    class(res) <- c("basestats", "kwtest")
    res
}

#' @export
pairwise_kwtest <- function(x, ...) UseMethod("pairwise_kwtest")

pairwise_kwtest.default <- function(x, g, data = NULL, ...) {
    if (!is.null(data)) {
        new_data <- select(data, value = x, group = g)
    } else {
        new_data <- data.frame(value = x, group = g)
    }

    out <- new_data |>
      mutate(
        r = rank(value),
        stat = kruskal.test(x = x, g = g)$statistic
      ) |>
      summarise(
        mean_r = mean(r),
        n = n(),
        stat = unique(stat),
        .by = group
      ) |>
      summarise(
        adhoc = list({
          groups <- group
          names(mean_r) <- groups
          N <- sum(n)
          n <- n
          names(n) <- groups
          k <- length(groups)
          T <- unique(stat)
          S2 <- N * (N + 1) / 12
          T1 <- (S2 * (N - 1 - T)) / (N - k)
          t_crit <- qt(1 - 0.05 / 2, N - k)

          # --- Pair names ---
          pairs <- expand_grid(x = groups, y = groups) |>
            filter(x < y)

          pair <- paste(
            gsub("x", "", pairs$x),
            "and",
            gsub("x", "", pairs$y)
          )

          # --- Mean Rank Difference ---
          diff <- map2_dbl(
            pairs$x,
            pairs$y,
            ~ (mean_r[.x] - mean_r[.y])
          )

          # --- Standard Error of Mean Rank Diff ---
          std_diff <- map2_dbl(
            pairs$x,
            pairs$y,
            ~ {
              t_crit * sqrt(T1) * sqrt(1 / n[.x] + 1 / n[.y])
            }
          )

          # --- T-statistic and Corresponding p-value ---
          test_stat <- diff / std_diff
          pval <- map_dbl(
            test_stat, ~ {
              p1 <- pnorm(.x, lower.tail = T)
              p2 <- pnorm(.x, lower.tail = F)
              min(p1, p2) * 2
            }
          )

          tibble(
            Methods = pair,
            diff = diff,
            std_diff = std_diff,
            test_stat = test_stat,
            pval = pval
          )
        })
      ) |>
      unnest(adhoc)

    res <- list(
        out = out,
        pairwise = TRUE
    )
    class(res) <- c("basestats", "kwtest")
    res
}

pairwise_kwtest.kwtest <- function(x, ...) {
    new_data <- x$out
}


print.kwtest <- function(x, ...) {
    data <- x$out |>
        rename(
            `Chi-Square statistic` = statistic
        )

    grouping_var <- x$grouping_var
    title <- "Kruskal-Wallis H Test"

    # has_group <- length(grouping_var) > 0
    cat("\n\n")
    table_default(
        x = data,
        title = title,
        n_space = 3,
        style_columns = list(
            pval = function (ctx) {
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
        center_table = TRUE,
        ...
    ) |>
        suppressWarnings()

}

box::register_S3_method("print", "kwtest")



