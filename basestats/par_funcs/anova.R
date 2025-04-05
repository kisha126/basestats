box::use(
    tibble[tibble, rownames_to_column],
    rlang[expr, quo],
    dplyr[select, summarise, group_vars, rename, filter, pick, group_split],
    tidyselect[everything, all_of],
    tidyr[pivot_longer, unnest],
    purrr[pluck, iwalk, pmap],
    stats[aov],
    cli,
    tabhelpers[table_default]
)

#' @export
anova_test <- function(.data, columns = NULL, formula = NULL, ...) {
    # Check if using formula approach
    if (is.null(columns) && !is.null(formula)) {
        out <- aov(formula = formula, data = .data, ...) |>
            broom::tidy()

        res <- list(
            out = out,
            grouping_var = group_vars(.data)
        )

        class(res) <- c("basestats", "ANOVA")
        return(res)
    }

    # For column-based approach
    new_data <- select(
        .data,
        all_of(tidyselect::eval_select(quo(columns), data = .data))
    )

    # Only worked with One-Way Anova
    out <- new_data |>
        summarise(
            aov_model = list({
                selected_cols <- pick(everything()) |>
                    select(where(is.numeric)) |>
                    pivot_longer(
                        cols = everything(),
                        names_to = "g",
                        values_to = "y"
                    )

                model <- aov(y ~ g, data = selected_cols, ...)

                broom::tidy(model)
            })
        ) |>
        unnest(aov_model)

    res <- list(
        out = out,
        grouping_var = group_vars(.data)
    )

    class(res) <- c("basestats", "ANOVA")
    res
}

print.ANOVA <- function(x, ...) {
    data <- x$out |>
        rename(
            ` ` = term,
            `F-value` = statistic,
            `Pr(>F)` = p.value
        )

    grouping_var <- x$grouping_var
    title <- "ANOVA Table"

    has_group <- length(grouping_var) > 0
    cat("\n\n")
    if (!has_group) {
        table_default(
            x = data,
            title = title,
            justify_cols = list(
                "1" = "right",
                "2" = "right",
                "3" = "right",
                "4" = "right",
                "5" = "right"
            ),
            vb = list(after = 1),
            n_space = 3,
            style_columns = list(
                `Pr(>F)` = function (ctx) {
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
    } else {
        group_names <- unique(data[[grouping_var]])

        data |>
            split(data[[grouping_var]]) |>
            unname() |>
            iwalk(
                function (res, i) {
                    gn <- group_names[i]
                    starter_text <- glue::glue("{i}. {gn}")
                    cat(starter_text, "\n\n")

                    table_default(
                        x = select(res, -1),
                        title = title,
                        justify_cols = list(
                            "1" = "right",
                            "2" = "right",
                            "3" = "right",
                            "4" = "right",
                            "5" = "right"
                        ),
                        vb = list(after = 1),
                        n_space = 3,
                        style_columns = list(
                            `Pr(>F)` = function (ctx) {
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
                    ) |> suppressWarnings()

                }
            )
    }
}

box::register_S3_method("print", "ANOVA")
