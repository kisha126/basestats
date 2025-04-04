box::use(
    dplyr[filter],
    tidyr[expand_grid],
    purrr[pmap]
)

inequality <- function(a, b, direction = "leq") {
    op <- switch(
        direction,
        "lt" = `<`,
        "leq" = `<=`,
        "gt" = `>`,
        "geq" = `>=`,
        "eq" = `==`,
        "neq" = `!=`,
        "all" = function(a, b) TRUE,
        stop("Invalid direction specified.")
    )
    op(a, b)
}

all_pairs <- function(x, direction = "leq", simplify = TRUE) {
    if (simplify) {
        expand_grid(x = x, y = x) |>
            filter(inequality(x, y, direction = direction))
    } else {
        expand_grid(x = x, y = x) |>
            filter(inequality(x, y, direction = direction)) |>
            pmap(~ list(.x, .y))
    }
}

format_pval <- function(pval, digits = 3, ...) {
    if (pval <= 0.001) {
        return("<0.001")
    } else if (pval <= 0.01) {
        return(format(round(pval, digits + 1), nsmall = digits + 1))
    } else {
        return(format(round(pval, digits), nsmall = digits))
    }
}
