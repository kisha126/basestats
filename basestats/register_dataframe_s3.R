box::use(
    tibble[tibble, as_tibble],
    data.table[data.table, as.data.table]
)

#' @export
as.data.frame.basestats <- function(x, ...) data.frame(x$out, ...)

#' @export
as_tibble.basestats <- function(x, ...) tibble(x$out, ...)

#' @export
as.data.table.basestats <- function(x, ...) data.table(x$out, ...)

# Register S3 methods
box::register_S3_method("as.data.frame", "basestats")
box::register_S3_method("as_tibble", "basestats")
box::register_S3_method("as.data.table", "basestats")

# box::register_S3_method("as.data.frame", "basestats", as.data.frame.basestats)
# box::register_S3_method("tibble::as_tibble", "basestats", as_tibble.basestats)
# box::register_S3_method("as.data.table", "basestats", as.data.table.basestats)

