box::use(
    tibble[tibble, as_tibble]
)

as.data.frame.basestats <- function(x, ...) data.frame(x$out, ...)
as_tibble.basestats <- function(x, ...) tibble(x$out, ...)

box::register_S3_method('as.data.frame', 'basestats')
box::register_S3_method('as_tibble', 'basestats')
