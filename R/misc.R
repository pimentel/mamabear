msg <- function(..., nl = TRUE) {
  message(..., appendLF = nl)
}

# shortcut for as.data.frame(x, stringsAsFactors = FALSE)
as_df <- function(x, ...) {
  as.data.frame(x, stringsAsFactors = FALSE, ...)
}

adf <- function(...) {
  data.frame(..., stringsAsFactors = FALSE)
}

dot <- function(count, max_dots = 50) {
  count <- count + 1
  new_line <- FALSE
  if ( (count %% max_dots) == 0 ) {
    new_line <- TRUE
  }

  msg('.', nl = new_line)

  count
}
