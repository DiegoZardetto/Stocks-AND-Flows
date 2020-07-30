mem <- function() {
  bit <- 8L * .Machine$sizeof.pointer
  if (!(bit == 32L || bit == 64L)) {
    stop("Unknown architecture", call. = FALSE)
  }
  
  node_size <- if (bit == 32L) 28L else 56L
  
  usage <- gc()
  sum(usage[, 1] * c(node_size, 8)) / (1024 ^ 2)
}

mem.change <- function(code) {
  start <- mem()
  
  expr <- substitute(code)
  eval(expr, parent.frame())
  rm(code, expr)
  
  round(mem() - start, 3)
}