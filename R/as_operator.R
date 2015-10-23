as.operator <- function(e1) {
  operator(alpha=e1@alpha,
           beta=e1@epsilon|-e1@delta,
           gamma=e1@delta|e1@epsilon,
           delta=e1@delta,
           epsilon=e1@epsilon)
}
