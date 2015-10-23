# methods for standard generic functions: show

# ------------------
# Class definition
# ------------------

mesh <- setClass("mesh", contains="numeric")

setValidity("mesh", function(object) {
  msg <- NULL
  valid <- TRUE
  if (length(object)<2) {
    valid <- FALSE
    msg <- c(msg,
             "Mesh should have at least two knot points")
  } else {
    if (min(diff(object))<=0) {
      valid <- FALSE
      msg <- c(msg,
               "Increments of mesh segments must be strictly positive")
    }
  }
  if (valid) TRUE else msg
})

# ---------
# Methods
# ---------

setMethod("show",
          signature="mesh",
          definition=function(object) {
            p <- length(object)
            cat("An object of class ",class(object),"\n",sep="")
            cat(" Partition of [",object[1],",",object[p],"] in ",p-1," segments.",
                " 5-number statistics of segment lengths:\n",sep="")
            print(summary(diff(object)))
            invisible(NULL)
          })
