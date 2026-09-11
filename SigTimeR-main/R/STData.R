setClass("STData",
  slots=c(
    data="matrix",
    sigs="matrix"
  )
)

#'@export
STData <- function(epoch_1, epoch_2, sigs){
  data = t(cbind(epoch_1, epoch_2))
  new("STData", data=data, sigs=sigs)
}

setGeneric("epoch_sum", function(x, ...) standardGeneric("epoch_sum"))
setMethod("epoch_sum", "STData", function(x, e){
  sum(x@data[e,])
})

setGeneric("nsigs", function(x) standardGeneric("nsigs"))
setMethod("nsigs", "STData", function(x){
  dim(x@sigs)[2]
})

setGeneric("stdata", function(x) standardGeneric("stdata"))
setMethod("stdata", "STData", function(x){
  x@data
})

setGeneric("stsigs", function(x) standardGeneric("stsigs"))
setMethod("stsigs", "STData", function(x){
  x@sigs
})

setMethod("show", "STData",
    function(object){
      cat(class(object), "instance with", epoch_sum(object, 1), "in epoch_1,",
          epoch_sum(object, 2), "in epoch_2, and", nsigs(object), "signatures.")
    }
)