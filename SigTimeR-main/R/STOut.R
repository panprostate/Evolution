setClass("STOut",
   slots=c(
     estimate="data.frame",
     metadata="character",
     raw="list"
   )
)

STOut <- function(estimate, metadata, raw){
  new("STOut", estimate=estimate, metadata=metadata, raw=raw)
}

setGeneric("stmethod", function(x) standardGeneric("stmethod"))
setMethod("stmethod", "STOut", function(x){
  x@estimate$method[[1]]
})

setMethod("show", "STOut",
  function(object){
    cat(class(object), "instance with", stmethod(object), "results.")
  }
)