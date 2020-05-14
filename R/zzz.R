##
## Load  metabolitesMapping into namespace
##

.onLoad <- function(libname, pkgname) {
    
    ns <- asNamespace( pkgname)
    ah <- AnnotationHub::AnnotationHub()
    metabolitesMapping <- ah[["AH79817"]]
        
    assign( "metabolitesMapping", metabolitesMapping, envir = ns)
    
    namespaceExport( ns, "metabolitesMapping")
    
}