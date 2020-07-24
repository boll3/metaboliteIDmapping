##
## Load  metabolitesMapping into namespace
##

.onLoad <- function(libname, pkgname) {
    
    ns <- asNamespace( pkgname)
    ah <- AnnotationHub::AnnotationHub()
    metabolitesMapping <- ah[["AH83115"]]
        
    assign( "metabolitesMapping", metabolitesMapping, envir = ns)
    
    namespaceExport( ns, "metabolitesMapping")
    
}