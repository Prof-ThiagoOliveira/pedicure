.onLoad <- function(libname, pkgname) {
  # don't load libs for man
  if(Sys.getenv("load_libs", unset="") == "" || Sys.getenv("load_libs") != "no") {
    if(Sys.info()["sysname"] == "Windows") {
      # library.dynam("api-ms-win-crt-convert-l1-1-0", pkgname, libname)
      # library.dynam("api-ms-win-crt-heap-l1-1-0", pkgname, libname)
      # library.dynam("api-ms-win-crt-runtime-l1-1-0", pkgname, libname)
      # library.dynam("api-ms-win-crt-stdio-l1-1-0", pkgname, libname)
      # library.dynam("api-ms-win-crt-string-l1-1-0", pkgname, libname)
      library.dynam("libiomp5md", pkgname, libname)
      library.dynam("vcruntime140", pkgname, libname)
      library.dynam("vcomp140", pkgname, libname)
      library.dynam("mkl_pedicure", pkgname, libname)
    }
    library.dynam("pedicure", pkgname, libname)
  }
}
