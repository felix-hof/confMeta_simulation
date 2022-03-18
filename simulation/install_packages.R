# get a vector of script files
script_dir <- "."
Rscripts <- paste0(script_dir, "/", list.files(script_dir))
Rscripts <- Rscripts[grepl(".R$", Rscripts)]

# extract packages used from all R scripts
libraries <- lapply(Rscripts, function(x){
  # read the scripts
  code <- readLines(x, warn = FALSE)
  # throw out comments
  code <- gsub("\\s*#.*", "", code)
  # throw out newlines and tabs
  code <- gsub("\\n|\\t", "", code)
  # concatenate lines
  code <- paste(code[code != ""], collapse = "")
  # insert newlines after function calls and split at newlines
  code <- unlist(strsplit(gsub("\\)", "\\)\\\n", code), "\\n"))
  # get libraries
  libs <- gsub("\\s*library\\((.+)\\)\\s*", "\\1", code[grepl("library\\(", code)])
  reqs <- gsub("\\s*require\\((.+)\\)\\s*", "\\1", code[grepl("require\\(", code)])
  libs <- c(libs, reqs)
  # get this into nice format
  libs <- unlist(strsplit(libs, ","), recursive = TRUE)
  libs <- gsub("\\W", "", libs)
})

# unlist
libraries <- unique(unlist(libraries))

# find installed/uninstalled packages
installed <- vapply(libraries, function(x){
  out <- tryCatch({
    find.package(x)
    TRUE
  },
  warning = function(cond){
    stop("A warning appeared when installing packages.\nCheck install_packages.R\n")
  },
  error = function(cond){
    return(FALSE)
  })
  return(out)
}, logical(1L), USE.NAMES = TRUE)

# install uninstalled ones CRAN packages
libs <- libraries[!installed]
if(length(libs) == 0){
	cat("All necessary libraries already installed.\n")
} else {
	cat("Installing packages:\n", paste(libs, collapse = "\n", "\n"))
	lapply(libs, function(x){
	  out <- tryCatch({install.packages(x); 0L},
		   warning = function(w) "NA on CRAN",
		   error = function(e) "NA on CRAN")
	  if(out == "NA on CRAN") devtools::install_github(paste0("felix-hof/", x))
	  return(NULL)
	})
}

