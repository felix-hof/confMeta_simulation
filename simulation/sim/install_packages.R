# get a vector of script files
script_dir <- "."
Rscripts <- paste0(script_dir, "/", list.files(script_dir))
Rscripts <- Rscripts[grepl(".R$", Rscripts)]
Rscripts <- Rscripts[!grepl("plot", Rscripts)] # remove plotting

# Function for extracting multiple matches
extract_matches <- function(data, pattern) {
    start <-  gregexpr(pattern, data)[[1]]
    stop  <-  start + attr(start, "match.length") - 1
    if(-1 %in% start) {
	""    ## **Note** you could return NULL if there are no matches 
    } else {
	mapply(substr, start, stop, MoreArgs = list(x = data))
    }
}

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
    # get all libraries
    libs <- c(libs, reqs)
    # get this into nice format
    libs <- unlist(strsplit(libs, ","), recursive = TRUE)
    libs <- gsub("\\W", "", libs)
    # Get namespaced calls
    ns_calls <- grep(":{2,3}", code, value = TRUE)
    ns <- lapply(ns_calls, extract_matches, pattern = "([[:alnum:]]+?):{2,3}")
    ns <- unique(sub(":{2,3}", "", do.call("c", ns)))
    # return all libraries
    unique(c(libs, ns))
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
install_libs <- libraries[!installed]
if(length(install_libs) == 0){
    cat("All necessary libraries already installed.\n")
} else {
    cat("Installing packages:\n", paste(install_libs, collapse = "\n", "\n"))
    lapply(install_libs, function(x){
	out <- tryCatch({install.packages(x, repos = "https://cloud.r-project.org/"); 0L},
	    warning = function(w) 1L,
	    error = function(e) 2L)
	if(out %in% 2L) remotes::install_github(paste0("felix-hof/", x))
	invisible(NULL)
    })
}

