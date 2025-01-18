# get a vector of script files
script_dir <- "."
Rscripts <- paste0(script_dir, "/", list.files(script_dir))
Rscripts <- Rscripts[grepl(".R$", Rscripts)]
Rscripts <- Rscripts[!grepl("plot", Rscripts)] # remove plotting

extract_all <- function(strings, pattern) {
    lapply(strings, function(string) {
        # Use gregexpr to find all matches and regmatches to extract them
        matches <- gregexpr(pattern, string, perl = TRUE)
        regmatches_list <- regmatches(string, matches)
        # Extract only the groups (inside parentheses) for each match
        unname(sapply(regmatches_list[[1]], function(match) {
            sub(pattern, "\\1", match, perl = TRUE)
        }))
    })
}

# extract packages used from all R scripts
libraries <- do.call(
    "c",
    lapply(Rscripts, function(x) {
        # read the scripts
        code <- readLines(x, warn = FALSE)
        # throw out comments
        code <- gsub("#.*$", "", code)
        # get all the calls to library or require
        lib_pat <- "\\s*library\\s*\\((.+?)\\)"
        req_pat <- "require\\s*\\((.+?)\\)"
        lib_lines <- grep(lib_pat, code, value = TRUE)
        req_lines <- grep(req_pat, code, value = TRUE)
        # Extract the libraries
        libs <- do.call("c", extract_all(lib_lines, lib_pat))
        reqs <- do.call("c", extract_all(req_lines, req_pat))
        # Get namespaced calls
        ns_calls <- grep(":{2,3}", code, value = TRUE)
        ns <- lapply(
            ns_calls,
            extract_matches,
            pattern = "([[:alnum:]]+?):{2,3}"
        )
        ns <- unique(sub(":{2,3}", "", do.call("c", ns)))
        # return all libraries
        pkgs <- unique(c(libs, reqs, ns))
        pkgs
    })
)

# find installed/uninstalled packages
installed <- vapply(libraries, function(x) {
    out <- tryCatch(
        {
            find.package(x)
            TRUE
        },
        warning = function(cond) {
            stop("A warning appeared when installing packages.\nCheck install_packages.R\n")
        },
        error = function(cond) {
            return(FALSE)
        }
    )
    return(out)
}, logical(1L), USE.NAMES = TRUE)

# Get uninstalled packages. Always add confMeta such that the newest version is
# installed.
install_libs <- unique(c(libraries[!installed], "confMeta"))

# install missing packages
if (length(install_libs) == 0) {
    cat("All necessary libraries already installed.\n")
} else {
    cat(
        paste0(
            "Installing packages:\n",
            paste(install_libs, collapse = "\n"),
            "\n"
        )
    )
    install <- lapply(install_libs, function(x) {
        out <- tryCatch(
            {
                install.packages(
                    x,
                    repos = "https://cloud.r-project.org/",
                    quiet = TRUE
                )
                list(status = 0L, message = NA_character_)
            },
            warning = function(w) {
                message <- paste0(conditionMessage(w), "\n")
                list(status = 1L, message = message)
            },
            error = function(e) {
                message <- paste0(conditionMessage(w), "\n")
                list(status = 2L, message = message)
            }
        )
        # If package installation threw a warning, print warning message
        if (out$status == 1L) {
            cat(
                paste0(
                    "Installation of package '", x,
                    "' threw the following warning:\n\n",
                    out$message, "\n"
                )
            )
            inst <- tryCatch(
                {
                    remotes::install_github(
                        paste0("felix-hof/", x),
                        ref = "dev"
                    )
                },
                error = function(e) {
                    msg <- conditionMessage(e)
                    cat(
                        paste0(
                            "Error during package installation: ",
                            x, "\n"
                        )
                    )
                }
            )
        }
        # If package installation threw an error, i.e. package not on CRAN, then
        # try to install from my github
        if (out$status == 2L) {
            inst <- tryCatch(
                {
                    remotes::install_github(
                        paste0("felix-hof/", x, ref = "dev")
                    )
                },
                error = function(e) {
                    msg <- conditionMessage(e)
                    cat(
                        paste0(
                            "Error during package installation:"
                        )
                    )
                }
            )
        }
        invisible(NULL)
    })
}
