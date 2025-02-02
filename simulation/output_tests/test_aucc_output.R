# Getting latest commit from master branch
latest_master_commit <- system(
    "git log -n 1 origin/master | grep commit | cut -d ' ' -f2-",
    intern = TRUE
)

current_commit <- system(
    "git log -n 1 | grep commit | cut -d ' ' -f2-",
    intern = TRUE
)

if (latest_master_commit != current_commit) {
    stop(
        paste0(
            "You do not have the latest commit from the master branch ",
            "on your repo. Try 'git checkout master' and 'git pull'."
        ),
        fill = TRUE
    )
}

test_completeness <- function(df, type) {
    is_df <- is.data.frame(df)
    correct_names <- all(names(df) == c(type, "method"))
    correct_types <- all(
        vapply(df, typeof, character(1L)) == c("double", "character")
    )
    no_nas <- all(is.finite(df[[type]]))
    c(
        "is_df" = is_df,
        "correct_names" = correct_names,
        "correct_types" = correct_types
    )
}

files <- c(
    "sampleSize_50-effect_0.2-I2_0.9-k_50-heterogeneity_additive-dist_Gaussian-bias_none-large_1.rds",
    "sampleSize_50-effect_0.2-I2_0.9-k_50-heterogeneity_additive-dist_Gaussian-bias_none-large_2.rds",
    "sampleSize_50-effect_0.2-I2_0.9-k_50-heterogeneity_additive-dist_snl-bias_none-large_0.rds",
    "sampleSize_50-effect_0.2-I2_0.9-k_50-heterogeneity_additive-dist_snl-bias_none-large_1.rds",
    "sampleSize_50-effect_0.2-I2_0.9-k_50-heterogeneity_additive-dist_snl-bias_none-large_2.rds"
)

res <- sapply(
    files,
    function(file) {
        dat <- readRDS(file = file.path("..", "RData", "CIs", file))
        aucc <- lapply(dat$cis, "[[", i = "aucc")
        aucc_ratio <- lapply(dat$cis, "[[", i = "aucc_ratio")
        aucc_tests <- sapply(aucc, test_completeness, type = "aucc")
        # idx <- which(!aucc_tests, arr.ind = TRUE)[, 2L]
        # aucc[idx]
        aucc_ratio_tests <- sapply(
            aucc_ratio,
            test_completeness,
            type = "aucc_ratio"
        )
        aucc_summary <- apply(aucc_tests, 1, all)
        aucc_ratio_summary <- apply(aucc_ratio_tests, 1, all)
    }
)

if (!all(res)) {
    err_idx <- which(res == FALSE, arr.ind = TRUE)
    dn <- dimnames(res)
    err_files <- unique(err_idx[, 2L])
    errors <- vapply(
        err_files,
        function(x) {
            paste0(names(err_idx[, 1L])[err_idx[, 2L] == x], collapse = ", ")
        },
        character(1L)
    )
    file_names <- dn[[2]][err_idx[, 2L]]
    msg <- paste0(
        paste0(file_names, ":\n", errors),
        collapse = "\n\n"
    )
    cat(
        paste0(
            "There have been the following errors:\n\n",
            msg
        ),
        fill = TRUE
    )
} else {
    cat("No unexpected errors found", fill = TRUE)
}
