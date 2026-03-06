
# Functions sent by Wolfram Grondwald (WG) on Tuesday, 3. March 2026 at 13:12,for reproducing the results from
# 'Analysis of human urine reveals metabolic changes related to the development
# of acute kidney injury following cardiac surgery' by Zacharias et al. (2012).

# Attempt a classification for Daenarys project

install_deps <- function() {
    install.packages("e1071")
    library(e1071)
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }

    BiocManager::install()
    # now install multtest for multiple testing adjustments
    BiocManager::install("multtest")
    library(multtest)

    # remove all old data
    rm(list = ls())
}

# load required functions 'fsel' and 'CrossVal'
# fsel
fsel <- function(x, y, n) # x data, y labels, n how many features will be used for predictions
{
    dat <- t(x) # transpose back to normal order only for SVM
    labels <- as.integer(y) - 1 # labels
    tscores <- multtest::mt.teststat(dat, labels, test = "t") # tscores
    sel <- order(-abs(tscores))[1:n] # select best
    sel # return selection
}

# CrossVal
# x = data
# y = labels
# n = number of selected features e.g. 15
CrossVal <- function(x, y, n)
{
    tot_acc <- as.numeric(10)
    k <- nrow(x) / 2 # number of test sets equals half number of patients
    sel_save <- matrix(1:10 * k * n, 10 * k, n)
    inpred <- numeric(nrow(x)) # number of patients
    inpredno <- numeric(nrow(x)) # number of patients
    m <- 0
    for (j in 1:10) # select 10 times different buckets
    {
        perm <- sample(seq_len(nrow(x))) # Create permutations from samples
        heaps <- matrix(perm, nrow = 2) # Split in sets a 2 samples
        k <- nrow(x) / 2 # number of sets
        acc <- numeric(k)
        for (i in 1:k) # over all test sets
        {
            test <- heaps[, i] # select 2 samples for testing
            m <- m + 1
            sel <- fsel(x[-test, ], y[-test], n) # exclude selected samples
            sel_save[m, ] <- sel
            SVM <- e1071::svm(x[-test, sel], y[-test], kernel = "linear", cost = 0.5) # train SVM
            p <- predict(SVM, x[test, sel]) # predict excluded samples
            acc[i] <- 100 * sum(p == y[test]) / length(test)
            tmp1 <- test[1]
            tmp2 <- test[2]
            inpredno[tmp1] <- inpredno[tmp1] + 1 # counts how often one patient was tested
            inpredno[tmp2] <- inpredno[tmp2] + 1 # counts how often one patient was tested
            inpred[tmp1] <- inpred[tmp1] + as.numeric(p[1]) - 1 # Sum of predictions
            inpred[tmp2] <- inpred[tmp2] + as.numeric(p[2]) - 1 # Sum of predictions
        }
        tot_acc[j] <- mean(acc)
        cat("repeat loop", j, "total accuracy", tot_acc[j], "\n")
        print(round(acc))
    }
    for (b in seq_len(nrow(x))) # Give individual predictions
    {
        cat(
            "ind pred patient", b, "no pred", inpredno[b], "sum of pred",
            inpred[b], "Ratio", inpred[b] / inpredno[b], "\n"
        )
    }
    list(avg.tot.acc = mean(tot_acc))
}

read_bt3 <- function() {
    # read urine bucket table, data already normalized relative to creatinine we do PQN in addition.
    bt <- read.table(pkg_file("akiwg/bucket_table"), as.is = TRUE, header = TRUE, row.names = 1)
    bt1 <- as.matrix(t(bt))
    kick.out <- apply(bt1, 1, function(z) {
        all(z == 0)
    })
    bt2 <- bt1[!kick.out, ]
    # exclude patients 47,54, 83 as diagnosis is unclear, Note Patient 42 already excluded as no sample available
    bt3 <- bt2[, -c(46, 53, 82)] # 106 samples with 701 buckets
    bt3
}

# Main function
wolframs_cv <- function(pqn = FALSE)
{
    bt3 <- read_bt3()

    # PQN normalization
    reference <- apply(bt3, 1, median)
    quotient <- bt3 / reference
    quotient.median <- apply(quotient, 2, median)
    bt3.pqn <- t(t(bt3) / quotient.median)

    l2 <- as.factor(l1) # for classification
    bt4 <- if (pqn) t(bt3.pqn) else t(bt3) # transpose back for SVM
    res <- CrossVal(bt4, l2, 10) # starting with selection of 10 features for predictions
    res
}

l1 <- c(
    0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1,
    1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0,
    0, 0, 1, 0, 0, 1, 0, 0, 0, 0
) # labels for samples
