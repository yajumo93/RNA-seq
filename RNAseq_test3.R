
# test

library(edgeR)

bcv <- 0.2
test_count <- matrix(rnbinom(40, size=1/bcv^2, mu=10), 20, 2)

test_count

test_y <- DGEList(counts = test_count, group = 1:2)

test_y

test_et <- exactTest(test_y, dispersion = bcv^2)
test_et

plotMDS(test_et)
