library(phylosignal)
context("Phylogenetic signal tests")

data(navic)

#adephylo::abouheif.moran(navic, method="Abouheif")$obs
test_that("Cmean is correct (cf. adephylo::abouheif.moran)", {
  expect_equal(as.numeric(phyloSignal(navic, methods="Cmean")$stat), 0.4791519, tol = 1e-3)
})

#adephylo::abouheif.moran(navic, method="patristic")$obs
test_that("I is correct (cf. adephylo::abouheif.moran)", {
  expect_equal(as.numeric(phyloSignal(navic, methods="I")$stat), 0.0428604, tol = 1e-3)
})

# picante::phylosignal(x = tipData(navic)$IPSS, phy = as(navic, "phylo"))$K
test_that("K is correct (cf. picante::phylosignal)", {
  expect_equal(as.numeric(phyloSignal(navic, methods="K")$stat), 0.7897245, tol = 1e-3)
})

# phytools::phylosig(tree = as(navic, "phylo"), x = tipData(navic)$IPSS, method = "lambda")$lambda
test_that("Lambda is correct (cf. phytools::phylosig)", {
  expect_equal(as.numeric(phyloSignal(navic, methods="Lambda")$stat), 0.9588398, tol = 1e-3)
})

