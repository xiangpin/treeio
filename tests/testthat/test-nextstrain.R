context("nextstrain input")

library(treeio)

file1 <- system.file("extdata/nextstrain.json", "minimal_v2.json", package = "treeio")

test_that("read.nextstrain.json returns a valid treedata object", {
    expect_silent(tr <- read.nextstrain.json(file1))
    expect_s4_class(tr, "treedata")
    expect_true(ape::is.rooted(tr@phylo))
    expect_equal(ape::Ntip(tr@phylo), 31L)
    expect_equal(ape::Nnode(tr@phylo), 25L)
    expect_equal(nrow(tr@data), 56L)
    expect_true(any(tr@phylo$tip.label == "Thailand/1610acTw"))
    expect_true("country" %in% names(tr@data))
})
