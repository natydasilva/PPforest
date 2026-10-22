test_that("mvote counts votes for a label when a lower label gets none", {
  votes <- matrix(c(1, 3, 3), nrow = 3)
  expect_equal(as.numeric(PPforest:::mvote(votes)), 3)
})

test_that("mvote picks the majority for each observation", {
  votes <- cbind(c(2, 2, 1), c(4, 1, 4), c(3, 3, 3))
  expect_equal(as.numeric(PPforest:::mvote(votes)), c(2, 4, 3))
})

test_that("mvoteoob keeps a column for every group and counts only out-of-bag votes", {
  votes <- cbind(c(1, 3, 3), c(2, 2, 3))
  oob <- cbind(c(1, 1, 1), c(1, 1, 0))
  result <- PPforest:::mvoteoob(votes, oob, 4)
  expect_equal(dim(result), c(2, 5))
  expect_equal(result[1, ], c(1, 0, 2, 0, 3))
  expect_equal(result[2, ], c(0, 2, 0, 0, 2))
})

test_that("PPforest trains with many overlapping groups", {
  # 32 normal groups with random means and covariances, which overlap enough
  # that some group gets no vote from any tree.
  set.seed(0)
  q <- 32
  p <- 4
  groups <- lapply(seq_len(q), function(i) {
    root <- matrix(stats::rnorm(p^2, 0, 25), nrow = p)
    x <- matrix(stats::rnorm(64 * p), ncol = p) %*% root
    sweep(x, 2, stats::rnorm(p, 100 * i, 25), "+")
  })
  data <- data.frame(do.call(rbind, groups), Type = factor(rep(seq_len(q), each = 64)))

  forest <- PPforest(data = data, y = "Type", std = "no", size.tr = 1, m = 20,
                     PPmethod = "PDA", lambda = 0.5, size.p = 1)
  expect_equal(ncol(forest$votes), q)
  expect_true(forest$oob.error.forest >= 0 && forest$oob.error.forest <= 1)
})

test_that("a tree predicts the same group codes for rows alone or in a larger batch", {
  set.seed(0)
  forest <- PPforest(data = iris, y = "Species", std = "no", size.tr = 2/3, m = 20,
                     PPmethod = "PDA", lambda = 0.5, size.p = 0.5)
  virginica <- iris[iris$Species == "virginica", 1:4][1:3, ]
  mixed <- rbind(iris[c(1, 51), 1:4], virginica)
  for (tree in forest$output.trees) {
    alone <- as.numeric(PPforest:::PPclassify(tree, virginica)[[2]])
    in_batch <- as.numeric(PPforest:::PPclassify(tree, mixed)[[2]])[3:5]
    expect_equal(alone, in_batch)
  }
})

test_that("predict() gives a row the same label alone or in a larger batch", {
  set.seed(0)
  forest <- PPforest(data = iris, y = "Species", std = "no", size.tr = 2/3, m = 50,
                     PPmethod = "PDA", lambda = 0.5, size.p = 0.5)
  rows <- iris[c(101, 1, 51), 1:4]
  together <- as.character(predict(forest, newdata = rows, parallel = FALSE)[[3]])
  alone <- vapply(seq_len(nrow(rows)), function(i) {
    as.character(predict(forest, newdata = rows[i, ], parallel = FALSE)[[3]])
  }, character(1))
  expect_equal(alone, together)
  expect_equal(together, c("virginica", "setosa", "versicolor"))
})

test_that("the out-of-bag error is the share of wrong out-of-bag predictions", {
  set.seed(0)
  forest <- PPforest(data = iris, y = "Species", std = "no", size.tr = 2/3, m = 50,
                     PPmethod = "PDA", lambda = 0.5, size.p = 0.5)
  expect_equal(forest$oob.error.forest,
               mean(as.character(forest$prediction.oob) != as.character(forest$train$Species)))
})

test_that("predict() handles a single row", {
  set.seed(0)
  forest <- PPforest(data = iris, y = "Species", std = "no", size.tr = 2/3, m = 20,
                     PPmethod = "PDA", lambda = 0.5, size.p = 0.5)
  expect_equal(as.character(predict(forest, newdata = iris[101, 1:4], parallel = FALSE)[[3]]), "virginica")
})
