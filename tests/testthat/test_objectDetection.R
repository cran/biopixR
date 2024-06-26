library(testthat)
library(biopixR)

test_that("objectDetection", {
  img <- beads
  res_objectDetection <- objectDetection(img, alpha = 0.75, sigma = 0.1)

  mat <- matrix(0, 4, 4)
  expect_error(objectDetection(mat))
  img_magick <- cimg2magick(img)
  expect_error(objectDetection(img_magick),
    regexp = "image must be of class 'cimg'"
  )

  expect_warning(cannyEdges(add.color(img)))
  expect_warning(objectDetection(droplet_beads),
    regexp = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S")," Image is from a luminescence channel and was converted into grayscale")
  )

  expect_equal(class(img)[1], "cimg")
  expect_equal(length(dim(img)), 4)
  expect_equal(dim(add.color(img))[4], 3)
  expect_equal(
    length(res_objectDetection$centers$value),
    length(res_objectDetection$center$size)
  )

  expect_type(res_objectDetection, "list")
  expect_length(res_objectDetection, 3)
  expect_equal(
    length(unique(res_objectDetection$coordinates$value)),
    length(res_objectDetection$centers$value)
  )
  expect_lt(nrow(res_objectDetection$centers), 15)
  expect_equal(res_objectDetection$centers$value[1], 1)
  expect_equal(
    seq_along(res_objectDetection$centers$value),
    seq_len(nrow(res_objectDetection$centers))
  )
})
