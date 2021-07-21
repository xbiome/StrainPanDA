expect_range <- function(obj, n) {
  # 1. Capture the object and label
  act <- quasi_label(rlang::enquo(obj), arg = "object")

  # 2. Call expert()
  expect(
    abs(act$val$rank - n) <= 1,
    sprintf("%s has value %i, not in the range(+/-1) of %i", act$lab, act$val$rank, n)
  )

  # 3. Invisibly return the value.
  invisible(act$val)
}
