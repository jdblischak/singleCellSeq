# Run all test files in tests/

test_files <- list.files(pattern = "^test")
stopifnot(file.exists(test_files))

for (test in test_files) {
  message("# Running ", test)
  source(test)
}
