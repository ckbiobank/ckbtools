# Safe versions of regression model function

These functions are wrappers around the regression model functions that
will capture side-effects such as errors and warnings and return a list.
This is particularly useful when using the model functions iteratively.

## Usage

``` r
linear_model_safely(...)

logistic_model_safely(...)

cox_model_safely(...)
```

## Arguments

- ...:

  Arguments are passed to the model function.

## Value

A named list:

- result: What is returned by the model function.

- errors: Errors produced by the model function as a character string.

- warnings: Warnings produced by the model function as a character
  string.

- messages: Messages produced by the model function as a character
  string.

- output: Output produced by the model function as a character string.

## See also

[`linear_model()`](https://ckbiobank.github.io/ckbtools/reference/linear_model.md),
[`logistic_model()`](https://ckbiobank.github.io/ckbtools/reference/logistic_model.md),
[`cox_model()`](https://ckbiobank.github.io/ckbtools/reference/cox_model.md)
