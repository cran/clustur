//
// Created by Gregory Johnson on 7/12/24.
//

#include <testthat.h>
#include "Tests/MatrixAdapterTestFixture.h"
#include "Adapters/MatrixAdapter.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("MatrixAdapter Test") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.

     test_that("Matrix Adapter can get its list vector") {
         MatrixAdapterTestFixture fixture;
         bool result = fixture.TestGetListVector(false, false);
         expect_true(result);
         result = fixture.TestGetListVector(true, true);
         expect_true(result);
         result = fixture.TestGetListVector(false, true);
         expect_false(result);
     }
     test_that("Matrix Adapter can create sparse matrices") {
         MatrixAdapterTestFixture fixture;
         bool result = fixture.TestCreateSparseMatrix(6);
         expect_true(result);
         result = fixture.TestCreateSparseMatrix(0);
         expect_false(result);
     }
 }
