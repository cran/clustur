//
// Created by Gregory Johnson on 9/16/24.
//

#include <testthat.h>
#include "Tests/ColumnReaderTestFixture.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test contexta. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

context("ColumnDistanceMatrixReader Test") {
    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Column Reader reads column files") {
        Rcpp::Environment pkg = Rcpp::Environment::namespace_env("testthat");
        Rcpp::Environment clustur = Rcpp::Environment::namespace_env("clustur");
        const Rcpp::Function test_path = pkg["test_path"];
        const Rcpp::Function read_count = clustur["read_count"];
        const std::string path = Rcpp::as<std::string>(test_path("extdata", "amazon_column.dist"));
        const std::string countTablePath = Rcpp::as<std::string>(test_path("extdata", "amazon.count_table"));
        const Rcpp::DataFrame df = read_count(countTablePath);
        ColumnDistanceMatrixReader reader(0.2, false);
        reader.CreateCountTableAdapter(df);
        reader.Read(path);
        ColumnReaderTestFixture fixture;
        bool result = fixture.TestReadColumnFile(path, df, reader.GetListVector()->getNumSeqs());
        expect_true(result);
        result = fixture.TestReadColumnFile("", df, 1);
        expect_false(result);
    }

    test_that("Column Reader properly creates a sparse distance matrix"){
        ColumnReaderTestFixture fixture;
         Rcpp::Environment pkg = Rcpp::Environment::namespace_env("testthat");
        Rcpp::Environment clustur = Rcpp::Environment::namespace_env("clustur");
        const Rcpp::Function test_path = pkg["test_path"];
        const Rcpp::Function read_count = clustur["read_count"];
        const std::string path = Rcpp::as<std::string>(test_path("extdata", "amazon_column.dist"));
        const std::string countTablePath = Rcpp::as<std::string>(test_path("extdata", "amazon.count_table"));
        const Rcpp::DataFrame df = read_count(countTablePath);
        ColumnDistanceMatrixReader reader(0.2, false);
        reader.CreateCountTableAdapter(df);
        reader.Read(path);
        bool result = fixture.TestGetDistanceMatrix(path, df, true);
        expect_true(result);
        result = fixture.TestGetDistanceMatrix("", df, false);
        expect_true(result);

    }
      test_that("Column Reader properly creates a list vector"){
        Rcpp::Environment pkg = Rcpp::Environment::namespace_env("testthat");
        Rcpp::Environment clustur = Rcpp::Environment::namespace_env("clustur");
        const Rcpp::Function test_path = pkg["test_path"];
        const Rcpp::Function read_count = clustur["read_count"];
        const std::string path = Rcpp::as<std::string>(test_path("extdata", "amazon_column.dist"));
        const std::string countTablePath = Rcpp::as<std::string>(test_path("extdata", "amazon.count_table"));
        const Rcpp::DataFrame df = read_count(countTablePath);
        ColumnDistanceMatrixReader reader(0.2, false);
        reader.CreateCountTableAdapter(df);
        reader.Read(path);
        ColumnReaderTestFixture fixture;
        bool result = fixture.TestGetListVector(path, df, reader.GetListVector()->getNumSeqs());
        expect_true(result);
        result = fixture.TestGetListVector("", df, 0);
        expect_true(result);
    }
}
