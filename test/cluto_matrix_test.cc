#include <string>

#include "gtest/gtest.h"
#include "sskkm/format/cluto.h"
#include "test/test_util.h"

using namespace sskkm;
using namespace sskkm::cluto;

TEST(ClutoMatrixTest, CanDetermineMatrixType) {
  std::string invalid_header =
      "3\n"
      "1 1 2\n"
      "3 5 8\n"
      "13 21 34\n";
  EXPECT_EQ(MatrixType::InvalidMatrix, DetermineMatrixType(invalid_header));

  std::string dense_matrix =
      "3 3\n"
      "1 1 2\n"
      "3 5 8\n"
      "13 21 34\n";
  EXPECT_EQ(MatrixType::DenseMatrix, DetermineMatrixType(dense_matrix));

  std::string sparse_matrix =
      "3 10 9\n"
      "1 1.0 3 1.0 5 1.0\n"
      "4 3.0 7 2.0 10 1.0\n"
      "2 2.0 4 4.0 6 6.0\n";
  EXPECT_EQ(MatrixType::SparseMatrix, DetermineMatrixType(sparse_matrix));
}

TEST(ClutoMatrixTest, CanParseDenseMatrix) {
  std::string matrix =
      "3 4\n"
      "1 1 2 3\n"
      "5 8 13 21\n"
      "34 55 89 144\n";
  DenseMatrix matrix_expected(3, 4);
  matrix_expected <<
      1, 1, 2, 3,
      5, 8, 13, 21,
      34, 55, 89, 144;

  DenseMatrix matrix_parsed = ParseDenseMatrix(matrix);
  EXPECT_EQ(matrix_expected, matrix_parsed);

  DenseMatrix matrix_parsed_transposed = ParseDenseMatrix(matrix, true);
  EXPECT_EQ(matrix_expected.transpose(), matrix_parsed_transposed);
}

TEST(ClutoMatrixTest, CanParseSparseMatrix) {
  std::string matrix =
      "3 10 9\n"
      "1 1.0 3 1.0 5 1.0\n"
      "4 3.0 7 2.0 10 1.0\n"
      "2 2.0 4 4.0 6 6.0\n";
  SparseMatrix matrix_expected(3, 10);
  matrix_expected.insert(0, 0) = 1.0;
  matrix_expected.insert(2, 1) = 2.0;
  matrix_expected.insert(0, 2) = 1.0;
  matrix_expected.insert(1, 3) = 3.0;
  matrix_expected.insert(2, 3) = 4.0;
  matrix_expected.insert(0, 4) = 1.0;
  matrix_expected.insert(2, 5) = 6.0;
  matrix_expected.insert(1, 6) = 2.0;
  matrix_expected.insert(1, 9) = 1.0;
  matrix_expected.finalize();

  SparseMatrix matrix_parsed = ParseSparseMatrix(matrix);
  /*
    XXX: Why |EXPECT_EQ(matrix_expected, matrix_parsed)| generates an error (no
    match for 'operator==') ?
   */
  EXPECT_TRUE(matrix_expected == matrix_parsed);

  SparseMatrix matrix_parsed_transposed = ParseSparseMatrix(matrix, true);
  EXPECT_TRUE(matrix_expected.transpose() == matrix_parsed_transposed);
}

TEST(ClutoMatrixTest, RaiseExceptionOnTailingNewLines) {
  std::string tailing_newline =
      "3 4\n"
      "1 1 2 3\n"
      "5 8 13 21\n"
      "34 55 89 144\n"
      "\n";
  EXPECT_THROW(ParseDenseMatrix(tailing_newline), InvalidFormat);
  EXPECT_THROW(ParseSparseMatrix(tailing_newline), InvalidFormat);
}

TEST(ClutoMatrixTest, RaiseExceptionOnEmptyInput) {
  std::string empty = "";
  EXPECT_THROW(ParseDenseMatrix(empty), InvalidFormat);
  EXPECT_THROW(ParseSparseMatrix(empty), InvalidFormat);
}

TEST(ClutoMatrixTest, RaiseExceptionOnInvalidHeader) {
  std::string invalid_header =
      "3\n"
      "1 1 2\n"
      "3 5 8\n"
      "13 21 34\n";
  EXPECT_THROW(ParseDenseMatrix(invalid_header), InvalidFormat);
  EXPECT_THROW(ParseSparseMatrix(invalid_header), InvalidFormat);
}

TEST(ClutoMatrixTest, RaiseExceptionOnInvalidBody) {
  std::string invalid_body1 =
      "3 4\n"
      "1 1 2 3\n"
      "5 8 13 \n"
      "34 55 89 144\n";
  EXPECT_THROW(ParseDenseMatrix(invalid_body1), InvalidFormat);

  std::string invalid_body2 =
      "3 4\n"
      "1 1 2 3\n"
      "8 13 21\n"
      "34 55 89 144\n";
  EXPECT_THROW(ParseDenseMatrix(invalid_body2), InvalidFormat);

  std::string invalid_body3 =
      "3 10 9\n"
      "1 1.0 3 1.0 5 1.0\n"
      "4 3.0 7 2.0 10\n"
      "2 2.0 4 4.0 6 6.0\n";
  EXPECT_THROW(ParseSparseMatrix(invalid_body3), InvalidFormat);
}

TEST(ClutoMatrixTest, RaiseExceptionOnWrongSizeInput) {
  std::string wrong_size1 =
      "3 3\n"
      "1 1 2 3\n"
      "5 8 13 21\n"
      "34 55 89 144\n";
  EXPECT_THROW(ParseDenseMatrix(wrong_size1), InvalidFormat);

  std::string wrong_size2 =
      "3 5\n"
      "1 1 2 3\n"
      "5 8 13 21\n"
      "34 55 89 144\n";
  EXPECT_THROW(ParseDenseMatrix(wrong_size2), InvalidFormat);

  std::string wrong_size3 =
      "2 4\n"
      "1 1 2 3\n"
      "5 8 13 21\n"
      "34 55 89 144\n";
  EXPECT_THROW(ParseDenseMatrix(wrong_size3), InvalidFormat);

  std::string wrong_size4 =
      "4 4\n"
      "1 1 2 3\n"
      "5 8 13 21\n"
      "34 55 89 144\n";
  EXPECT_THROW(ParseDenseMatrix(wrong_size4), InvalidFormat);

  std::string wrong_size5 =
      "3 9 9\n"
      "1 1.0 3 1.0 5 1.0\n"
      "4 3.0 7 2.0 10 1.0\n"
      "2 2.0 4 4.0 6 6.0\n";
  EXPECT_THROW(ParseSparseMatrix(wrong_size5), InvalidFormat);

  std::string wrong_size6 =
      "2 10 9\n"
      "1 1.0 3 1.0 5 1.0\n"
      "4 3.0 7 2.0 10 1.0\n"
      "2 2.0 4 4.0 6 6.0\n";
  EXPECT_THROW(ParseSparseMatrix(wrong_size6), InvalidFormat);

  std::string wrong_size7 =
      "4 10 9\n"
      "1 1.0 3 1.0 5 1.0\n"
      "4 3.0 7 2.0 10 1.0\n"
      "2 2.0 4 4.0 6 6.0\n";
  EXPECT_THROW(ParseSparseMatrix(wrong_size7), InvalidFormat);
}

TEST(ClutoMatrixTest, RaiseExceptionOnWrongNumberOfNonZeros) {
  std::string wrong_number_of_nonzeros =
      "3 10 7\n"
      "1 1.0 3 1.0 5 1.0\n"
      "4 3.0 7 2.0 10 1.0\n"
      "2 2.0 4 4.0 6 6.0\n";
  EXPECT_THROW(
      ParseSparseMatrix(wrong_number_of_nonzeros), InvalidFormat);
  EXPECT_THROW(
      ParseSparseMatrix(wrong_number_of_nonzeros, true), InvalidFormat);
}
