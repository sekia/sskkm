#ifndef SSKKM_FORMAT_CLUTO_H_
#define SSKKM_FORMAT_CLUTO_H_

#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>

#include "sskkm/base.h"
#include "sskkm/format/base.h"

namespace sskkm {

namespace cluto {

enum MatrixType { kInvalidMatrix, kDenseMatrix, kSparseMatrix };

namespace internal {

/**
   Parses a header and sets results to the given pointers. |iter| will be bumped
   forward to right next at the end of the header.
*/
inline void ReadMatrixHeader(
    boost::tokenizer<>::iterator &iter,
    const boost::tokenizer<>::iterator &iter_end,
    unsigned *num_rows,
    unsigned *num_cols,
    unsigned *num_nonzeros = 0) throw (InvalidFormat) {
  if (iter == iter_end) { throw InvalidFormat("No header"); }
  try {
    *num_rows = boost::lexical_cast<unsigned>(*iter);
    ++iter;
    if (iter == iter_end) { throw InvalidFormat("No number of rows"); }
    *num_cols = boost::lexical_cast<unsigned>(*iter);
    ++iter;
    if (iter == iter_end) { throw InvalidFormat("No number of columns"); }
    if (num_nonzeros != 0) {
      *num_nonzeros = boost::lexical_cast<unsigned>(*iter);
      ++iter;
      if (iter == iter_end) {
        throw InvalidFormat("No number of non-zero coefficients");
      }
    }
    if (*iter != "\n") {
      throw InvalidFormat("Header includes not expected properties");
    }
    ++iter;
  } catch (const boost::bad_lexical_cast &e) {
    throw InvalidFormat("Not looks like a integer");
  }
}

inline DenseMatrix ReadDenseMatrix(
    boost::tokenizer<>::iterator &iter,
    const boost::tokenizer<>::iterator &iter_end,
    unsigned num_rows,
    unsigned num_cols,
    bool transposed) throw (InvalidFormat) {
  unsigned p = num_rows, q = num_cols;
  if (transposed) { std::swap(p, q); }
  DenseMatrix matrix(num_rows, num_cols);
  for (unsigned i = 0; i < p; ++i, ++iter) {
    for (unsigned j = 0; j < q; ++j, ++iter) {
      if (iter == iter_end) {
        throw InvalidFormat("Not enough coefficients");
      }
      try {
        if (transposed) {
          matrix(j, i) = boost::lexical_cast<double>(*iter);
        } else {
          matrix(i, j) = boost::lexical_cast<double>(*iter);
        }
      } catch (const boost::bad_lexical_cast &e) {
        throw InvalidFormat(
            "Coefficient does not look a like real number");
      }
    }
    if (iter != iter_end && *iter != "\n") {
      throw InvalidFormat(
          transposed ? "Row index out of range" : "Column index out of range");
    }
  }

  if (iter != iter_end) {
    throw InvalidFormat(
        transposed ? "Column index out of range" : "Row index out of range");
  }

  return matrix;
}

inline SparseMatrix ReadSparseMatrix(
    boost::tokenizer<>::iterator &iter,
    const boost::tokenizer<>::iterator &iter_end,
    unsigned num_rows,
    unsigned num_cols,
    unsigned num_nonzeros,
    bool transposed) throw (InvalidFormat) {
  unsigned p = num_rows, q = num_cols;
  if (transposed) { std::swap(p, q); }
  std::vector<SparseMatrixCoefficient> matrix_values;
  matrix_values.reserve(num_nonzeros);
  unsigned count_nonzeros = 0;
  unsigned i;
  bool no_tailing_newline = false;
  for (i = 0; iter != iter_end; ++i, ++iter) {
    if (i >= p) {
      throw InvalidFormat(
          transposed ? "Column index out of range" : "Row index out of range");
    }
    while (iter != iter_end && *iter != "\n") {
      unsigned j;
      double value;
      try {
        j = boost::lexical_cast<unsigned>(*iter);
        ++iter;
        value = boost::lexical_cast<double>(*iter);
        ++iter;
      } catch (const boost::bad_lexical_cast &e) {
        throw InvalidFormat(
            transposed ?
            "Row index or coefficient value do not look like a number" :
            "Column index or coefficient value do not look like a number");
      }

      // In CLUTO format, vector indices start with 1.
      if (j < 1 || j > q) {
        throw InvalidFormat(
            transposed ?
            "Row index out of range" :
            "Column index out of range");
      }

      --j;  // In Eigen, vector indices start with 0, off course.
      matrix_values.push_back(
          transposed ?
          SparseMatrixCoefficient(j, i, value) :
          SparseMatrixCoefficient(i, j, value));
      ++count_nonzeros;
    }
    if (iter == iter_end) {
      no_tailing_newline = true;
      break;
    }
  }
  if (i != (transposed ? num_cols : num_rows) - (no_tailing_newline ? 1 : 0)) {
    throw InvalidFormat(
        transposed ?
        "Number of columns does not match to declared" :
        "Number of rows does not match to declared");
  }
  if (count_nonzeros != num_nonzeros) {
    throw InvalidFormat(
        "Number of non-zero coefficients does not match to declared");
  }
  SparseMatrix matrix(num_rows, num_cols);
  matrix.setFromTriplets(matrix_values.begin(), matrix_values.end());
  return matrix;
}

}  // namespace internal

inline MatrixType DetermineMatrixType(const std::string &src) throw () {
  boost::tokenizer<> tokenizer(src, ::sskkm::internal::separator_);
  unsigned num_header_elements = 0;
  for (boost::tokenizer<>::iterator iter = tokenizer.begin();
       iter != tokenizer.end() && *iter != "\n";
       ++iter) {
    try {
      boost::lexical_cast<unsigned>(*iter);
    } catch (const boost::bad_lexical_cast &e) {
      return kInvalidMatrix;
    }
    ++num_header_elements;
  }

  switch (num_header_elements) {
    case 2:
      return kDenseMatrix;
    case 3:
      return kSparseMatrix;
    default:
      return kInvalidMatrix;
  }
}

inline DenseMatrix ParseDenseMatrix(
    const std::string &src,
    bool transposed = false) throw (InvalidFormat) {
  if (DetermineMatrixType(src) != kDenseMatrix) {
    throw InvalidFormat("Invalid dense matrix foramt");
  }

  boost::tokenizer<> tokenizer(src, ::sskkm::internal::separator_);
  boost::tokenizer<>::iterator iter = tokenizer.begin();
  boost::tokenizer<>::iterator iter_end = tokenizer.end();
  unsigned num_rows, num_cols;
  internal::ReadMatrixHeader(iter, iter_end, &num_rows, &num_cols);
  if (transposed) { std::swap(num_rows, num_cols); }
  return internal::ReadDenseMatrix(
      iter, iter_end, num_rows, num_cols, transposed);
}

inline SparseMatrix ParseSparseMatrix(
    const std::string &src,
    bool transposed = false) throw (InvalidFormat) {
  if (DetermineMatrixType(src) != kSparseMatrix) {
    throw InvalidFormat("Invalid sparse matrix foramt");
  }

  boost::tokenizer<> tokenizer(src, ::sskkm::internal::separator_);
  boost::tokenizer<>::iterator iter = tokenizer.begin();
  boost::tokenizer<>::iterator iter_end = tokenizer.end();
  unsigned num_rows, num_cols, num_nonzeros;
  internal::ReadMatrixHeader(
      iter, iter_end, &num_rows, &num_cols, &num_nonzeros);
  if (transposed) { std::swap(num_rows, num_cols); }
  return internal::ReadSparseMatrix(
      iter, iter_end, num_rows, num_cols, num_nonzeros, transposed);
}

}  // namespace cluto

}  // namespace sskkm

#endif  // SSKKM_FORMAT_CLUTO_H_
