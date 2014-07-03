#ifndef SSKKM_FORMAT_BASE_H_
#define SSKKM_FORMAT_BASE_H_

#include <boost/tokenizer.hpp>

#include "sskkm/base.h"

namespace sskkm {

class InvalidFormat : public std::invalid_argument {
 public:
  explicit InvalidFormat(const std::string& what_arg) :
      std::invalid_argument(what_arg) {}
};

namespace internal {

boost::char_delimiters_separator<char> separator_(true, "\n", " ");

}  // namespace internal

}  // namespace sskkm

#endif  // SSKKM_FORMAT_BASE_H_
