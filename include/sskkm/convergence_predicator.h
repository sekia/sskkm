#ifndef SSKKM_CONVERGENCE_PREDICATOR_H_
#define SSKKM_CONVERGENCE_PREDICATOR_H_

namespace sskkm {

class ConvergencePredicator {
 public:
  ConvergencePredicator(
      unsigned num_max_iterations,
      unsigned num_least_iterations,
      double epsilon) throw () :
      epsilon_(epsilon),
      iteration_count_(0),
      num_least_iterations_(num_least_iterations),
      num_max_iterations_(num_max_iterations) {}

  // Use default copy constructor and assignment operator.

  bool operator()(const double prev_score, const double curr_score) throw () {
    ++iteration_count_;
    if (iteration_count_ <= num_least_iterations_) { return false; }
    if (iteration_count_ >= num_max_iterations_) { return true; }
    if (abs(prev_score - curr_score) / curr_score < epsilon_) { return true; }
    return false;
  }

 private:
  ConvergencePredicator();

  double epsilon_;
  unsigned iteration_count_;
  unsigned num_max_iterations_;
  unsigned num_least_iterations_;
};

}  // namespace sskkm

#endif  // SSKKM_CONVERGENCE_PREDICATOR_H_
