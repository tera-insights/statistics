#ifndef _CLUS_EM_H_
#define _CLUS_EM_H_

#include "distribution.h"

#include <cstddef>
#include <cmath>
#include <array>
#include <armadillo>

namespace CLUS {

template<class T_Distribution>
class ExpectationMaximizer {
 public:
  using em_type = ExpectationMaximizer<T_Distribution>;
	using size_type = std::size_t;
	using dist_type = T_Distribution;
	using dist_array = std::array<dist_type, 2>;

	using value_type = double;
	using vector_type = arma::vec;

	struct ConfigType {
		size_type n_disc_rounds;
		size_type n_disc_iter;
		size_type n_learn_iter;
		size_type max_trials;
	};

 private:
	enum class State {
		build_overall,
		discover_children_init,
		discover_children_iter,
		build_children,
		done
	};

	// Current state of EM state machine
	State state;

	// The number of input dimensions
	const size_type inDim;

	// Overall distribution of the data
	// Built in state `build_overall` or given in constructor
	dist_type overallDist;

	// Prospective child distributions
	dist_array childDists;

	// Current best child distributions
	dist_array bestDists;

	// Used to determine a good set of child distributions.
	double likelihood;
	double best_likelihood;

	// Current iteration/restart number
	size_type iteration;

	// Number of random restarts to perform, to find a good set of children
	const size_type num_discovery_rounds;
	size_type discovery_round;

	// Number of iterations to perform per random restart
	const size_type num_discovery_iter;

	// Number of iterations to perform to teach the distributions
	const size_type num_iterations;

	// Maximum number of trials to perform. If during a trial, we find that
	// one of the child distributions has become invalid, we will attempt to
	// restart as long as num_trials < max_trials. If we fail to find good
	// distributions after max_trials attempts, we will abandon the search,
	// essentially declaring that the current node should become a leaf.
	const size_type max_trials;
	size_type num_trials;

	// Keeps track of wheter or not the child distributions are valid.
	// They may be invalid if we failed to find good distributions on all
	// trials.
	bool distsValid;

 public:
	// With no overall distribution given
	ExpectationMaximizer(size_type _inDim, const ConfigType& config)
      : state(State::build_overall),
        inDim(_inDim),
        overallDist(_inDim),
        childDists({{_inDim, _inDim}}),
        bestDists({{_inDim, _inDim}}),
        likelihood(0),
        best_likelihood(-arma::datum::inf),
        iteration(0),
        num_discovery_rounds(config.n_disc_rounds),
        discovery_round(0),
        num_discovery_iter(config.n_disc_iter),
        num_iterations(config.n_learn_iter),
        max_trials(config.max_trials),
        num_trials(0),
        distsValid(false) {
  }

	// Copy constructor
	ExpectationMaximizer(const em_type& other)
      : state(other.state),
        inDim(other.inDim),
        overallDist(other.overallDist),
        childDists(other.childDists),
        bestDists(other.bestDists),
        likelihood(other.likelihood),
        best_likelihood(other.best_likelihood),
        iteration(other.iteration),
        num_discovery_rounds(other.num_discovery_rounds),
        discovery_round(other.discovery_round),
        num_discovery_iter(other.num_discovery_iter),
        num_iterations(other.num_iterations),
        max_trials(other.max_trials),
        num_trials(other.num_trials),
        distsValid(other.distsValid) {
	}

	void start_round() {
    std::cout << "Starting EM state " << state_name() << std::endl;
		switch(state) {
			case State::build_overall:
				break;
			case State::discover_children_init:
				/* fallthrough */
			case State::discover_children_iter:
				likelihood = 0.0;
				iteration++;
				break;
			case State::build_children:
				likelihood = 0.0;
				iteration++;
				break;
			case State::done: /* fallthrough */
			default:
				// do nothing
				break;
		}
	}

	void learn(const vector_type& data, double prob = 1) {
		switch (state) {
			case State::build_overall:
				overallDist.Update(data, prob);
				break;
			case State::discover_children_init:
			case State::discover_children_iter:
				update_children(data, prob, childDists);
				break;
			case State::build_children:
				update_children(data, prob, bestDists);
				break;
			case State::done: /* fallthrough */
			default:
				break;
		}
	}

  void Merge(const em_type& other) {
    switch (state) {
      case State::build_overall:
        overallDist.Merge(other.overallDist);
        break;
      case State::discover_children_init:
			case State::discover_children_iter:
        childDists[0].Merge(other.childDists[0]);
        childDists[1].Merge(other.childDists[1]);
        break;
      case State::build_children:
        bestDists[0].Merge(other.bestDists[0]);
        bestDists[1].Merge(other.bestDists[1]);
        break;
      case State::done:
      default:
        break;
    }
  }

	void end_round(double convergenceLimit ) {
		double convFactor, c0, c1;

		State next_state = state;

		switch (state) {
			case State::build_overall:
				overallDist.Estimate();
				next_state = State::discover_children_init;
				iteration = 0;
        ResetChildren();
				break;
			case State::discover_children_init:
				next_state = State::discover_children_iter;
				/* fallthrough*/
			case State::discover_children_iter:
				if (iteration >= num_discovery_iter) {
					// Completed iterating on this set of children for now
          std::cout << "Comparing likelihoods, best: " << best_likelihood << " curr: " << likelihood << std::endl;
					if (likelihood > best_likelihood) {
            std::cout << "Improved likelihood" << std::endl;
						best_likelihood = likelihood;
						bestDists = childDists;
					}

					discovery_round++;
					iteration = 0;
					next_state = State::discover_children_init;
          ResetChildren();
				}

				// If we have performed enough restarts, continue on to the
				// child building phase
				if (discovery_round >= num_discovery_rounds) {
					iteration = 0;
					next_state = State::build_children;

					// Ensure the best children have updated parameters
					bestDists[0].Estimate();
					bestDists[1].Estimate();
				}
				break;
			case State::build_children:
				c0 = bestDists[0].Estimate();
				c1 = bestDists[1].Estimate();
				convFactor = (c0 + c1) / 2.0;

				if (bestDists[0].HasZeroWeight() || bestDists[1].HasZeroWeight()) {
					// This trial has failed
					num_trials++;
					if (num_trials < max_trials) {
						// Start a new trial
						iteration = 0;
						discovery_round = 0;
						best_likelihood = -arma::datum::inf;

						next_state = State::discover_children_init;
            ResetChildren();
					} else {
						// Failed completely, finish with invalid children.
            std::cout << "Children are invalid. No splitting." << std::endl;
						distsValid = false;

						next_state = State::done;
					}
				} else {
					if (convFactor <= convergenceLimit || iteration >= num_iterations) {
						// Trial has successfully completed.
						bestDists[0].DenormalizeParameters(overallDist);
						bestDists[1].DenormalizeParameters(overallDist);
						distsValid = true;

						next_state = State::done;
					}

					// Otherwise continue with the current set of iterations.
				}
			case State::done: /* fallthrough */
			default:
				break;
		}

		state = next_state;
	}

	// Whether or not the EM is done learning, and can be used to predict.
	bool inline done() const {
		return state == State::done;
	}

	bool inline should_split() const {
		return distsValid;
	}

	double inline prob_left(const vector_type& data, double p) {
		return bestDists[0].PDF(data);
	}

	double inline prob_right(const vector_type& data, double p) {
		return bestDists[1].PDF(data);
	}

 private:
	void update_children(const vector_type& data, double p, dist_array& children) {
		double p0, p1, p01;

		// Normalize the data wrt the overall distribution
		vector_type norm_data = overallDist.Normalize(data);

		// Probability that each distribution contains the given point.
		p0 = children[0].PDF(norm_data);
		p1 = children[1].PDF(norm_data);
		p01 = p0 + p1;
    if (p01 == 0) {
      p0 = p1 = 0.5;
      p01 = 1;
    }

		if (std::isfinite(p01) && std::isfinite(likelihood)) {
			likelihood += std::log(p01);
      if (p01 / p <= 0) {
        std::cout << std::endl;
        std::cout << "Log domain error occurred. ";
        std::cout << " p01: " << p01 << " p " << p;
        std::cout << " p1: " << p1 << " p0 " << p0 << std::endl;
        std::cout << "Distribution0: weightIsZero: " << children[0].HasZeroWeight() << std::endl;
        std::cout << "Distribution1: weightIsZero: " << children[1].HasZeroWeight() << std::endl;
        std::cout << "Mu: " << overallDist.mu.t() << "Cholesky: " << std::endl << overallDist.chol
                  << "Cholesky inv" << std::endl << overallDist.chol.i();
        std::cout << "Input vector" << data.t() << "normalized" << norm_data.t() << std::endl;
      }
			children[0].Update(norm_data, p * p0 / p01);
			children[1].Update(norm_data, p * p1 / p01);
		}
	}

  void ResetChildren() {
		childDists[0].RandomDistribution(2);
		childDists[1].RandomDistribution(2);
  }

  const char* state_name() const {
    switch(state) {
      case State::build_overall:
        return "build_overall";
      case State::discover_children_init:
        return "discover_children_init";
      case State::discover_children_iter:
        return "discover_chilren_iter";
      case State::build_children:
        return "build_children";
      case State::done:
        return "done";
      default:
        return "invalid state";
    }
  }
};

}

#endif // _CLUS_EM_H_
