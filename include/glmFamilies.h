#ifndef _GLM_FAMILIES_H_
#define _GLM_FAMILIES_H_

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>

#define PI (boost::math::double_constants::pi)

//LINK FUNCTIONS

/* Identity link */
inline double _identity(double x){ return x; }

/* Inverse link - negative inverse */
inline double _inverse(double x){ return 1.0/x; }

/* Inverse squared link */
inline double _inverseSquared(double x){ return 1.0/(x * x); }

/* Logarithmic link - natural log */
inline double _log(double x){ return std::log(x); }

/* Logit link - logarithmic odds */
inline double _logit(double x){ return std::log(x/(1 - x)); }

/* Cloglog link - complementary log log */
inline double _cloglog(double x){ return std::log(-std::log(1 - x)); }

/* Sqrt link */
inline double _sqrt(double x){ return std::sqrt(x); }

/* Probit link - normal distribution quantile function */
inline double _probit(double x){ return std::sqrt(2) * boost::math::erf_inv(2*x - 1); }

/* Cauchit link - quantile function, good for outliers that cause the predictor to be quite large */
inline double _cauchit(double x){ return std::tan(PI *(x - 0.5));}
//INVERSE LINK FUNCTIONS

/* Identity link */
inline double inverse_identity(double x){ return x; }

/* Inverse link - negative inverse */
inline double inverse_inverse(double x){ return 1.0/x; }

/* Inverse squared link */
inline double inverse_inverseSquared(double x){ return 1.0/(std::sqrt(x)); }

/* Logarithmic link - natural log */
inline double inverse_log(double x){ return std::exp(x); }

/* Logit link - logarithmic odds */
inline double inverse_logit(double x){ return 1.0/(1.0 + std::exp(-x)); }

/* Cloglog link - complementary log log */
inline double inverse_cloglog(double x){ return 1.0 - std::exp(-std::exp(x)); }

/* Sqrt link */
inline double inverse_sqrt(double x){ return x * x; }

/* Probit link - normal distribution quantile function */
inline double inverse_probit(double x){ return (boost::math::erf(x / std::sqrt(2)) + 1.0)/2.0; }

/* Cauchit link */
inline double inverse_cauchit(double x){ return std::atan(x) / PI + 0.5; }

//DERIVATIVE LINK FUNCTIONS

/* Identity link */
inline double deriv_identity(double x){ return 1; }

/* Inverse link - negative inverse */
inline double deriv_inverse(double x){ return -1.0/(x * x); }

/* Inverse squared link */
inline double deriv_inverseSquared(double x){ return -2.0/(x * x * x); }

/* Logarithmic link - natural log */
inline double deriv_log(double x){ return 1.0/x; }

/* Logit link - logarithmic odds */
inline double deriv_logit(double x){ return 1.0 /(x * (1.0 - x)); }


/* Cloglog link - complementary log log */
inline double deriv_cloglog(double x){ return 1.0/((1 - x) * log (1 - x)); }

/* Sqrt link */
inline double deriv_sqrt(double x){ return 1.0/(2 *  std::sqrt(x)); }

/* Probit link - normal distribution quantile function */
inline double deriv_probit(double x){ return  sqrt(2 * PI) * exp(boost::math::erf_inv(pow(2 * x - 1, 2))); }

/* Cauchi link */
inline double deriv_cauchit(double x){ return PI * std::pow(1/std::cos(PI * (x - 0.5)), 2); }
//VARIANCE FUNCTIONS

/* Guassian Variance - input: stdDev x */
inline double var_gaussian(double x){ return 1; }

/* Binomial Variance - input: prob x */
inline double var_binomial(double x){ return x * (1 - x); }

/* Negative Binomial Variance - input: prob x, numFail n */
inline double var_negativeBinomial(double x, int n){ return x + x * x / n;}

/* Poisson Variance - input: parameter x */
inline double var_poisson(double x){ return x; }

/* Exponential Variance - input: parameter x */
inline double var_exponential(double x){ return x * x; }

/* Gamma Variance - input: scale x */
inline double var_gamma(double x){ return x * x; }

/* Inverse Gaussian Function - input: parameter x, parameter y */
inline double var_inverseGaussian(double x){ return x * x * x; }

/* Negative Bionimial - input: probability x, failues y */
inline double var_negativeBinomial(double x){ return x / ((1 - x) * (1 - x)); }

//DISPERSION FUNCTIONS

/* Binomial Dispersion - input: numTrials n */
inline double dispersion_binomial(){ return 1.0; }

/* Poisson Dispersion - constant 1 */
inline double dispersion_poisson() { return 1; }

/* Gamma Dispersion - shape parameter */
inline double dispersion_gamma(double x){ return 1.0 / x; }

/* Inverse Gaussian Dispersion - shape parameter (sigma) */
inline double dispersion_inverseGaussian(double x) { return x * x; }

/* Gaussian Dispersion - shape paramter (sigma) */
inline double dispersion_gaussian(double x) { return x * x; }

/* expoential dispersion - constant 1 */
inline double dispersion_exponential(){ return 1;}

//DEVIANCE FUNCTIONS - all take y as response and u as predicted mean.

/*Gaussian Deviance - residual sum of squares*/
inline double deviance_gaussian(double y, double mu){ return std::pow(y - mu, 2);}

/*Poisson Deviance */
inline double deviance_poisson(double y, double mu){return 2*(y * std::log(y / mu) - (y - mu));}

/*Binomial Deviance - n is count*/
inline double deviance_binomial(double y, double mu){return 2 * (y * std::log(y / mu) + (1 - y) * std::log((1 - y) / (1 - mu)));}

/*Gamma Deviance */
inline double deviance_gamma(double y, double mu){return 2 * (-std::log(y / mu) + (y - mu) / mu); }

/*Inverse Gaussian Deviance */
inline double deviance_inverseGaussian(double y, double mu){return std::pow(y - mu, 2) / (mu * mu * y); }

#endif // _GLM_FAMILIES_H_
