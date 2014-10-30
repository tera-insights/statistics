/* Alin: This code needs to be moved into a global location */

// These functions are used to serialize Armadillo object using JSON. Although
// columns and rows are essentially specialized matrices in Armadillo, they are
// serialized so that there is a clear distinction between matrices and vectors.
// However, rows and columns are treated identically. For example, for vector v,
// the serialization of v and v.t() would be indistinguishable.

#ifndef _ARMA_JSON_
#define _ARMA_JSON_

#include <armadillo>
#include <jsoncpp/json/json.h>

template<class T>
void ToJson(const arma::Col<T>& src, Json::Value& dest) {
  dest["__type__"] = "vector";
  dest["n_elem"] = src.n_elem;
  Json::Value content(Json::arrayValue);
  for (int i = 0; i < src.n_elem; i++)
    content[i] = src(i);
  dest["data"] = content;
}

template<class T>
void ToJson(const arma::Row<T>& src, Json::Value& dest) {
  dest["__type__"] = "vector";
  dest["n_elem"] = src.n_elem;
  Json::Value content(Json::arrayValue);
  for (int i = 0; i < src.n_elem; i++)
    content[i] = src(i);
  dest["data"] = content;
}

template<class T>
void ToJson(const arma::Mat<T>& src, Json::Value& dest) {
  dest["__type__"] = "matrix";
  dest["n_rows"] = src.n_rows;
  dest["n_cols"] = src.n_cols;
  Json::Value content(Json::arrayValue);
  for (int i = 0; i < src.n_rows; i++)
    for (int j = 0; j < src.n_cols; j++)
	    content[i * src.n_cols + j] = src(i, j);
  dest["data"] = content;
}

// Copied from the ostream operator << definition in Armadillo source. Probably
// converts src to one of the basic types (Mat<T>, etc) and uses above. This is
// necessary due to the templating engine of Armadillo not recognizing properly
// recognizing that the above templated functions can be used on some of the
// more complicated templated types made by Armadillo.
template<typename eT, typename T1>
void ToJson(const arma::Base<eT,T1>& src, Json::Value & dest ){
    const arma::unwrap<T1> tmp(src.get_ref());
    ToJson(tmp.M, dest);
}

#endif //  _ARMA_JSON_
