// This class is used to compress OpenCV Decision Trees. Due to their versaility

#ifndef _CvDTree_
#define _CvDTree_

#include <cstddef>
#include <memory>

#include <opencv/cv.h>
#include <opencv/ml.h>



struct GiDTreeSplit {
  // The index of the split variable.
  int var_idx;

  // Whether the categorical split is inversed.
  int inversed;

  // The split point.
  union {
    int subset[2];
    struct {
      float c;
      int split_point;
    } ord;
  };

  GiDTreeSplit(CvDTreeSplit* copy)
      : var_idx(copy->var_idx),
        inversed(copy->inversed) {
    subset[0] = copy->subset[0];
    subset[1] = copy->subset[1];
  }
};

struct GiDTreeNode {
  using ChildPtr = std::unique_ptr<GiDTreeNode>;
  using SplitPtr = std::unique_ptr<GiDTreeSplit>;

  // The predicted value at this node.
  double value;

  // The split for this node.
  SplitPtr split;

  // The children of this node.
  ChildPtr left, right;

  GiDTreeNode(const CvDTreeNode* copy)
      : value(copy->value) {
    // If non-leaf, children and split are duplicated.
    // Otherwise, they left as the default values, null pointers.
    if (copy->left != nullptr) {
      split.reset(new GiDTreeSplit(copy->split));
      left.reset(new GiDTreeNode(copy->left));
      right.reset(new GiDTreeNode(copy->right));
    }
  }
};

class GiDTree {
 public:
  using RootPtr = std::unique_ptr<GiDTreeNode>;

  GiDTree(CvDTree* copy);

  ~GiDTree();

  double predict(const arma::vec& sample) const;

 private:
  // Root of this binary tree.
  RootPtr root;

  // Information regarding type specifications.
  arma::ivec var_idx, cat_map, cat_ofs, var_type, cat_count;

  int is_buf_16u;
};

GiDTree::GiDTree(CvDTree* copy)
    : root(new GiDTreeNode(copy->get_root())),
      is_buf_16u(copy->get_data()->is_buf_16u) {
  // CvMat converted to Mat to deal with possible null pointers.
  cv::Mat var_idx_mat(  copy->get_data()->var_idx);
  cv::Mat cat_map_mat(  copy->get_data()->cat_map);
  cv::Mat cat_ofs_mat(  copy->get_data()->cat_ofs);
  cv::Mat var_type_mat( copy->get_data()->var_type);
  cv::Mat cat_count_mat(copy->get_data()->cat_count);

  // Casting to const pointers forces the inner data to be copied.
  var_idx   = arma::ivec((const int*) var_idx_mat.data,   var_idx_mat.cols);
  cat_map   = arma::ivec((const int*) cat_map_mat.data,   cat_map_mat.cols);
  cat_ofs   = arma::ivec((const int*) cat_ofs_mat.data,   cat_ofs_mat.cols);
  var_type  = arma::ivec((const int*) var_type_mat.data,  var_type_mat.cols);
  cat_count = arma::ivec((const int*) cat_count_mat.data, cat_count_mat.cols);
}

GiDTree::~GiDTree() {
  var_idx.reset();
  cat_map.reset();
  cat_ofs.reset();
  var_type.reset();
  cat_count.reset();
}

double GiDTree::predict(const arma::vec& sample) const {
  GiDTreeNode* node(root.get());
  arma::vec catbuf(cat_count.n_elem);
  catbuf.fill(-1);

  while (node->left != nullptr) {
    GiDTreeSplit* split(node->split.get());
    int dir = 0;
    int vi = split->var_idx; // Index of splitting variable.
    int ci = var_type(vi); // Type of splitting variable.
    double val = sample(vi); // Value of splitting variable

    if (ci < 0) { // Numerical split, simple bounds check.
      dir = (val <= split->ord.c) ? -1 : 1;
    } else {
      int c = -1;
      if (c < 0) {
        int a = c = cat_ofs(ci);
        int b = (ci + 1 >= cat_ofs.n_elem) ? cat_map.n_elem : cat_ofs(ci + 1);

        int ival = (int) std::round(val);
        int sh = 0;

        while (a < b) {
          sh++;
          c = (a + b) >> 1;
          if (ival < cat_map(c))
            b = c;
          else if (ival > cat_map(c))
            a = c + 1;
          else
            break;
        }

        if (c < 0 || ival != cat_map(c))
          continue;

        catbuf(ci) = c -= cat_ofs(ci);
      }
      c = ((c == 65535) && is_buf_16u) ? -1 : c;
      dir = (2 * ((split->subset[c >> 5] & (1 << (c & 31))) == 0) - 1);
    }

    if (split->inversed)
      dir = -dir;

    node = (dir < 0) ? node->left.get() : node->right.get();
  }

  return node->value;
}

#endif
