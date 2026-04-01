// src/pval_correct.cpp
#include <cmath>
#include <cpp11.hpp>
#include <limits>

[[cpp11::register]]
cpp11::writable::doubles pval_correct_cpp(cpp11::doubles pval_matrix) {
  cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(pval_matrix.attr("dim"));
  if (dim.size() != 2)
    cpp11::stop("Input must be a matrix");
  int nrow = dim[0];
  int ncol = dim[1];
  if (nrow != ncol)
    cpp11::stop("Input matrix must be square");
  int p = nrow;

  auto idx = [&](int row, int col) -> int {
    // 1-based row,col to 0-based index in column-major storage
    return (row - 1) + (col - 1) * p;
  };

  cpp11::writable::doubles corrected(p * p);
  // corrected[p, ] <- pval_matrix[p, p:1]
  for (int col = 1; col <= p; ++col) {
    int src_col = p - (col - 1);
    corrected[idx(p, col)] = pval_matrix[idx(p, src_col)];
  }

  for (int var = 1; var <= p; ++var) {
    int k0 = var;
    int orig0 = p - ((k0 - 1) % p);
    double pval_var = pval_matrix[idx(p, orig0)];
    int fine = var;

    for (int riga = p - 1; riga >= 1; --riga) {
      ++fine;
      double seg_max = std::numeric_limits<double>::quiet_NaN();
      bool seg_has_value = false;

      for (int k = var; k <= fine; ++k) {
        int orig_idx = p - ((k - 1) % p);
        double v = pval_matrix[idx(riga, orig_idx)];
        if (!std::isnan(v)) {
          if (!seg_has_value) {
            seg_max = v;
            seg_has_value = true;
          } else if (v > seg_max) {
            seg_max = v;
          }
        }
      }

      if (seg_has_value) {
        if (std::isnan(pval_var)) {
          pval_var = seg_max;
        } else if (seg_max > pval_var) {
          pval_var = seg_max;
        }
      }
      // if seg_has_value is false, pval_var unchanged (matches na.rm = TRUE
      // behaviour)
      corrected[idx(riga, var)] = pval_var;
    }
  }

  // Return corrected[, p:1] (reverse columns)
  cpp11::writable::doubles out(p * p);
  for (int col = 1; col <= p; ++col) {
    int src_col = p - (col - 1);
    for (int row = 1; row <= p; ++row) {
      out[idx(row, col)] = corrected[idx(row, src_col)];
    }
  }

  out.attr("dim") = dim;
  return out;
}