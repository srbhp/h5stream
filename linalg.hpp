#pragma once
#include <cmath>
#include <complex>
#include <iostream>

#include <cblas.h>
#include <lapacke.h>
// Uncomment bellow lines if you want to use
// MKL
//#define MKL_Complex16 std::complex<double>
//#include <mkl.h>
//#include <mkl_lapacke.h>
#include <string>
#include <tuple>
#include <vector>
typedef std::vector<double> vec;
namespace linalg {
vec dot(const vec &a, const vec &b, double alpha = 1.);
std::vector<double> kroneckerDot(const vec &a, const vec &b,
                                 double alpha = 1.0);
vec id_matrix(int n);
void disp_matrix(const std::string &st = " ", const vec &mat = vec(0));
void disp_array(const std::string &st = " ", const vec &a = vec(0));
std::tuple<vec, vec, vec, vec> nonsys_diag(vec &a);
void unitary_transform(vec &a, const vec &U);
double sum(const vec &a) {
  double aa{0};
  for (auto &a1 : a) {
    aa += a1;
  }
  return aa;
};
} // namespace linalg
void linalg::unitary_transform(vec &a, const vec &eigen_vector) {
  if (a.size() != eigen_vector.size()) {
    std::cout << "Error :  Matrix sizes are different of Uni.Tran.!"
              << std::endl;
  }
  unsigned nsize = std::sqrt(eigen_vector.size());
  vec rho_final(eigen_vector.size());
  vec rho_tmp(eigen_vector.size());

  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nsize, nsize, nsize, 1.0,
              eigen_vector.data(), nsize, a.data(), nsize, 0, rho_tmp.data(),
              nsize);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nsize, nsize, nsize,
              1.0, rho_tmp.data(), nsize, eigen_vector.data(), nsize, 0,
              rho_final.data(), nsize);

  a = std::move(rho_final);
}

vec linalg::dot(const vec &a, const vec &b, double alpha) {
  // This is done for square matrix
  // TODO: call lapack for large matrix
  if (a.size() != b.size()) {
    std::cout << "Matrix dont have same size !" << std::endl;
  }
  unsigned int n = std::sqrt(a.size());
  vec c(n * n, 0);
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < n; j++) {
      for (unsigned int k = 0; k < n; k++) {
        c[i * n + k] += alpha * a[i * n + j] * b[j * n + k];
      }
    }
  }
  return c;
}

std::tuple<vec, vec, vec, vec> linalg::nonsys_diag(vec &a) {
  auto n = static_cast<int>(std::sqrt(a.size()));
  vec wr(n, 0);
  vec wi(n, 0);
  vec vl(n * n, 0);
  vec vr(n * n, 0);
  auto info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'V', 'V', n, a.data(), n,
                            wr.data(), wi.data(), vl.data(), n, vr.data(), n);
  if (info > 0) {
    std::cout << "The algorithm failed to compute eigenvalues.\n";
    exit(1);
  }
  return {vl, vr, wr, wi};
}
vec linalg::id_matrix(int n) {
  vec id(n * n, 0);
  for (int i = 0; i < n; i++) {
    id[i * n + i] = 1.0;
  }
  return id;
}
vec linalg::add(const vec &a, const vec &b) {
  vec c(a.size());
  if (a.size() == b.size()) {
    for (unsigned i = 0; i < a.size(); i++) {
      c[i] = a[i] + b[i];
    }
  } else {
    // TODO use std::optional
    std::cout << "Error: Matrix sizes are different! A:" << a.size()
              << " B: " << b.size() << std::endl;
  }
  return c;
}
std::vector<double> linalg::kroneckerDot(const std::vector<double> &a,
                                         const std::vector<double> &b,
                                         double alpha) {
  // https://en.wikipedia.org/wiki/Kronecker_product
  long long unsigned m = std::sqrt(a.size());
  long long unsigned p = std::sqrt(b.size());
  long long unsigned dim = m * p;
  std::vector<double> c(dim * dim, 0);

  for (unsigned r = 0; r < m; r++) {
    for (unsigned s = 0; s < m; s++) {
      for (unsigned v = 0; v < p; v++) {
        for (unsigned w = 0; w < p; w++) {
          c[(r * p + v) * dim + (s * p + w)] +=
              a[r * m + s] * b[v * p + w] * alpha;
        }
      }
    }
  }
  return c;
}

void linalg::disp_array(const std::string &st, const vec &a) {
  std::cout << st << std::endl;
  unsigned n = a.size();
  for (unsigned i = 0; i < n; i++) {
    std::cout << a[i] << " ";
  }
  std::cout << std::endl;
}
void linalg::disp_matrix(const std::string &st, const vec &a) {
  std::cout << st << std::endl;
  unsigned n = static_cast<int>(std::sqrt(a.size()));
  for (unsigned i = 0; i < n; i++) {
    for (unsigned j = 0; j < n; j++) {
      std::cout << a[i * n + j] << " ";
    }
    std::cout << std::endl;
  }
}
