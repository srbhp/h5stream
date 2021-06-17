#pragma once
#include <cblas.h>
#include <complex>
#include <initializer_list>
#include <iostream>
#include <lapacke.h>
// #include <mkl.h>
// #include <mkl_lapacke.h>
#include <omp.h>
#include <ostream>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

template <typename NType1>
std::ostream &operator<<(std::ostream &out, const std::vector<NType1> &val) {
  for (size_t i = 0; i < val.size(); ++i) {
    out << val[i] << " ";
  }
  out << "\n";
  return out;
}

template <typename NType1>
std::ostream &operator<<(std::ostream &out,
                         const std::vector<std::vector<NType1>> &val) {
  out << "\n";
  for (auto x : val) {
    for (auto y : x) {
      out << y << " ";
    }
    out << "\n";
  }
  return out;
}
template <typename T> void LOGGER(std::string name, T x) {
  std::cout << "## " << name << " : " << x << std::endl;
}
#define logger(name) LOGGER(#name, (name))
// This is a matrix class quamtum(q) operators
// This is also consistent with Lapack_row_MAJOR
using cm_vec = std::vector<std::complex<double>>;
template <class T = double> class qmatrix {
  std::vector<T> mat{0};

public:
  size_t row{0}, column{0}, dim{0};
  qmatrix(std::initializer_list<T> inputVec, size_t _row = 0,
          size_t _column = 0)
      : mat(inputVec) {
    // check if this is a quare matrix
    // std::cout << "constructed with a " << inputVec.size() << "-element
    // list\n";
    if (_row == 0 && _column == 0) {
      this->row = std::sqrt(inputVec.size());
      this->column = std::sqrt(inputVec.size());
      this->dim = inputVec.size();
      if (this->dim != this->row * this->column) {
        throw std::runtime_error(
            "initializer_list:vector is not a square matrix! ");
      }
    } else {
      this->row = _row;
      this->column = _column;
      this->dim = _column * _row;
    }
    this->mat.shrink_to_fit();
  }

  qmatrix(std::vector<T> &inputVec, size_t _row = 0, size_t _column = 0) {
    // check if this is a quare matrix
    if (_row == 0 && _column == 0) {
      this->row = std::sqrt(inputVec.size());
      this->column = std::sqrt(inputVec.size());
      this->dim = inputVec.size();
      if (this->dim != this->row * this->column) {
        throw std::runtime_error("vector is not a square matrix! ");
      }
    } else {
      this->row = _row;
      this->column = _column;
      this->dim = _column * _row;
    }
    this->mat = std::move(inputVec);
    this->mat.shrink_to_fit();
  }
  qmatrix(size_t _row, size_t _column, T populate) {
    this->row = _row;
    this->column = _column;
    this->dim = _column * _row;
    this->mat = std::vector<T>(_row * _column, populate);
    this->mat.shrink_to_fit();
    // TODO: restrict memory usage
  }
  // This is for square matrix
  qmatrix(size_t N, T populate) { qmatrix(N, N, populate); }
  qmatrix() { qmatrix(0, 0, 0); }
  [[nodiscard]] T &operator()(size_t i) { return this->mat[i]; }
  [[nodiscard]] T operator()(size_t i) const { return this->mat[i]; }
  [[nodiscard]] T &at(size_t i) { return this->mat[i]; }
  [[nodiscard]] T at(size_t i) const { return this->mat[i]; }

  [[nodiscard]] T &operator()(size_t i, size_t j) {
    return this->mat[i * column + j];
  }
  [[nodiscard]] T operator()(size_t i, size_t j) const {
    return this->mat[i * column + j];
  }
  // add a at operator
  [[nodiscard]] T at(size_t i, size_t j) { return this->mat[i * column + j]; }
  [[nodiscard]] T at(size_t i, size_t j) const {
    return this->mat[i * column + j];
  }

  // TODO: arithmetic operator
  // TODO: use stl
  [[nodiscard]] T trace() const {
    if (this->column == this->row) {
      T result{0};
      for (size_t i = 0; i < this->row; i++) {
        result += this->mat[i * column + i];
      }
      return result;
    } else {
      throw std::runtime_error("Matrix is not square matrix");
    }
  }
  [[nodiscard]] qmatrix<T> id(size_t _row = 0) const {
    qmatrix<T> result(_row, _row, 0); // reverse the row and column
    for (size_t i = 0; i < _row; i++) {
      result(i, i) = 1.0;
    }
    return result;
  }
  [[nodiscard]] qmatrix<double> real() const {
    qmatrix<double> result(this->column, this->row,
                           0); // reverse the row and column
    for (size_t i = 0; i < this->row; i++) {
      result(i) = this->at(i).real();
    }
    return result;
  }
  [[nodiscard]] qmatrix<double> imag() const {
    qmatrix<double> result(this->column, this->row,
                           0); // reverse the row and column
    for (size_t i = 0; i < this->row; i++) {
      result(i) = this->at(i).imag();
    }
    return result;
  }

  [[nodiscard]] qmatrix<T> transpose() const {
    qmatrix<T> result(this->column, this->row, 0); // reverse the row and column
    for (size_t i = 0; i < this->row; i++) {
      for (size_t j = 0; j < this->column; j++) {
        result(j, i) = this->mat[i * column + j];
      }
    }
    return result;
  }

  [[nodiscard]] qmatrix<T> operator+(const T &x) const {
    qmatrix result(this->row, this->column, 0);
#pragma omp parallel for
    for (size_t i = 0; i < this->dim; i++) {
      result(i) = this->mat.at(i) + x;
    }
    return result;
  }
  [[nodiscard]] qmatrix operator-(const T &x) const {
    qmatrix result(this->row, this->column, 0);
#pragma omp parallel for
    for (size_t i = 0; i < this->dim; i++) {
      result(i) = this->mat.at(i) - x;
    }
    return result;
  }
  [[nodiscard]] qmatrix operator*(const T &x) const {
    qmatrix result(this->row, this->column, 0);
#pragma omp parallel for
    for (size_t i = 0; i < this->dim; i++) {
      result(i) = this->mat.at(i) * x;
    }
    return result;
  }
  [[nodiscard]] qmatrix operator/(const T &x) const {
    qmatrix result(this->row, this->column, 0);
#pragma omp parallel for
    for (size_t i = 0; i < this->dim; i++) {
      result(i) = this->mat.at(i) / x;
    }
    return result;
  }
  //  void reset(const T &x = 0) {
  //    for (auto &aa : this->mat) {
  //      aa = x;
  //    }
  //  }
  [[nodiscard]] size_t size() const { return this->column * this->row; }
  auto data() { return this->mat.data(); }
  [[nodiscard]] auto begin() const { return this->mat.begin(); }
  [[nodiscard]] auto end() const { return this->mat.end(); }
  //
  //
  void display() {}
  [[nodiscard]] qmatrix operator+(const qmatrix<T> &rhs) const {
    // std::cout << "Started operator+";
    if (this->row == rhs.row && this->column == rhs.column) {
      qmatrix result(this->row, this->column, 0);
#pragma omp parallel for
      for (size_t i = 0; i < this->dim; i++) {
        result(i) = this->mat[i] + rhs(i);
      }
      // std::cout << "End of operator+" << std::endl;
      return result;
    } else {
      throw std::runtime_error("qmatrix have different size for operator+");
    }
  }
  // Substract
  [[nodiscard]] qmatrix operator-(const qmatrix<T> &rhs) const {
    if (this->row == rhs.row && this->column == rhs.column) {
      qmatrix result(this->row, this->column, 0);
#pragma omp parallel for
      for (size_t i = 0; i < this->column; i++) {
        result(i) = this->mat[i] - rhs(i);
      }
      return result;
    } else {
      throw std::runtime_error("qmatrix have different size for operator -\n");
    }
  }
  [[nodiscard]] qmatrix operator*(qmatrix<T> &rhs) { return this->dot(rhs); }
  [[nodiscard]] qmatrix dot(qmatrix<T> &rhs) {
    if (this->column == rhs.row) {
      qmatrix result(this->row, rhs.column, 0);
      int m = this->row;
      int k = this->column;
      int n = rhs.column;
      const double alpha = 1;
      const double beta = 0;
      if constexpr (std::is_same_v<T, double>) {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha,
                    this->data(), k, rhs.data(), n, beta, result.data(), n);
      }
      if constexpr (std::is_same_v<T, std::complex<double>>) {
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, &alpha,
                    this->data(), k, rhs.data(), n, &beta, result.data(), n);
      }
      return result;
    } else {
      throw std::runtime_error("dot:qmatrix have different size for dot\n");
    }
  }

  [[nodiscard]] std::vector<double> diag() {
    // This function diagonalizes a symmetric/harmitian matrix.
    // So eigen value are always real(double).
    // If the matrix is not symmetric the call nonsys_diag
    if (this->row != this->column) {
      throw std::runtime_error("Error: Matrix is not a square matrix! \n");
    }
    std::vector<double> w(this->row, 0);
    int n = w.size();

    int info = -1;
    if constexpr (std::is_same_v<T, double>) {
      info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', n, this->data(), n,
                            w.data());
    }
    if constexpr (std::is_same_v<T, std::complex<double>>) {
      info = LAPACKE_zheevd(
          LAPACK_ROW_MAJOR, 'V', 'U', n,
          reinterpret_cast<__complex__ double *>(this->data()), n, w.data());
    }
    // int info= LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', n, a, n, w );
    if (info > 0) {
      std::cout << "Error:Not able to solve Eigen value problem." << std::endl;
    }
    return w;
  }

  [[nodiscard]] std::tuple<qmatrix<T>, qmatrix<T>, cm_vec> nonsys_diag() {
    if (this->row != this->column) {
      throw std::runtime_error("Error: Matrix is not a square matrix! \n");
    }
    std::vector<T> w(this->row, 0);
    int n = w.size();
    qmatrix<T> vl(n, n, 0);
    qmatrix<T> vr(n, n, 0);

    if constexpr (std::is_same_v<T, std::complex<double>>) {
      auto info =
          LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'V', 'V', n,
                        // recast the complex pointer
                        reinterpret_cast<__complex__ double *>(this->data()), n,
                        // recast the complex pointer
                        reinterpret_cast<__complex__ double *>(w.data()),
                        // recast the complex pointer
                        reinterpret_cast<__complex__ double *>(vl.data()), n,
                        // recast the complex pointer
                        reinterpret_cast<__complex__ double *>(vr.data()), n);
      /* Check for convergence */
      if (info > 0) {
        throw std::runtime_error(
            "The algorithm LAPACKE_zgeev failed to compute eigenvalues.\n");
      }
    }
    return {vl, vr, w};
  }

  friend std::ostream &operator<<(std::ostream &out, const qmatrix<T> &val) {
    out << "\n";
    for (size_t i = 0; i < val.row; ++i) {
      for (size_t j = 0; j < val.column; ++j) {
        out << val(i, j) << " ";
      }
      out << "\n";
    }
    return out;
  }

  qmatrix<T> krDot(const qmatrix<T> &rhs, double alpha = 1) {
    // https://en.wikipedia.org/wiki/Kronecker_product
    size_t m = this->row;
    size_t n = this->column;
    size_t p = rhs.row;
    size_t q = rhs.column;
    qmatrix<T> result(this->row * rhs.row, this->column * rhs.column, 0);
#pragma omp parallel for
    for (size_t r = 0; r < m; r++) {
      for (size_t s = 0; s < n; s++) {
        for (size_t v = 0; v < p; v++) {
          for (size_t w = 0; w < q; w++) {
            result((r * p + v), (s * q + w)) =
                this->at(r, s) * rhs(v, w) * alpha;
          }
        }
      }
    }
    return result;
  }
};
