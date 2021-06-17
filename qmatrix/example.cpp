#include "qmatrix.hpp"

int main() {
  qmatrix<double> mat(
      3, 3,
      2.2); // This defines a 3x3 matrix wiill all the elements is set to zero
  mat - 1.2;
  std::cout << mat;
  auto mat2 = mat + 2;
  std::cout << mat2 << std::endl;
  return 0;
}
