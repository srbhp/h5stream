#include "h5stream.hpp"

#include "testfn.hpp"
#include <array>
#include <iostream>
#include <string>
#include <vector>
int main() {
  // Create File
  h5stream::h5stream file("sample.h5", "tr");

  // Create a vector
  std::vector<double> matrix{1, 2, 3282, 932};

  // Write to the file
  file.create_group("/data");
  file.write<double>("/data/matrix", matrix);
  file.write<double>("/data/matrix2 ", matrix.data(),
                     matrix.size()); // Save raw pointer storage

  testfn();
  // Write Attributes( Metadata) to the to the same data space
  auto dspace = file.get_dataspace("/data/matrix");
  dspace.write_atr<double>(1.2, "Units");

  // Read data from the file into an std::vector or std::array
  std::vector<double> xx;
  file.read<double>(xx, "/data/matrix");
  // OR
  // file.read<double, std::vector>(xx, "matrix");
  // Read Attribute
  double x = 0;
  dspace.read_atr<double>(x, "Units");
  std::cout << "Attribute : " << x << std::endl;
  std::cout << "HDF file size (MB): " << file.file_size() << std::endl;

  return 0;
}
