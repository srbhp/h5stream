#include "h5stream.hpp"
#include <iostream>
#include <string>
#include <vector>

int main()
{
  // Create File
  h5stream::h5stream file("sample.h5", "tr");

  // Create a vector
  std::vector<double> matrix { 1, 2, 3282, 932 };
  // Write to the file
  file.create_group("/data");
  file.write<double, std::vector>(
      matrix, "/data/matrix"); // Call write_vector
                               //  file.write<double, std::vector>(matrix,
                               //  "/matrix/Mat"); // Call write_vector
  // Write Attributes( Metadata) to the to the same data space
  auto dspace = file.get_dataspace("/data/matrix");
  dspace.write_atr<double>(1.2, "Units");

  // Read data from the file
  auto xx = file.read<double, std::vector>("/data/matrix");
  // OR
  // file.read<double, std::vector>(xx, "matrix");
  // Read Attribute
  double x = 0;
  dspace.read_atr<double>(x, "Units");
  std::cout << "Attribute : " << x << std::endl;
  std::cout << "HDF file size (MB): " << file.file_size() << std::endl;

  return 0;
}
