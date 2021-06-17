#include "h5stream.hpp"

#include "testfn.hpp"
void testfn() {

  h5stream::h5stream file("sample.h5", "rw");
  std::vector<double> matrix{1, 2, 3282, 932};

  // Write to the file
  file.write<double>("/data/testMatrix", matrix);
}
