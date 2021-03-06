# [h5stream](https://github.com/srbhp/h5stream)
C++ Header only library for simple HDF5 input/output

## How to use 

Just include the `h5stream.hpp` into your your main file.

### Compile 

```
g++ -lhdf5 -lhdf5_cpp -std=c++1z example.cpp
```

### Example

#### Create File with a mode.


- "tr":   Create file, truncate if exists, Default
- "r":    Readonly, file must exist
- "rw": Read/write, file must exist
- "x":   Create file, fail if exists


```
h5stream::h5stream file("sample.h5", "tr");
// or 
h5stream::h5stream file("sample.h5");
```

#### write and read `std::vector` 

Create a vector and write it to the file


```
std::vector<double> matrix { 1, 2, 3282, 932 };
file.write<double>(matrix, "matrix"); 
```


#### write and read Metadata

Write Attributes( Metadata) to the to the same data space


```
auto dspace = file.get_dataspace("matrix");
dspace.write_atr<double>(1.2, "Units");
```


#### Read data from the file


```
auto xx = file.read_vector<double>("matrix");
//OR
file.read<double>(xx, "matrix");
```


#### Read Attribute (Metadata)
```
double x = 0;
dspace.read_atr<double>(x, "Units");
std::cout << "Attribute : " << x << std::endl;
std::cout << "HDF file size (MB): " << file.file_size() << std::endl;
```

