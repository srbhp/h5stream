/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "H5Cpp.h"
#include <iostream>
#include <string>
#include <vector>
//#include <vector>

// Maps C++ type to HDF5 type
template <typename T> inline const H5::PredType &get_datatype_for_hdf5();

// Reference:
// https://www.hdfgroup.org/HDF5/doc/cpplus_RM/class_h5_1_1_pred_type.html

template <> inline const H5::PredType &get_datatype_for_hdf5<char>() {
  return H5::PredType::NATIVE_CHAR;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<unsigned char>() {
  return H5::PredType::NATIVE_UCHAR;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<short>() {
  return H5::PredType::NATIVE_SHORT;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<unsigned short>() {
  return H5::PredType::NATIVE_USHORT;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<int>() {
  return H5::PredType::NATIVE_INT;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<unsigned int>() {
  return H5::PredType::NATIVE_UINT;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<long>() {
  return H5::PredType::NATIVE_LONG;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<unsigned long>() {
  return H5::PredType::NATIVE_ULONG;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<long long>() {
  return H5::PredType::NATIVE_LLONG;
}
template <>
inline const H5::PredType &get_datatype_for_hdf5<unsigned long long>() {
  return H5::PredType::NATIVE_ULLONG;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<float>() {
  return H5::PredType::NATIVE_FLOAT;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<double>() {
  return H5::PredType::NATIVE_DOUBLE;
}
template <> inline const H5::PredType &get_datatype_for_hdf5<long double>() {
  return H5::PredType::NATIVE_LDOUBLE;
}
//---------------------------------------------------------
namespace h5stream {
class dspace {
public:
  H5::DataSet dataset;
  explicit dspace(H5::DataSet datasetx) : dataset(datasetx){};
  template <typename T>
  void write_atr(const T data, const H5std_string &dataname) {
    auto type = get_datatype_for_hdf5<T>();
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    H5::Attribute attribute =
        dataset.createAttribute(dataname, type, attr_dataspace);
    attribute.write(type, &data);
  }
  template <typename T> void read_atr(T &data, const H5std_string &dataname) {
    // auto type = get_datatype_for_hdf5<T>();
    H5::Attribute attribute = dataset.openAttribute(dataname);
    H5::DataType type = attribute.getDataType();
    attribute.read(type, &data);
  }
};
} // namespace h5stream
//----------------------------------------------------------
namespace h5stream {
class h5stream {
public:
  H5std_string hdf5FileName;
  H5::H5File hdf5File;

  h5stream();
  h5stream(const std::string &fileName, const std::string &rw);
  h5stream(const std::string &fileName);
  void setFileName(const H5std_string &fileName);
  void setFileName(const H5std_string &fileName, const std::string &rw);

  template <typename T = double, template <typename...> class vec = std::vector>
  void write(const vec<T> &data, const H5std_string &datasetName);

  template <typename T = double, template <typename...> class vec = std::vector>
  void write(const H5std_string &datasetName, const vec<T> &data);
  // Write raw pointer
  template <typename T = double>
  void write(const H5std_string &datasetName, T *data, unsigned data_size);

  template <typename T = double, template <typename...> class vec = std::vector>
  vec<T> read(const H5std_string &datasetName);

  template <typename T = double, template <typename...> class vec = std::vector>
  void read(vec<T> &data, const H5std_string &datasetName) {
    data = read<T, vec>(datasetName);
  }
  void close() { hdf5File.close(); }

  double file_size() { return hdf5File.getFileSize() / (1024 * 1024.); }
  dspace get_dataspace(const H5std_string &dataset_name) {
    return dspace(hdf5File.openDataSet(dataset_name));
  }
  H5::Group create_group(const H5std_string &group_name) {
    return hdf5File.createGroup(group_name);
  }
};
} // namespace h5stream
template <typename T, template <typename...> class vec>
void h5stream::h5stream::write(const H5std_string &datasetName,
                               const vec<T> &data) {
  write(data, datasetName);
}

void h5stream::h5stream::setFileName(const H5std_string &fileName,
                                     const std::string &rw) {
  hdf5FileName = fileName;
  try {
    H5::Exception::dontPrint();
    if (rw.compare("r") == 0)
      hdf5File = H5::H5File(fileName, H5F_ACC_RDONLY);
    if (rw == "rw")
      hdf5File = H5::H5File(fileName, H5F_ACC_RDWR);
    if (rw == "x")
      hdf5File = H5::H5File(fileName, H5F_ACC_EXCL);
    if (rw == "tr")
      hdf5File = H5::H5File(fileName, H5F_ACC_TRUNC);
  } catch (...) {
    std::cout << " Error ::FileIException! (1 ) " << fileName << std::endl;
  }
}

void h5stream::h5stream::setFileName(const H5std_string &fileName) {
  setFileName(fileName, "tr");
}
h5stream::h5stream::h5stream() {}
h5stream::h5stream::h5stream(const std::string &fileName,
                             const std::string &rw) {
  setFileName(fileName, rw);
}
h5stream::h5stream::h5stream(const std::string &fileName) {
  setFileName(fileName, "tr");
}
template <typename T>
void h5stream::h5stream::write(const H5std_string &datasetName, T *data,
                               unsigned data_size) {
  try {
    H5::Exception::dontPrint();
    const int RANK = 1;
    auto type = get_datatype_for_hdf5<T>();
    hsize_t dimsf[1];     // dataset dimensions
    dimsf[0] = data_size; //
    H5::DataSpace dataspace(RANK, dimsf);
    H5::DataSet dataset = hdf5File.createDataSet(datasetName, type, dataspace);
    dataset.write(data, type);
  } catch (const H5::FileIException &error) {
    // error.printErrorStack();
    std::cout << " Error ::FileIException! (Write )" << std::endl;
  } catch (const H5::DataSetIException &error) {
    //    error.printErrorStack();
    std::cout << " Error ::DataSetIException ! " << std::endl;
  }
}

template <typename T, template <typename...> class vec>
void h5stream::h5stream::write(const vec<T> &data,
                               const H5std_string &datasetName) {
  try {
    H5::Exception::dontPrint();
    const int RANK = 1;
    auto type = get_datatype_for_hdf5<T>();
    hsize_t dimsf[1]; // dataset dimensions
    dimsf[0] = data.size();
    H5::DataSpace dataspace(RANK, dimsf);
    H5::DataSet dataset = hdf5File.createDataSet(datasetName, type, dataspace);
    dataset.write(data.data(), type);
  }

  catch (const H5::FileIException &error) {
    // error.printErrorStack();
    std::cout << " Error ::FileIException! (Write )" << std::endl;
  }

  catch (const H5::DataSetIException &error) {
    //    error.printErrorStack();
    std::cout << " Error ::DataSetIException ! " << std::endl;
  }
}

template <typename T, template <typename...> class vec>
vec<T> h5stream::h5stream::read(const H5std_string &dataset_name) {
  vec<T> data;
  try {
    H5::Exception::dontPrint();
    auto type = get_datatype_for_hdf5<T>();
    H5::DataSet dataset = hdf5File.openDataSet(dataset_name);
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t dim[1];
    dataspace.getSimpleExtentDims(dim, NULL);
    data.resize(dim[0]);
    dataset.read(data.data(), type, dataspace, dataspace);
  } catch (...) {
    //    H5::FileIException::printErrorStack();
    std::cout << " Error ::FileIException! (Read) " << std::endl;
  }

  return data;
}
