#include "io/output.hpp"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#define ASSERT_TRUE(cond) \
  do { \
    if (!(cond)) { \
      std::cerr << "ASSERTION FAILED (" #cond ") at " << __FILE__ << ":" << __LINE__ << "\n"; \
      std::exit(1); \
    } \
  } while (false)

#define ASSERT_FALSE(cond) ASSERT_TRUE(!(cond))

#define ASSERT_EQ(a, b) \
  do { \
    if ((a) != (b)) { \
      std::cerr << "ASSERTION FAILED (" #a " == " #b ") [" << (a) << " != " << (b) << "] at " \
                << __FILE__ << ":" << __LINE__ << "\n"; \
      std::exit(1); \
    } \
  } while (false)

#define ASSERT_THROWS(expr, exc_type) \
  do { \
    bool caught = false; \
    try { \
      expr; \
    } catch (const exc_type&) { \
      caught = true; \
    } catch (...) {} \
    if (!caught) { \
      std::cerr << "Expected exception " #exc_type " was not thrown at " << __FILE__ << ":" << __LINE__ << "\n"; \
      std::exit(1); \
    } \
  } while (false)

using namespace cppwrf::io;

void test_json_writer() {
  // Test escaping
  std::string escaped = JsonWriter::escape_string("quote: \" backslash: \\ newline: \n tab: \t ctrl: \x01");
  ASSERT_EQ(escaped, "\"quote: \\\" backslash: \\\\ newline: \\n tab: \\t ctrl: \\u0001\"");

  // Test standard structure
  JsonWriter jw;
  jw.start_object();
  jw.key("name");
  jw.write_string("cppwrf");
  jw.key("count");
  jw.write_int(42);
  jw.key("pi");
  jw.write_double(3.141592653589793);
  jw.key("valid");
  jw.write_bool(true);
  jw.key("none");
  jw.write_null();
  jw.key("list");
  jw.start_array();
  jw.write_int(1);
  jw.write_int(2);
  jw.end_array();
  jw.end_object();

  std::string s = jw.str();
  ASSERT_TRUE(s.find("\"name\":\"cppwrf\"") != std::string::npos);
  ASSERT_TRUE(s.find("\"count\":42") != std::string::npos);
  ASSERT_TRUE(s.find("\"valid\":true") != std::string::npos);
  ASSERT_TRUE(s.find("\"none\":null") != std::string::npos);
  ASSERT_TRUE(s.find("\"list\":[1,2]") != std::string::npos);

  // NaN in JSON -> error
  JsonWriter jw_nan;
  jw_nan.start_object();
  jw_nan.key("val");
  ASSERT_THROWS(jw_nan.write_double(std::numeric_limits<double>::quiet_NaN()), std::runtime_error);

  // Inf in JSON -> error
  JsonWriter jw_inf;
  jw_inf.start_object();
  jw_inf.key("val");
  ASSERT_THROWS(jw_inf.write_double(std::numeric_limits<double>::infinity()), std::runtime_error);
  ASSERT_THROWS(jw_inf.write_double(-std::numeric_limits<double>::infinity()), std::runtime_error);
}

void test_diag_writer() {
  std::filesystem::path test_dir = std::filesystem::temp_directory_path() / "cppwrf_test_diag";
  std::filesystem::create_directories(test_dir);
  std::filesystem::path csv_path = test_dir / "diag.csv";

  {
    DiagWriter dw(csv_path.string(), {"step", "time", "energy", "temp"});
    dw.row({0.0, 0.0, 1.23456789012, 300.15});
    dw.row({1.0, 0.1, 1.23456789099, 300.25});
    ASSERT_THROWS(dw.row({2.0, 0.2}), std::runtime_error); // size mismatch
  }

  // Read back
  std::ifstream in(csv_path);
  ASSERT_TRUE(in.is_open());
  std::string header, r1, r2;
  std::getline(in, header);
  std::getline(in, r1);
  std::getline(in, r2);

  ASSERT_EQ(header, "step,time,energy,temp");
  // Check %.9g precision for 1.23456789012 -> 1.23456789
  ASSERT_TRUE(r1.find("1.23456789") != std::string::npos);
  ASSERT_TRUE(r1.find("300.15") != std::string::npos);

  std::filesystem::remove_all(test_dir);
}

void test_write_text_file_atomic() {
  std::filesystem::path test_file = std::filesystem::temp_directory_path() / "cppwrf_test_atomic" / "sub" / "file.txt";
  write_text_file_atomic(test_file.string(), "atomic content 123");

  std::ifstream in(test_file);
  ASSERT_TRUE(in.is_open());
  std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  ASSERT_EQ(content, "atomic content 123");

  // Overwrite
  write_text_file_atomic(test_file.string(), "updated content");
  std::ifstream in2(test_file);
  std::string content2((std::istreambuf_iterator<char>(in2)), std::istreambuf_iterator<char>());
  ASSERT_EQ(content2, "updated content");

  std::filesystem::remove_all(std::filesystem::temp_directory_path() / "cppwrf_test_atomic");
}

void test_snapshot_roundtrip() {
  std::filesystem::path out_dir = std::filesystem::temp_directory_path() / "cppwrf_test_snapshot";
  std::filesystem::remove_all(out_dir);

  SnapshotWriter writer(out_dir.string());

  // Set up strided 3D array with halos
  // Logical: nx = 4, ny = 3, nz = 2
  // Halo: 2 in x, 2 in y, 1 in z
  const int nx = 4, ny = 3, nz = 2;
  const int halo_x = 2, halo_y = 2, halo_z = 1;
  const int alloc_nx = nx + 2 * halo_x; // 8
  const int alloc_ny = ny + 2 * halo_y; // 7
  const int alloc_nz = nz + 2 * halo_z; // 4

  const std::ptrdiff_t sj = alloc_nx;                 // 8
  const std::ptrdiff_t sk = alloc_nx * alloc_ny;      // 56
  const size_t total_alloc = sk * alloc_nz;           // 224

  // Buffer 1: F32 field
  std::vector<float> buf_f32(total_alloc, -999.0f); // fill halo with garbage
  float* origin_f32 = buf_f32.data() + (halo_x + halo_y * sj + halo_z * sk);

  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        origin_f32[i + j * sj + k * sk] = static_cast<float>(100.0f + i + 10.0f * j + 100.0f * k);
      }
    }
  }

  // Buffer 2: F64 field
  std::vector<double> buf_f64(total_alloc, -999.0);
  double* origin_f64 = buf_f64.data() + (halo_x + halo_y * sj + halo_z * sk);

  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        origin_f64[i + j * sj + k * sk] = static_cast<double>(200.0 + i * 0.1 + j * 0.01 + k * 0.001);
      }
    }
  }

  // Buffer 3: 2D F32 field (nz = 1)
  const int nz_2d = 1;
  const size_t total_alloc_2d = sj * alloc_ny;
  std::vector<float> buf_2d(total_alloc_2d, -999.0f);
  float* origin_2d = buf_2d.data() + (halo_x + halo_y * sj);
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      origin_2d[i + j * sj] = static_cast<float>(50.0f + i + 5.0f * j);
    }
  }

  std::vector<VarView> vars = {
    {"u", 'U', nx, ny, nz, origin_f32, DType::F32, sj, sk},
    {"temp", 'M', nx, ny, nz, origin_f64, DType::F64, sj, sk},
    {"psfc", 'M', nx, ny, nz_2d, origin_2d, DType::F32, sj, 0}
  };

  writer.write(42, 10.5, vars);

  // Verify snap_000042.json and snap_000042.bin exist
  std::filesystem::path json_file = out_dir / "snap_000042.json";
  std::filesystem::path bin_file = out_dir / "snap_000042.bin";

  ASSERT_TRUE(std::filesystem::exists(json_file));
  ASSERT_TRUE(std::filesystem::exists(bin_file));

  // Check binary file size
  // u: 4 * 3 * 2 * 4 = 96 bytes
  // temp: 4 * 3 * 2 * 8 = 192 bytes
  // psfc: 4 * 3 * 1 * 4 = 48 bytes
  // total = 336 bytes
  size_t bin_size = std::filesystem::file_size(bin_file);
  ASSERT_EQ(bin_size, 336UL);

  // Read .bin back and compare bit-exactly
  std::ifstream bin_in(bin_file, std::ios::binary);
  ASSERT_TRUE(bin_in.is_open());
  std::vector<char> raw_bytes(bin_size);
  bin_in.read(raw_bytes.data(), bin_size);
  ASSERT_EQ(bin_in.gcount(), static_cast<std::streamsize>(bin_size));

  // Verify variable 1 (u): offset 0, 96 bytes
  const float* read_u = reinterpret_cast<const float*>(raw_bytes.data() + 0);
  size_t u_idx = 0;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        float expected = origin_f32[i + j * sj + k * sk];
        float actual = read_u[u_idx++];
        ASSERT_TRUE(std::memcmp(&expected, &actual, sizeof(float)) == 0);
      }
    }
  }

  // Verify variable 2 (temp): offset 96, 192 bytes
  const double* read_temp = reinterpret_cast<const double*>(raw_bytes.data() + 96);
  size_t temp_idx = 0;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        double expected = origin_f64[i + j * sj + k * sk];
        double actual = read_temp[temp_idx++];
        ASSERT_TRUE(std::memcmp(&expected, &actual, sizeof(double)) == 0);
      }
    }
  }

  // Verify variable 3 (psfc): offset 288, 48 bytes
  const float* read_psfc = reinterpret_cast<const float*>(raw_bytes.data() + 288);
  size_t psfc_idx = 0;
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      float expected = origin_2d[i + j * sj];
      float actual = read_psfc[psfc_idx++];
      ASSERT_TRUE(std::memcmp(&expected, &actual, sizeof(float)) == 0);
    }
  }

  std::filesystem::remove_all(out_dir);
}

int main() {
  std::cout << "Running test_output...\n";
  test_json_writer();
  test_diag_writer();
  test_write_text_file_atomic();
  test_snapshot_roundtrip();
  std::cout << "test_output: ALL TESTS PASSED\n";
  return 0;
}
