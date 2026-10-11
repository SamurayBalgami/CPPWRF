#ifndef CPPWRF_IO_OUTPUT_HPP
#define CPPWRF_IO_OUTPUT_HPP

#include <cstddef>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

namespace cppwrf::io {

enum class DType { F32, F64 };

struct VarView {
  std::string name;
  char stagger;            // 'M', 'U', 'V', 'W'
  int nx, ny, nz;          // logical extents of THIS variable (already staggered), no halo
  const void* data;        // host pointer to element (0,0,0)
  DType dtype;
  std::ptrdiff_t sj, sk;   // element strides for j and k; the i stride is always 1
};                         // element (i,j,k) = data[i + j*sj + k*sk]; 2D fields: nz = 1

class JsonWriter {
 public:
  JsonWriter();

  void start_object();
  void end_object();

  void start_array();
  void end_array();

  void key(std::string_view k);

  void write_null();
  void write_bool(bool b);
  void write_int(long long v);
  void write_double(double v);
  void write_string(std::string_view s);

  // Convenience methods
  void key_value(std::string_view k, std::string_view v);
  void key_value(std::string_view k, long long v);
  void key_value(std::string_view k, double v);
  void key_value(std::string_view k, bool v);

  static std::string escape_string(std::string_view s);

  std::string str() const;

 private:
  struct Frame {
    enum Type { Object, Array } type;
    bool has_items = false;
    bool in_value = false;
  };

  std::vector<Frame> stack_;
  std::string buf_;

  void prepare_value();
};

class SnapshotWriter {
 public:
  explicit SnapshotWriter(std::string out_dir);
  void write(long step, double time, const std::vector<VarView>& vars) const;

 private:
  std::string out_dir_;
};

class DiagWriter {
 public:
  DiagWriter(const std::string& path, std::vector<std::string> columns);
  ~DiagWriter();
  void row(const std::vector<double>& values);

 private:
  std::string path_;
  std::vector<std::string> columns_;
  std::ofstream file_;
};

void write_text_file_atomic(const std::string& path, const std::string& content);

} // namespace cppwrf::io

#endif // CPPWRF_IO_OUTPUT_HPP
