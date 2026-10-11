#include "io/output.hpp"

#include <atomic>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <unistd.h>

namespace cppwrf::io {

namespace {

std::atomic<uint64_t> g_temp_counter{0};

std::string make_temp_path(const std::filesystem::path& target) {
  auto dir = target.parent_path();
  std::string filename = target.filename().string();
  pid_t pid = ::getpid();
  uint64_t count = ++g_temp_counter;
  std::string temp_name = filename + ".tmp." + std::to_string(pid) + "_" + std::to_string(count);
  return (dir.empty() ? std::filesystem::path(temp_name) : dir / temp_name).string();
}

} // namespace

// ==================== JsonWriter ====================

JsonWriter::JsonWriter() = default;

void JsonWriter::prepare_value() {
  if (stack_.empty()) {
    if (!buf_.empty()) {
      throw std::runtime_error("JsonWriter: cannot write multiple root values");
    }
    return;
  }
  auto& top = stack_.back();
  if (top.type == Frame::Object) {
    if (!top.in_value) {
      throw std::runtime_error("JsonWriter: value in object without a key");
    }
    top.in_value = false;
  } else { // Array
    if (top.has_items) {
      buf_.push_back(',');
    }
    top.has_items = true;
  }
}

void JsonWriter::start_object() {
  prepare_value();
  buf_.push_back('{');
  stack_.push_back({Frame::Object, false, false});
}

void JsonWriter::end_object() {
  if (stack_.empty() || stack_.back().type != Frame::Object || stack_.back().in_value) {
    throw std::runtime_error("JsonWriter: mismatched end_object");
  }
  stack_.pop_back();
  buf_.push_back('}');
}

void JsonWriter::start_array() {
  prepare_value();
  buf_.push_back('[');
  stack_.push_back({Frame::Array, false, false});
}

void JsonWriter::end_array() {
  if (stack_.empty() || stack_.back().type != Frame::Array) {
    throw std::runtime_error("JsonWriter: mismatched end_array");
  }
  stack_.pop_back();
  buf_.push_back(']');
}

void JsonWriter::key(std::string_view k) {
  if (stack_.empty() || stack_.back().type != Frame::Object || stack_.back().in_value) {
    throw std::runtime_error("JsonWriter: key called in invalid state");
  }
  auto& top = stack_.back();
  if (top.has_items) {
    buf_.push_back(',');
  }
  top.has_items = true;
  buf_ += escape_string(k);
  buf_.push_back(':');
  top.in_value = true;
}

void JsonWriter::write_null() {
  prepare_value();
  buf_ += "null";
}

void JsonWriter::write_bool(bool b) {
  prepare_value();
  buf_ += (b ? "true" : "false");
}

void JsonWriter::write_int(long long v) {
  prepare_value();
  char b[32];
  auto [ptr, ec] = std::to_chars(b, b + sizeof(b), v);
  buf_.append(b, ptr - b);
}

void JsonWriter::write_double(double v) {
  if (std::isnan(v) || std::isinf(v)) {
    throw std::runtime_error("JsonWriter: NaN or Infinity cannot be written to JSON");
  }
  prepare_value();
  char b[64];
  auto [ptr, ec] = std::to_chars(b, b + sizeof(b), v, std::chars_format::general, std::numeric_limits<double>::max_digits10);
  buf_.append(b, ptr - b);
}

void JsonWriter::write_string(std::string_view s) {
  prepare_value();
  buf_ += escape_string(s);
}

void JsonWriter::key_value(std::string_view k, std::string_view v) {
  key(k);
  write_string(v);
}

void JsonWriter::key_value(std::string_view k, long long v) {
  key(k);
  write_int(v);
}

void JsonWriter::key_value(std::string_view k, double v) {
  key(k);
  write_double(v);
}

void JsonWriter::key_value(std::string_view k, bool v) {
  key(k);
  write_bool(v);
}

std::string JsonWriter::escape_string(std::string_view s) {
  std::string res;
  res.reserve(s.size() + 2);
  res.push_back('"');
  for (unsigned char c : s) {
    switch (c) {
      case '"': res += "\\\""; break;
      case '\\': res += "\\\\"; break;
      case '\b': res += "\\b"; break;
      case '\f': res += "\\f"; break;
      case '\n': res += "\\n"; break;
      case '\r': res += "\\r"; break;
      case '\t': res += "\\t"; break;
      default:
        if (c < 0x20) {
          char hex[8];
          std::snprintf(hex, sizeof(hex), "\\u%04x", c);
          res += hex;
        } else {
          res.push_back(static_cast<char>(c));
        }
        break;
    }
  }
  res.push_back('"');
  return res;
}

std::string JsonWriter::str() const {
  if (!stack_.empty()) {
    throw std::runtime_error("JsonWriter: unclosed containers in JSON");
  }
  return buf_;
}

// ==================== write_text_file_atomic ====================

void write_text_file_atomic(const std::string& path, const std::string& content) {
  std::filesystem::path target(path);
  if (target.has_parent_path()) {
    std::filesystem::create_directories(target.parent_path());
  }

  std::filesystem::path tmp_path = make_temp_path(target);
  {
    std::ofstream out(tmp_path, std::ios::out | std::ios::trunc);
    if (!out.is_open()) {
      throw std::runtime_error("write_text_file_atomic: failed to open temp file " + tmp_path.string());
    }
    out.write(content.data(), static_cast<std::streamsize>(content.size()));
    out.close();
    if (!out) {
      std::error_code ec;
      std::filesystem::remove(tmp_path, ec);
      throw std::runtime_error("write_text_file_atomic: failed to write to temp file " + tmp_path.string());
    }
  }

  std::error_code ec;
  std::filesystem::rename(tmp_path, target, ec);
  if (ec) {
    std::filesystem::remove(tmp_path, ec);
    throw std::runtime_error("write_text_file_atomic: failed to rename " + tmp_path.string() + " to " + target.string());
  }
}

// ==================== SnapshotWriter ====================

SnapshotWriter::SnapshotWriter(std::string out_dir) : out_dir_(std::move(out_dir)) {
  if (!out_dir_.empty()) {
    std::filesystem::create_directories(out_dir_);
  }
}

void SnapshotWriter::write(long step, double time, const std::vector<VarView>& vars) const {
  char step_buf[32];
  std::snprintf(step_buf, sizeof(step_buf), "snap_%06ld", step);
  std::string base_name = step_buf;

  std::filesystem::path dir(out_dir_);
  if (!dir.empty()) {
    std::filesystem::create_directories(dir);
  }

  std::filesystem::path bin_path = dir / (base_name + ".bin");
  std::filesystem::path json_path = dir / (base_name + ".json");

  // Write binary file atomically
  std::filesystem::path tmp_bin_path = make_temp_path(bin_path);
  struct VarMeta {
    std::string name;
    char stagger;
    int nx, ny, nz;
    std::string dtype_str;
    uint64_t offset;
    uint64_t nbytes;
  };
  std::vector<VarMeta> meta_list;
  uint64_t current_offset = 0;

  {
    std::ofstream bin_out(tmp_bin_path, std::ios::binary | std::ios::out | std::ios::trunc);
    if (!bin_out.is_open()) {
      throw std::runtime_error("SnapshotWriter: failed to open temp file " + tmp_bin_path.string());
    }

    for (const auto& var : vars) {
      size_t elem_size = (var.dtype == DType::F32 ? sizeof(float) : sizeof(double));
      std::string dtype_str = (var.dtype == DType::F32 ? "float32" : "float64");
      size_t n_elements = static_cast<size_t>(var.nx) * static_cast<size_t>(var.ny) * static_cast<size_t>(var.nz);
      size_t nbytes = n_elements * elem_size;

      if (n_elements > 0 && var.data == nullptr) {
        bin_out.close();
        std::error_code ec;
        std::filesystem::remove(tmp_bin_path, ec);
        throw std::invalid_argument("SnapshotWriter: null data pointer for variable " + var.name);
      }

      // Pack contiguous [nz][ny][nx]
      if (var.dtype == DType::F32) {
        const float* src = static_cast<const float*>(var.data);
        for (int k = 0; k < var.nz; ++k) {
          for (int j = 0; j < var.ny; ++j) {
            const float* row_ptr = src + (j * var.sj + k * var.sk);
            bin_out.write(reinterpret_cast<const char*>(row_ptr), static_cast<std::streamsize>(var.nx * sizeof(float)));
          }
        }
      } else {
        const double* src = static_cast<const double*>(var.data);
        for (int k = 0; k < var.nz; ++k) {
          for (int j = 0; j < var.ny; ++j) {
            const double* row_ptr = src + (j * var.sj + k * var.sk);
            bin_out.write(reinterpret_cast<const char*>(row_ptr), static_cast<std::streamsize>(var.nx * sizeof(double)));
          }
        }
      }

      meta_list.push_back({var.name, var.stagger, var.nx, var.ny, var.nz, dtype_str, current_offset, nbytes});
      current_offset += nbytes;
    }

    bin_out.close();
    if (!bin_out) {
      std::error_code ec;
      std::filesystem::remove(tmp_bin_path, ec);
      throw std::runtime_error("SnapshotWriter: failed to write binary data to " + tmp_bin_path.string());
    }
  }

  std::error_code ec;
  std::filesystem::rename(tmp_bin_path, bin_path, ec);
  if (ec) {
    std::filesystem::remove(tmp_bin_path, ec);
    throw std::runtime_error("SnapshotWriter: failed to rename " + tmp_bin_path.string() + " to " + bin_path.string());
  }

  // Build JSON metadata
  JsonWriter jw;
  jw.start_object();
  jw.key("step");
  jw.write_int(step);
  jw.key("time");
  jw.write_double(time);
  jw.key("variables");
  jw.start_array();
  for (const auto& vm : meta_list) {
    jw.start_object();
    jw.key("name");
    jw.write_string(vm.name);
    jw.key("stagger");
    jw.write_string(std::string(1, vm.stagger));
    jw.key("shape");
    jw.start_array();
    jw.write_int(vm.nz);
    jw.write_int(vm.ny);
    jw.write_int(vm.nx);
    jw.end_array();
    jw.key("dtype");
    jw.write_string(vm.dtype_str);
    jw.key("offset");
    jw.write_int(static_cast<long long>(vm.offset));
    jw.key("nbytes");
    jw.write_int(static_cast<long long>(vm.nbytes));
    jw.end_object();
  }
  jw.end_array();
  jw.end_object();

  write_text_file_atomic(json_path.string(), jw.str());
}

// ==================== DiagWriter ====================

DiagWriter::DiagWriter(const std::string& path, std::vector<std::string> columns)
    : path_(path), columns_(std::move(columns)) {
  std::filesystem::path target(path_);
  if (target.has_parent_path()) {
    std::filesystem::create_directories(target.parent_path());
  }

  file_.open(path_, std::ios::out | std::ios::trunc);
  if (!file_.is_open()) {
    throw std::runtime_error("DiagWriter: failed to open " + path_);
  }

  for (size_t i = 0; i < columns_.size(); ++i) {
    if (i > 0) file_ << ",";
    file_ << columns_[i];
  }
  file_ << "\n";
  file_.flush();
}

DiagWriter::~DiagWriter() {
  if (file_.is_open()) {
    file_.close();
  }
}

void DiagWriter::row(const std::vector<double>& values) {
  if (values.size() != columns_.size()) {
    throw std::runtime_error("DiagWriter: row size mismatch (got " + std::to_string(values.size()) +
                             ", expected " + std::to_string(columns_.size()) + ")");
  }

  for (size_t i = 0; i < values.size(); ++i) {
    if (i > 0) file_ << ",";
    char buf[64];
    auto [ptr, ec] = std::to_chars(buf, buf + sizeof(buf), values[i], std::chars_format::general, 9);
    file_.write(buf, ptr - buf);
  }
  file_ << "\n";
  file_.flush();
}

} // namespace cppwrf::io
