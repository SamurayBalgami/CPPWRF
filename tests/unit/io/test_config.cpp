#include "io/config.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
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

#define ASSERT_THROWS_CONTAINS(expr, exc_type, substr) \
  do { \
    bool caught = false; \
    try { \
      expr; \
    } catch (const exc_type& e) { \
      caught = true; \
      std::string msg = e.what(); \
      if (msg.find(substr) == std::string::npos) { \
        std::cerr << "Exception did not contain expected substring '" << (substr) \
                  << "'. Got: '" << msg << "' at " << __FILE__ << ":" << __LINE__ << "\n"; \
        std::exit(1); \
      } \
    } catch (const std::exception& e) { \
      std::cerr << "Unexpected exception type: " << e.what() << " at " << __FILE__ << ":" << __LINE__ << "\n"; \
      std::exit(1); \
    } catch (...) { \
      std::cerr << "Unknown exception type at " << __FILE__ << ":" << __LINE__ << "\n"; \
      std::exit(1); \
    } \
    if (!caught) { \
      std::cerr << "Expected exception " #exc_type " was not thrown at " << __FILE__ << ":" << __LINE__ << "\n"; \
      std::exit(1); \
    } \
  } while (false)

using namespace cppwrf::io;

void test_basic_parsing() {
  const std::string text = R"(
# comment            ; comment too
[domain]
nx = 512
dx = 100.0           # inline comment allowed after whitespace
eta_levels = 1.0, 0.98, 0.95
[run]
case = straka
use_cuda = true
)";

  Config cfg = Config::from_string(text);
  ASSERT_TRUE(cfg.has("domain.nx"));
  ASSERT_EQ(cfg.get_int("domain.nx"), 512LL);
  ASSERT_TRUE(cfg.has("domain.dx"));
  ASSERT_EQ(cfg.get_double("domain.dx"), 100.0);
  ASSERT_TRUE(cfg.has("run.case"));
  ASSERT_EQ(cfg.get_string("run.case"), "straka");
  ASSERT_TRUE(cfg.has("run.use_cuda"));
  ASSERT_TRUE(cfg.get_bool("run.use_cuda"));

  auto eta = cfg.get_double_list("domain.eta_levels");
  ASSERT_EQ(eta.size(), 3UL);
  ASSERT_EQ(eta[0], 1.0);
  ASSERT_EQ(eta[1], 0.98);
  ASSERT_EQ(eta[2], 0.95);

  auto keys = cfg.keys();
  ASSERT_EQ(keys.size(), 5UL);
  ASSERT_EQ(keys[0], "domain.dx");
  ASSERT_EQ(keys[1], "domain.eta_levels");
  ASSERT_EQ(keys[2], "domain.nx");
  ASSERT_EQ(keys[3], "run.case");
  ASSERT_EQ(keys[4], "run.use_cuda");
}

void test_comments_and_whitespace() {
  const std::string text = R"(
   # Leading spaces comment
   ; Semicolon comment

   [sec1]   # section comment
   key1   =   val1   ; inline semicolon
   key2 = val2
)";
  Config cfg = Config::from_string(text);
  ASSERT_EQ(cfg.get_string("sec1.key1"), "val1");
  ASSERT_EQ(cfg.get_string("sec1.key2"), "val2");
}

void test_syntax_errors_with_line_numbers() {
  // Key before section header at line 2
  const std::string text_no_sec = "\nkey = val\n[sec]\n";
  ASSERT_THROWS_CONTAINS(Config::from_string(text_no_sec), ConfigError, ":2:");

  // Duplicate key at line 4
  const std::string text_dup = "[sec]\nk1 = 1\nk2 = 2\nk1 = 3\n";
  ASSERT_THROWS_CONTAINS(Config::from_string(text_dup), ConfigError, ":4:");
  ASSERT_THROWS_CONTAINS(Config::from_string(text_dup), ConfigError, "sec.k1");

  // Missing '=' delimiter at line 3
  const std::string text_bad_syntax = "[sec]\nk1 = 1\nbad_line_here\n";
  ASSERT_THROWS_CONTAINS(Config::from_string(text_bad_syntax), ConfigError, ":3:");

  // Unclosed bracket at line 2
  const std::string text_unclosed = "\n[unclosed\nk = 1\n";
  ASSERT_THROWS_CONTAINS(Config::from_string(text_unclosed), ConfigError, ":2:");

  // Invalid section name (uppercase/hyphen) at line 1
  ASSERT_THROWS_CONTAINS(Config::from_string("[Sec-1]\nk=1\n"), ConfigError, ":1:");
  // Invalid key name at line 2
  ASSERT_THROWS_CONTAINS(Config::from_string("[sec]\nBad-Key=1\n"), ConfigError, ":2:");
}

void test_overrides() {
  const std::string text = R"(
[domain]
nx = 512
[run]
case = straka
)";
  Config cfg = Config::from_string(text, {"domain.nx=1024", "run.case=warm_bubble", "physics.enable=yes"});
  ASSERT_EQ(cfg.get_int("domain.nx"), 1024LL);
  ASSERT_EQ(cfg.get_string("run.case"), "warm_bubble");
  ASSERT_TRUE(cfg.has("physics.enable"));
  ASSERT_TRUE(cfg.get_bool("physics.enable"));

  // Bad override syntax
  ASSERT_THROWS_CONTAINS(Config::from_string(text, {"no_equal_sign"}), ConfigError, "missing '='");
  ASSERT_THROWS_CONTAINS(Config::from_string(text, {"no_section=1"}), ConfigError, "must be section.key");
  ASSERT_THROWS_CONTAINS(Config::from_string(text, {"Bad-Sec.k=1"}), ConfigError, "invalid characters");
}

void test_getters_and_conversions() {
  const std::string text = R"(
[test]
str_val = hello world
int_val = -42
int_pos = +100
int_bad = 12abc
dbl_val = 3.14159
dbl_sci = -2.5e-3
dbl_bad = 12.3.4
bool_t1 = true
bool_t2 = 1
bool_t3 = yes
bool_t4 = YES
bool_f1 = false
bool_f2 = 0
bool_f3 = no
bool_f4 = NO
bool_bad = maybe
list_val = 1.0, 2.5, -3.0
list_bad = 1.0, bad, 3.0
list_empty_item = 1.0, , 3.0
empty_val =
)";

  Config cfg = Config::from_string(text);

  // String
  ASSERT_EQ(cfg.get_string("test.str_val"), "hello world");
  ASSERT_EQ(cfg.get_string("missing.key", "default"), "default");
  ASSERT_THROWS_CONTAINS(cfg.get_string("missing.key"), ConfigError, "missing.key");

  // Int
  ASSERT_EQ(cfg.get_int("test.int_val"), -42LL);
  ASSERT_EQ(cfg.get_int("test.int_pos"), 100LL);
  ASSERT_EQ(cfg.get_int("missing.key", 99LL), 99LL);
  ASSERT_THROWS_CONTAINS(cfg.get_int("test.int_bad"), ConfigError, "test.int_bad");
  ASSERT_THROWS_CONTAINS(cfg.get_int("test.dbl_val"), ConfigError, "test.dbl_val");
  ASSERT_THROWS_CONTAINS(cfg.get_int("test.empty_val"), ConfigError, "test.empty_val");
  ASSERT_THROWS_CONTAINS(cfg.get_int("missing.key"), ConfigError, "missing.key");

  // Double
  ASSERT_EQ(cfg.get_double("test.dbl_val"), 3.14159);
  ASSERT_EQ(cfg.get_double("test.dbl_sci"), -2.5e-3);
  ASSERT_EQ(cfg.get_double("missing.key", 1.23), 1.23);
  ASSERT_THROWS_CONTAINS(cfg.get_double("test.dbl_bad"), ConfigError, "test.dbl_bad");
  ASSERT_THROWS_CONTAINS(cfg.get_double("test.str_val"), ConfigError, "test.str_val");
  ASSERT_THROWS_CONTAINS(cfg.get_double("test.empty_val"), ConfigError, "test.empty_val");
  ASSERT_THROWS_CONTAINS(cfg.get_double("missing.key"), ConfigError, "missing.key");

  // Bool
  ASSERT_TRUE(cfg.get_bool("test.bool_t1"));
  ASSERT_TRUE(cfg.get_bool("test.bool_t2"));
  ASSERT_TRUE(cfg.get_bool("test.bool_t3"));
  ASSERT_TRUE(cfg.get_bool("test.bool_t4"));
  ASSERT_FALSE(cfg.get_bool("test.bool_f1"));
  ASSERT_FALSE(cfg.get_bool("test.bool_f2"));
  ASSERT_FALSE(cfg.get_bool("test.bool_f3"));
  ASSERT_FALSE(cfg.get_bool("test.bool_f4"));
  ASSERT_TRUE(cfg.get_bool("missing.key", true));
  ASSERT_FALSE(cfg.get_bool("missing.key", false));
  ASSERT_THROWS_CONTAINS(cfg.get_bool("test.bool_bad"), ConfigError, "test.bool_bad");
  ASSERT_THROWS_CONTAINS(cfg.get_bool("missing.key"), ConfigError, "missing.key");

  // Lists
  auto list = cfg.get_double_list("test.list_val");
  ASSERT_EQ(list.size(), 3UL);
  ASSERT_EQ(list[0], 1.0);
  ASSERT_EQ(list[1], 2.5);
  ASSERT_EQ(list[2], -3.0);
  ASSERT_THROWS_CONTAINS(cfg.get_double_list("test.list_bad"), ConfigError, "test.list_bad");
  ASSERT_THROWS_CONTAINS(cfg.get_double_list("test.list_empty_item"), ConfigError, "test.list_empty_item");
  ASSERT_THROWS_CONTAINS(cfg.get_double_list("missing.key"), ConfigError, "missing.key");
  auto empty_list = cfg.get_double_list("test.empty_val");
  ASSERT_EQ(empty_list.size(), 0UL);
}

void test_to_json() {
  const std::string text = R"(
[domain]
nx = 512
dx = 100.0
[run]
case = straka
)";
  Config cfg = Config::from_string(text);
  std::string json_str = cfg.to_json();
  const std::string expected = "{\"domain\":{\"dx\":\"100.0\",\"nx\":\"512\"},\"run\":{\"case\":\"straka\"}}";
  ASSERT_EQ(json_str, expected);
}

int main() {
  std::cout << "Running test_config...\n";
  test_basic_parsing();
  test_comments_and_whitespace();
  test_syntax_errors_with_line_numbers();
  test_overrides();
  test_getters_and_conversions();
  test_to_json();
  std::cout << "test_config: ALL TESTS PASSED\n";
  return 0;
}
