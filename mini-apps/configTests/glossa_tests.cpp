/*
 * Unit tests for the glossa config-expression parser (glossa.hpp).
 *
 * Run directly: ./glossa_tests
 * Exit code 0 = all passed, 1 = at least one failure.
 *
 * These tests exercise glossa::evaluate_config() in-process. For the full
 * config-file -> Readparameters pipeline see config_reader_test.cpp.
 */
#include "glossa.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_set>

static int g_pass = 0;
static int g_fail = 0;

static bool contains(const std::string& haystack, const std::string& needle) {
   return haystack.find(needle) != std::string::npos;
}

static void check_bool(const char* test, bool cond) {
   if (cond) {
      g_pass++;
      std::printf("  PASS %s\n", test);
   } else {
      g_fail++;
      std::printf("  FAIL %s\n", test);
   }
}

// Evaluate src; true if the output contains expectedLine. A thrown exception is a failure.
static bool expect_line(const char* test, const std::string& src, const std::string& expectedLine) {
   std::string out;
   try {
      out = glossa::evaluate_config(src, {});
   } catch (const std::exception& e) {
      g_fail++;
      std::printf("  FAIL %s (threw: %s)\n", test, e.what());
      return false;
   }
   if (!contains(out, expectedLine)) {
      g_fail++;
      std::printf("  FAIL %s\n    wanted output containing: [%s]\n    got:\n%s\n", test, expectedLine.c_str(), out.c_str());
      return false;
   }
   g_pass++;
   std::printf("  PASS %s\n", test);
   return true;
}

// Evaluate src; true if the output does NOT contain forbidden.
static bool expect_no_line(const char* test, const std::string& src, const std::string& forbidden) {
   std::string out = glossa::evaluate_config(src, {});
   if (contains(out, forbidden)) {
      g_fail++;
      std::printf("  FAIL %s\n    did not want output containing: [%s]\n    got:\n%s\n", test, forbidden.c_str(), out.c_str());
      return false;
   }
   g_pass++;
   std::printf("  PASS %s\n", test);
   return true;
}

int main() {
   std::printf("glossa_tests: expression parser coverage\n");

   std::printf("[numeric literals]\n");
   expect_line("integer literal", "[s]\nx = 42\n", "x = 42");
   expect_line("float with leading dot (.5)", "[s]\nx = .5 * 10\n", "x = 5");
   expect_line("plain float", "[s]\nx = 3.75\n", "x = 3.75");
   expect_line("scientific notation 40e3 (the PR example)", "[s]\nx = 40e3\n", "x = 40000");
   expect_line("scientific with explicit plus exponent", "[s]\nx = 1e+2\n", "x = 100");
   expect_line("scientific with negative exponent", "[s]\nx = 1.5e-3\n", "x = 0.0015");
   expect_line("trailing-dot number", "[s]\nx = 7.\n", "x = 7");
   expect_line("negative literal via unary minus", "[s]\nx = -5\n", "x = -5");

   std::printf("[unary minus / plus and precedence]\n");
   expect_line("PR example verbatim: -.5 * vx_length * 4 * velres",
               "[s]\nvx_length = 26\nvelres = 40e3\nvx_min = -.5 * vx_length * 4 * velres\n",
               "vx_min = -2080000");
   expect_line("double unary minus", "[s]\nx = --5\n", "x = 5");
   expect_line("unary plus ignored", "[s]\nx = +5\n", "x = 5");
   expect_line("minus binds tighter than plus: -1+2 = 1", "[s]\nx = -1 + 2\n", "x = 1");
   expect_line("mul before add: 1+2*3 = 7", "[s]\nx = 1 + 2 * 3\n", "x = 7");
   expect_line("parentheses override precedence: (1+2)*3 = 9", "[s]\nx = (1 + 2) * 3\n", "x = 9");
   expect_line("left-assoc division: 8/4/2 = 1", "[s]\nx = 8 / 4 / 2\n", "x = 1");
   expect_line("grouped division: 8/(4/2) = 4", "[s]\nx = 8 / (4 / 2)\n", "x = 4");
   expect_line("subtraction left-associative: 10-3-2 = 5", "[s]\nx = 10 - 3 - 2\n", "x = 5");
   expect_line("unary minus inside parens: (-2)*(-3) = 6", "[s]\nx = (-2) * (-3)\n", "x = 6");

   std::printf("[variables and local/global semantics]\n");
   expect_line("variable reuse within a section", "[s]\na = 10\nb = a + 5\n", "b = 15");
   expect_line("self-update uses the OLD value (a=a+1 then z=2a)", "[s]\na = 1\na = a + 1\nz = a * 2\n", "z = 4");
   expect_line("local helper does not leak into later sections' math",
               "[s1]\nlocal tmp = 7\nkeep = tmp\n[s2]\nother = keep + 1\n",
               "other = 8");
   {
      // After the section that declared it, the local name is gone; use of it fails cleanly.
      std::string out = glossa::evaluate_config("[s1]\nlocal tmp = 7\n[s2]\ny = tmp * 2\n", {});
      check_bool("local var unusable after its section (line passes through verbatim)", contains(out, "y = tmp * 2"));
   }
   expect_no_line("local variable is suppressed from emitted option lines (PR pattern)",
                  "[s]\nlocal velres = 40e3\nvx_min = -.5 * velres\n", "velres");
   expect_line("global variable survives a section change",
               "[s1]\nglobal g = 21\nx = g * 2\n[s2]\ny = g + 1\n", "y = 22");
   expect_no_line("global helper is not emitted as an option line either",
                  "[s1]\nglobal g = 21\nx = g * 2\n", "g = ");
   expect_line("cross-section reference with dotted name (section.key)",
               "[sec_a]\nx = 30\n[sec_b]\ny = sec_a.x / 6\n", "y = 5");

   std::printf("[functions]\n");
   expect_line("sqrt", "[s]\nx = sqrt(144)\n", "x = 12");
   expect_line("pow with expression args", "[s]\nx = pow(2, 10) + 1\n", "x = 1025");
   expect_line("min/max nested", "[s]\nx = min(max(3, 9), 7)\n", "x = 7");
   expect_line("clamp to upper bound", "[s]\nx = clamp(15, 0, 10)\n", "x = 10");
   expect_line("abs of a negative expression", "[s]\nx = abs(3 - 8)\n", "x = 5");
   expect_line("int() truncation", "[s]\nx = int(7.9)\n", "x = 7");
   expect_line("floor", "[s]\nf = floor(2.1)\n", "f = 2");
   expect_line("ceil", "[s]\nc = ceil(2.1)\n", "c = 3");
   expect_line("sin(0) = 0", "[s]\nx = sin(0)\n", "x = 0");
   expect_line("cos(0) = 1", "[s]\nx = cos(0)\n", "x = 1");
   expect_line("atan2 scaled to pi", "[s]\nx = atan2(1, 1) * 4\n", "x = 3.14159265358979");
   expect_line("fract fractional part", "[s]\nx = fract(2.75)\n", "x = 0.75");

   std::printf("[predefined physical constants]\n");
   expect_line("EPS_0 available as a global", "[s]\nx = EPS_0\n", "x = 8.85418782e-12");
   expect_line("MU_0*EPS_0 = 1/c^2 sanity value", "[s]\nx = MU_0 * EPS_0\n", "x = 1.11265005508126e-17");
   expect_line("CHARGE / MASS_PROTON ~ 9.58e7 C/kg", "[s]\nz = CHARGE / MASS_PROTON\n", "z = 95788345.0242224");
   {
      std::string out = glossa::evaluate_config("[s]\nEPS_0 = 12345\nx = EPS_0\n", {});
      check_bool("redefining a reserved constant is refused; builtin keeps winning", contains(out, "x = 8.85418782e-12"));
   }

   std::printf("[strings]\n");
   expect_line("string literal round-trips unquoted into the value slot", "[s]\nsay = \"hello world\"\n", "say = hello world");
   expect_line("string concatenation with +", "[s]\na = \"foo\"\nb = \"bar\"\nc = a + b\n", "c = foobar");
   {
      // Numeric + string is not supported: line must pass through verbatim (no silent garbage).
      std::string out = glossa::evaluate_config("[s]\nn = 42\nmsg = \"the answer is \" + n\n", {});
      check_bool("mixed string+number concat rejected, passes through verbatim", contains(out, "msg = \"the answer is \" + n"));
   }
   {
      // Escaped quotes inside a string survive lexing.
      std::string out = glossa::evaluate_config("[s]\nsay = \"a\\\"quoted\\\" word\"\n", {});
      check_bool("escaped quotes inside strings", contains(out, "say = a\"quoted\" word"));
   }

   std::printf("[comments]\n");
   {
      // The PR body itself uses a # comment right after a local: it must not break evaluation.
      std::string out = glossa::evaluate_config("[s]\nlocal velres = 40e3 # 40 km/s resolution\nvx_min = -.5 * velres\n", {});
      check_bool("# comment after value is stripped before evaluating", contains(out, "vx_min = -20000"));
   }
   expect_line("trailing comment on a plain arithmetic line", "[s]\nx = 1 + 2 # plus three\ny = x * 2\n", "y = 6");
   {
      std::string out = glossa::evaluate_config("[s]\nsay = \"hash # is not a comment\"\n", {});
      check_bool("# inside a quoted string does not start a comment", contains(out, "say = hash # is not a comment"));
   }
   expect_line("full-line comment passes through untouched", "[s]\n# this is a comment\nx = 1\n", "# this is a comment");

   std::printf("[pass-through / backwards compatibility]\n");
   expect_line("plain integer value unchanged (backwards compat)", "[s]\nx_length = 26\n", "x_length = 26");
   {
      std::string out = glossa::evaluate_config("[s]\nname = this is a plain string\n", {});
      check_bool("plain string value passes through unchanged", contains(out, "name = this is a plain string"));
   }
   {
      std::string out = glossa::evaluate_config("[s]\nfreeform line without equals\nx = 1\n", {});
      check_bool("line without '=' passes through verbatim", contains(out, "freeform line without equals"));
   }
   {
      // Trailing operator used to produce a null-child Expr that segfaulted at eval time.
      std::string out = glossa::evaluate_config("[s]\nx = 1 +\n", {});
      check_bool("incomplete expression (trailing '+') fails cleanly, no crash", contains(out, "x = 1 +"));
   }
   {
      std::string out = glossa::evaluate_config("[s]\nx = ((1 + 2\n", {});
      check_bool("unbalanced parenthesis fails cleanly, no crash", contains(out, "x = ((1 + 2"));
   }

   std::printf("[error handling: clean failures, no silent corruption]\n");
   {
      // Forward reference: variable used before it is defined.
      std::string out = glossa::evaluate_config("[s]\nz = q * 2\nq = 10\n", {});
      check_bool("forward reference fails cleanly (line passes through, later def still works)",
                 contains(out, "z = q * 2") && contains(out, "q = 10"));
   }
   {
      std::string out = glossa::evaluate_config("[s]\nz = frobnicate(3)\n", {});
      check_bool("unknown function fails cleanly and passes through", contains(out, "z = frobnicate(3)"));
   }
   {
      std::string out = glossa::evaluate_config("[s]\nz = 1 % 2\n", {});
      check_bool("invalid operator char fails cleanly and passes through", contains(out, "z = 1 % 2"));
   }
   {
      // The old (pre-fix) code turned this exact pattern into a stack overflow: with live
      // expression trees, a <-> b aliasing plus re-evaluation recursed forever. With the
      // value-snapshot commit semantics each line is evaluated against previously fixed
      // values, so mutual updates resolve to OLD values and later uses are safe: no hang,
      // no crash, fully deterministic.
      std::string out = glossa::evaluate_config("[s]\na = 1\nb = 2\na = b + 0.5\nb = a + 0.5\nz = a * 4\n", {});
      check_bool("mutual reassignment (former stack-overflow pattern) resolves on old values, no hang",
                 contains(out, "z = 10"));
   }
   {
      // Self-cycle through an alias: b = a then a = b + 0.5; later use of a is still safe.
      std::string out = glossa::evaluate_config("[s]\na = 1\nb = a\na = b + 0.5\nz = a * 2\n", {});
      check_bool("self-referential update chain stays finite (old values win)", contains(out, "z = 3"));
   }
   {
      // Direct exercise of the in-flight cycle detector at eval() level: hand-crafted Vars
      // whose expression trees genuinely reference each other must raise a clean error
      // instead of recursing forever. (evaluate_config itself only ever commits value
      // snapshots, so this guards the lower-level API.)
      const std::size_t N = 65536;
      void* mem = malloc(N);
      glossa::BumpAllocator arena(mem, N);
      glossa::Vars vars;
      vars["x"] = glossa::new_variable(arena, "y");
      vars["y"] = glossa::new_variable(arena, "x");
      bool threw = false;
      std::string msg;
      try {
         std::unordered_set<std::string> inflight;
         (void)glossa::eval(vars["x"], vars, &inflight);
      } catch (const std::exception& e) {
         threw = true;
         msg = e.what();
      }
      free(mem);
      check_bool("eval() raises a clean 'circular variable dependency' on adversarial var graphs",
                 threw && contains(msg, "circular"));
   }

   std::printf("\nglossa_tests: %d passed, %d failed\n", g_pass, g_fail);
   return g_fail == 0 ? 0 : 1;
}
