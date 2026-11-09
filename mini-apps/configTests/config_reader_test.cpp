/*
 * End-to-end tests for the config-file reading pipeline:
 *
 *   .cfg file  ->  glossa::evaluate_config (expression evaluation)
 *              ->  Readparameters::parse  (option registry, typed assignment)
 *
 * The binary is self-executing:
 *   - parent mode (default): writes temp .cfg files, re-runs itself as a child with
 *     --cftest-child <args...>, and asserts on exit code / stdout / stderr. Failure
 *     modes that call exit(1) or abort() inside Readparameters::parse() can only be
 *     observed from a separate process, hence the self-exec design.
 *   - child mode (--cftest-child): registers a fixed set of typed options in section
 *     "cftest", parses the given config file, and prints every option's resulting
 *     value as `OUT <name>=<value>` lines (UNSET if not touched).
 *
 * Run: mpirun -np 1 ./config_reader_test   (or plain ./config_reader_test)
 * Exit code 0 = all passed, 1 = at least one failure.
 */
#include "readparameters.h"

#include <mpi.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

// ---------------------------------------------------------------------------
// Options under test, registered in child mode.
// ---------------------------------------------------------------------------
static std::string cftest_str = "default-string";
static double cftest_dbl = 1.5;
static double cftest_dbl2 = -999.0;
static long long cftest_ll = -7;
static int cftest_int = 3;
static bool cftest_flag = false;
static std::vector<long long> cftest_vec;

static const char* OPTION_NAMES[] = {"cftest.str", "cftest.dbl", "cftest.dbl2", "cftest.ll",
                                     "cftest.int", "cftest.flag",  "cftest.vec"};

// ---------------------------------------------------------------------------
// Child mode
// ---------------------------------------------------------------------------
static int childMain(int argc, char** argv) {
   // Readparameters::make_tokens() skips argv[0] as the "program name"; synthesize that slot
   // because our dispatch consumed the --cftest-child marker before getting here.
   std::vector<char*> fake;
   static std::string prog = "config_reader_test";
   fake.push_back(prog.data());
   for (int i = 0; i < argc; ++i) {
      fake.push_back(argv[i]);
   }

   MPI_Init(&argc, &argv);
   const int cmdArgc = argc; // captured before any further mangling

   Readparameters params(cmdArgc + 1, fake.data());
   Readparameters::add<std::string>("cftest.str", "a string option", cftest_str, std::string("default-string"));
   Readparameters::add<double>("cftest.dbl", "a double option", cftest_dbl, 1.5);
   Readparameters::add<double>("cftest.dbl2", "second double option", cftest_dbl2, -999.0);
   Readparameters::add<long long>("cftest.ll", "a long long option", cftest_ll, -7LL);
   Readparameters::add<int>("cftest.int", "an int option", cftest_int, 3);
   Readparameters::addFlag("cftest.flag", "a flag option", cftest_flag);
   Readparameters::add<std::vector<long long>>("cftest.vec", "a vector option", cftest_vec);

   std::vector<std::string> invalid;
   std::vector<std::string> filenames;
   params.parse(invalid, filenames, false /* extras: be strict */, std::unordered_map<std::string, double>{});

   auto printVec = [](const std::vector<long long>& v) {
      std::string s;
      for (std::size_t i = 0; i < v.size(); ++i) {
         if (i) s += ",";
         s += std::to_string(v[i]);
      }
      return s;
   };

   printf("OUT cftest.str=[%s]\n", Readparameters::isSet("cftest.str") ? cftest_str.c_str() : "UNSET");
   auto f15 = [](double v) {
      std::ostringstream ss;
      ss << std::setprecision(15) << v; // same formatting glossa uses when writing values back out
      return ss.str();
   };
   printf("OUT cftest.dbl=%s\n", Readparameters::isSet("cftest.dbl") ? f15(cftest_dbl).c_str() : "UNSET");
   printf("OUT cftest.dbl2=%s\n", Readparameters::isSet("cftest.dbl2") ? f15(cftest_dbl2).c_str() : "UNSET");
   printf("OUT cftest.ll=%s\n", Readparameters::isSet("cftest.ll") ? std::to_string(cftest_ll).c_str() : "UNSET");
   printf("OUT cftest.int=%s\n", Readparameters::isSet("cftest.int") ? std::to_string(cftest_int).c_str() : "UNSET");
   printf("OUT cftest.flag=%d\n", cftest_flag ? 1 : 0);
   printf("OUT cftest.vec=[%s]\n", Readparameters::isSet("cftest.vec") ? printVec(cftest_vec).c_str() : "UNSET");
   fflush(stdout);

   MPI_Finalize();
   return 0;
}

// ---------------------------------------------------------------------------
// Test framework + self-exec harness (parent mode)
// ---------------------------------------------------------------------------
static int g_pass = 0;
static int g_fail = 0;

static void report(const char* test, bool cond, const std::string& detail = "") {
   if (cond) {
      g_pass++;
      printf("  PASS %s\n", test);
   } else {
      g_fail++;
      printf("  FAIL %s%s%s\n", test, detail.empty() ? "" : " :: ", detail.c_str());
   }
}

static bool contains(const std::string& hay, const std::string& needle) { return hay.find(needle) != std::string::npos; }

// Extract the `OUT <name>=<value>` payload for one option from child stdout.
static std::string outValue(const std::string& stdoutStr, const std::string& name) {
   const std::string prefix = "OUT " + name + "=";
   auto pos = stdoutStr.find(prefix);
   if (pos == std::string::npos) return "<NO SUCH LINE>";
   pos += prefix.size();
   auto end = stdoutStr.find('\n', pos);
   return stdoutStr.substr(pos, end == std::string::npos ? std::string::npos : end - pos);
}

struct ChildResult {
   int rc = -1;
   std::string out;
   std::string err;
};

// Run this binary as a child: --cftest-child <args...>. Captures stdout/stderr/rc.
static ChildResult runChild(const std::vector<std::string>& args, unsigned timeoutSec = 30) {
   ChildResult res;
   int outPipe[2], errPipe[2];
   if (pipe(outPipe) != 0 || pipe(errPipe) != 0) return res;

   pid_t pid = fork();
   if (pid < 0) return res;
   if (pid == 0) {
      // child: wire pipes, guard against hangs with alarm(), re-exec self.
      dup2(outPipe[1], STDOUT_FILENO);
      dup2(errPipe[1], STDERR_FILENO);
      close(outPipe[0]); close(outPipe[1]); close(errPipe[0]); close(errPipe[1]);
      alarm(timeoutSec);

      char exePath[4096];
      ssize_t n = readlink("/proc/self/exe", exePath, sizeof(exePath) - 1);
      if (n <= 0) _exit(127);
      exePath[n] = '\0';

      std::vector<char*> argv;
      argv.push_back(exePath);
      static std::string childFlag = "--cftest-child";
      argv.push_back(childFlag.data());
      for (auto& a : args) {
         // mutable copies: execv requires char*
         auto keep = new std::string(a);
         argv.push_back(keep->data());
      }
      argv.push_back(nullptr);
      execv(exePath, argv.data());
      _exit(127);
   }

   close(outPipe[1]); close(errPipe[1]);
   auto drain = [](int fd, std::string& dst) {
      char buf[4096];
      ssize_t r;
      while ((r = read(fd, buf, sizeof(buf))) > 0) dst.append(buf, r);
      close(fd);
   };
   drain(outPipe[0], res.out);
   drain(errPipe[0], res.err);

   int status = 0;
   waitpid(pid, &status, 0);
   if (WIFEXITED(status)) res.rc = WEXITSTATUS(status);
   else if (WIFSIGNALED(status)) res.rc = 128 + WTERMSIG(status);
   return res;
}

// Write a cfg file into the given directory; returns its full path.
static std::string writeCfg(const std::string& dir, const std::string& fname, const std::string& content) {
   std::string path = dir + "/" + fname;
   std::ofstream f(path);
   f << content;
   f.close();
   return path;
}

int parentMain(int argc, char** argv) {
   (void)argc; (void)argv;
   printf("config_reader_test: end-to-end config pipeline (glossa + Readparameters)\n");

   char dirTpl[] = "/tmp/cftestXXXXXX";
   char* dir = mkdtemp(dirTpl);
   if (!dir) { perror("mkdtemp"); return 1; }

   // ------------------------------------------------------------------
   printf("[basic typed values]\n");
   {
      std::string cfg = writeCfg(dir, "basic.cfg",
         "[cftest]\n"
         "str = hello world\n"
         "dbl = 2.5\n"
         "ll = 9007199254740993\n"   // beyond int64-safe range of a 32-bit int
         "int = 42\n"
         "vec = [1, 2, 3]\n");
      ChildResult r = runChild({"--run_config=" + cfg});
      std::string why;
      if (r.rc != 0) why = "rc=" + std::to_string(r.rc) + " err: " + r.err;
      report("plain values parse with rc=0", r.rc == 0, why);
      report("string value round-trips (note: config lines are space-stripped, so no internal spaces)",
             outValue(r.out, "cftest.str") == "[helloworld]", outValue(r.out, "cftest.str"));
      report("double value assigned", outValue(r.out, "cftest.dbl") == "2.5", outValue(r.out, "cftest.dbl"));
      // Long long options keep their full magnitude (2^53 is exactly representable; plain
      // decimal emission, no scientific notation that stoll would truncate).
      report("long long keeps 2^53 magnitude through the double round trip",
             outValue(r.out, "cftest.ll") == "9007199254740992", outValue(r.out, "cftest.ll"));
      report("int value assigned", outValue(r.out, "cftest.int") == "42", outValue(r.out, "cftest.int"));
      report("vector [1,2,3] parsed elementwise", outValue(r.out, "cftest.vec") == "[1,2,3]", outValue(r.out, "cftest.vec"));
      report("flag defaults to false when absent", outValue(r.out, "cftest.flag") == "0", outValue(r.out, "cftest.flag"));
   }

   // ------------------------------------------------------------------
   printf("[expressions via glossa]\n");
   {
      // The PR's motivating example, verbatim (scaled to our test section).
      std::string cfg = writeCfg(dir, "expr.cfg",
         "[cftest]\n"
         "local velres = 40e3 # 40 km/s resolution\n"
         "dbl = -.5 * 26 * 4 * velres\n"
         "dbl2 = (1 + 2) * 3 / 2\n"
         "ll = pow(2, 32) - 1\n"
         "int = sqrt(900)\n"
         "str = \"a quoted literal\"\n");
      ChildResult r = runChild({"--run_config=" + cfg});
      std::string why;
      if (r.rc != 0) why = "rc=" + std::to_string(r.rc) + " err: " + r.err;
      report("PR-style expression config parses with rc=0", r.rc == 0, why);
      report("-.5 * 26 * 4 * velres -> -2080000", outValue(r.out, "cftest.dbl") == "-2080000", outValue(r.out, "cftest.dbl"));
      report("parenthesised arithmetic (1+2)*3/2 -> 4.5", outValue(r.out, "cftest.dbl2") == "4.5", outValue(r.out, "cftest.dbl2"));
      report("pow(2,32)-1 fits long long exactly", outValue(r.out, "cftest.ll") == "4294967295", outValue(r.out, "cftest.ll"));
      report("sqrt(900) -> 30 for int option", outValue(r.out, "cftest.int") == "30", outValue(r.out, "cftest.int"));
      report("quoted string literal assigned to string option (space-stripped, quotes consumed)",
             outValue(r.out, "cftest.str") == "[aquotedliteral]", outValue(r.out, "cftest.str"));
      // The local helper must not have been attempted as an option.
      report("local helper does not reach the registry (no 'velres' complaint)", !contains(r.err, "velres"), r.err);
   }

   // ------------------------------------------------------------------
   printf("[predefined constants]\n");
   {
      std::string cfg = writeCfg(dir, "const.cfg",
         "[cftest]\n"
         "dbl = CHARGE / MASS_PROTON\n");
      ChildResult r = runChild({"--run_config=" + cfg});
      report("physical constant expression evaluates", outValue(r.out, "cftest.dbl") == "95788345.0242224", outValue(r.out, "cftest.dbl"));
   }

   // ------------------------------------------------------------------
   printf("[misspelled / unknown names must fail]\n");
   {
      std::string cfg = writeCfg(dir, "misspell.cfg",
         "[cftest]\n"
         "strb = hello\n");           // 'strb' does not exist (should be cftest.str)
      ChildResult r = runChild({"--run_config=" + cfg});
      report("unknown option in config -> nonzero exit", r.rc == 1, "rc=" + std::to_string(r.rc));
      report("error message lists the offending name", contains(r.err, "cftest.strb"), r.err);
      report("error mentions 'invalid'", contains(r.err, "invalid"), r.err);
   }
   {
      std::string cfg = writeCfg(dir, "misspell2.cfg",
         "[cftes]\n"                   // section name itself misspelled
         "dbl = 1\n");
      ChildResult r = runChild({"--run_config=" + cfg});
      report("misspelled SECTION is reported as an invalid option", r.rc == 1 && contains(r.err, "cftes.dbl"),
             "rc=" + std::to_string(r.rc) + " err: " + r.err);
   }
   {
      std::string cfg = writeCfg(dir, "missingval.cfg",
         "[cftest]\nstr = ok\n"
         "int");                       // option name without a value at all
      ChildResult r = runChild({"--run_config=" + cfg});
      report("option line with no '=' is skipped, not fatal", r.rc == 0 && outValue(r.out, "cftest.str") == "[ok]",
             "rc=" + std::to_string(r.rc) + " err: " + r.err);
   }

   // ------------------------------------------------------------------
   printf("[expression errors surface cleanly]\n");
   {
      std::string cfg = writeCfg(dir, "undefvar.cfg",
         "[cftest]\ndbl = nonexistent_var * 2\n");
      ChildResult r = runChild({"--run_config=" + cfg});
      report("undefined variable in expression -> nonzero exit (no silent garbage)", r.rc != 0, "rc=" + std::to_string(r.rc));
      report("glossa reported the undefined variable", contains(r.err, "nonexistent_var") && contains(r.err, "glossa"), r.err);
   }
   {
      // Former stack-overflow pattern: mutual reassignment of two options followed by reuse.
      // Under value-snapshot semantics every line reads the values fixed by earlier lines,
      // so this must resolve cleanly and deterministically (dbl=1.5, dbl2=2, ll=6) within
      // the timeout — no hang, no segfault.
      std::string cfg = writeCfg(dir, "cycle.cfg",
         "[cftest]\ndbl = 1\ndbl2 = dbl\ndbl = dbl2 + 0.5\ndbl2 = dbl + 0.5\nll = dbl * 4\n");
      ChildResult r = runChild({"--run_config=" + cfg}, 30);
      report("mutual reassignment (former segfault pattern) resolves on old values, no hang",
             r.rc == 0 && outValue(r.out, "cftest.ll") == "6" && outValue(r.out, "cftest.dbl2") == "2",
             "rc=" + std::to_string(r.rc) + " err: " + r.err);
   }
   {
      // Integer option fed with a non-integral double: stoll takes the leading integer
      // part. Document current semantics (silent truncation) rather than hide them.
      std::string cfg = writeCfg(dir, "trunc.cfg",
         "[cftest]\nint = 5.5\n");
      ChildResult r = runChild({"--run_config=" + cfg});
      report("documented quirk: int option <- 5.5 truncates to 5 (no error)", r.rc == 0 && outValue(r.out, "cftest.int") == "5",
             "rc=" + std::to_string(r.rc) + " val=" + outValue(r.out, "cftest.int"));
   }

   // ------------------------------------------------------------------
   printf("[command line overrides config]\n");
   {
      std::string cfg = writeCfg(dir, "override.cfg",
         "[cftest]\ndbl = 2.5\nint = 42\n");
      ChildResult r = runChild({"--run_config=" + cfg, "--cftest.dbl=9.75"});
      report("CLI value beats config-file value for same option", outValue(r.out, "cftest.dbl") == "9.75", outValue(r.out, "cftest.dbl"));
      report("options without CLI override keep the config value", outValue(r.out, "cftest.int") == "42", outValue(r.out, "cftest.int"));
   }

   // ------------------------------------------------------------------
   printf("[flags and bools]\n");
   {
      std::string cfg = writeCfg(dir, "flag.cfg",
         "[cftest]\nflag = true\n");
      ChildResult r = runChild({"--run_config=" + cfg});
      report("flag set from config ('true')", outValue(r.out, "cftest.flag") == "1", outValue(r.out, "cftest.flag"));

      ChildResult r2 = runChild({});   // no config file at all -> default "config.cfg" is missing; everything stays default
      report("missing config file does not crash (defaults survive)", r2.rc == 0, "rc=" + std::to_string(r2.rc));
      // parse() enters via resetAll(): every option is restored to its registered
      // default value AND marked wasSet=false. Options then set from config/CLI get
      // wasSet=true. So with nothing supplied, isSet() is false for everything even
      // though the underlying values hold the registered defaults.
      report("nothing supplied -> all options back at default and reported UNSET by isSet", outValue(r2.out, "cftest.str") == "[UNSET]" && outValue(r2.out, "cftest.dbl") == "UNSET" && outValue(r2.out, "cftest.int") == "UNSET",
             "str=" + outValue(r2.out, "cftest.str") + " dbl=" + outValue(r2.out, "cftest.dbl") + " int=" + outValue(r2.out, "cftest.int"));
   }

   printf("\nconfig_reader_test: %d passed, %d failed\n", g_pass, g_fail);
   return g_fail == 0 ? 0 : 1;
}

// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
   // Parent mode must not need a real MPI world (fork+exec children each init their own).
   for (int i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == "--cftest-child") {
         // Everything after the marker is the child's own argument vector.
         return childMain(argc - i - 1, argv + i + 1);
      }
   }
   MPI_Init(&argc, &argv); // keep the binary honest when run directly under mpirun
   int rc = parentMain(argc, argv);
   MPI_Finalize();
   return rc;
}
