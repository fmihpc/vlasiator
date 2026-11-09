// Little sleep, lots of coffee and Vlasiatoring at GH200 Hackathon @CSC
// Glossa means tongue in Greek :)
#pragma once
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static const std::string LOCALKW = "local";
static const std::string GLOBALKW = "global";
static std::unordered_map<std::string, double> predefined_globals = {
    {"EPS_0", 8.85418782e-12},
    {"MU_0", 1.25663706e-6},
    {"K_B", 1.3806503e-23},
    {"CHARGE", 1.60217653e-19},
    {"MASS_ELECTRON", 9.10938188e-31},
    {"MASS_PROTON", 1.67262158e-27},
    {"R_E", 6.3712e6},
};

static std::unordered_map<std::string, double> external_globals;

#define GLOSSA_FATAL(msg) throw std::runtime_error(msg)

namespace glossa {
   struct BumpAllocator {
      void* mem = nullptr;
      std::size_t sp = 0;
      std::size_t cap = 0;
      BumpAllocator(void* buf, std::size_t bytes) : mem(buf), sp(0), cap(bytes) {}
      BumpAllocator(const BumpAllocator&) = delete;
      BumpAllocator& operator=(const BumpAllocator&) = delete;
      template <class T> T* allocate(std::size_t n, int align_force = -1) {
         if (n == 0) {
            return nullptr;
         }
         std::size_t need = n * sizeof(T);
         std::size_t al = (align_force > 0) ? (std::size_t)align_force : std::max<std::size_t>(alignof(T), 8);
         std::size_t base = (std::size_t)((char*)mem + sp);
         std::size_t pad = (base % al == 0) ? 0 : (al - (base % al));
         if (sp + pad + need > cap) {
            GLOSSA_FATAL("OOM");
            return nullptr;
         }
         void* p = (char*)mem + sp + pad;
         sp += pad + need;
         return (T*)p;
      }
      // realloc esssentially?
      template <class T> void unsafe_extend_allocation(std::size_t extraElems) {
         std::size_t extra = extraElems * sizeof(T);
         if (sp + extra > cap) {
            GLOSSA_FATAL("OOM");
         };
         sp += extra;
      }
      void release() { sp = 0; }
   };

#undef GLOSSA_FATAL

   enum class BinaryOp { Add, Sub, Mul, Div, Invalid };

   struct Expr;

   struct Expr {
      enum class Kind { Num, String, Var, Bin, Callable } kind;
      double number = 0;
      std::string text;
      Expr *lhs = nullptr, *rhs = nullptr;
      BinaryOp op = BinaryOp::Invalid;
      std::vector<Expr*> args;
   };

   inline Expr* alloc_expr(BumpAllocator& arena) { return new (arena.allocate<Expr>(1)) Expr(); }

   inline Expr* new_number(BumpAllocator& arena, double x) {
      Expr* e = alloc_expr(arena);
      e->kind = Expr::Kind::Num;
      e->number = x;
      return e;
   }

   inline Expr* new_string(BumpAllocator& arena, std::string s) {
      Expr* e = alloc_expr(arena);
      e->kind = Expr::Kind::String;
      e->text = s;
      return e;
   }

   inline Expr* new_variable(BumpAllocator& arena, std::string name) {
      Expr* e = alloc_expr(arena);
      e->kind = Expr::Kind::Var;
      e->text = name;
      return e;
   }

   inline Expr* new_binary(BumpAllocator& arena, Expr* l, BinaryOp op, Expr* r) {
      Expr* e = alloc_expr(arena);
      e->kind = Expr::Kind::Bin;
      e->lhs = l;
      e->op = op;
      e->rhs = r;
      return e;
   }

   inline Expr* new_call(BumpAllocator& arena, std::string name, std::vector<Expr*> args) {
      Expr* e = alloc_expr(arena);
      e->kind = Expr::Kind::Callable;
      e->text = name;
      e->args = args;
      return e;
   }

   struct Value {
      enum class Kind { Number, String } kind;
      double number = 0;
      std::string text;

      double to_number() const {
         if (kind != Kind::Number) {
            throw std::runtime_error("ERROR: expected number");
         }
         return number;
      }
   };

   inline Value numeric(double x) { return Value{.kind = Value::Kind::Number, .number = x, .text = ""}; }
   inline Value str(std::string s) { return Value{.kind = Value::Kind::String, .number = 0, .text = s}; }

   inline std::string to_string(const Value& v) {
      if (v.kind == Value::Kind::Number) {
         std::ostringstream oss;
         oss << v.number;
         return oss.str();
      }
      return "\"" + v.text + "\"";
   }

   struct Token {
      enum class Kind { Number, String, Id, Plus, Minus, Star, Slash, Equal, LParen, RParen, Comma, _EOF } kind;
      double number = 0;
      std::string text;
   };

   inline std::vector<Token> lex(const std::string& input) {
      std::vector<Token> tokens;
      size_t i = 0, n = input.size();
      while (i < n) {
         char c = input[i];
         if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
            i++;
         } else if (std::isdigit((unsigned char)c) || c == '.') {
            std::string s;
            while (i < n && (std::isdigit((unsigned char)input[i]) || input[i] == '.'))
               s += input[i++];
            if (i < n && (input[i] == 'e' || input[i] == 'E')) {
               size_t save = i;
               std::string exp;
               exp += input[i++];
               if (i < n && (input[i] == '+' || input[i] == '-')) {
                  exp += input[i++];
               }
               if (i < n && std::isdigit((unsigned char)input[i])) {
                  while (i < n && std::isdigit((unsigned char)input[i]))
                     exp += input[i++];
                  s += exp;
               } else {
                  i = save;
               }
            }
            tokens.push_back({.kind = Token::Kind::Number, .number = std::stod(s), .text = ""});
         } else if (c == '"') {
            i++;
            std::string s;
            while (i < n && input[i] != '"') {
               if (input[i] == '\\' && i + 1 < n) {
                  char esc = input[i + 1];
                  i += 2;
                  switch (esc) {
                  case 'n':
                     s += '\n';
                     break;
                  case 't':
                     s += '\t';
                     break;
                  case 'r':
                     s += '\r';
                     break;
                  case '0':
                     s += '\0';
                     break;
                  default:
                     s += esc;
                     break;
                  }
                  continue;
               }
               s += input[i++];
            }
            if (i < n) {
               i++;
            }
            tokens.push_back({.kind = Token::Kind::String, .number = 0, .text = s});
         } else if (std::isalpha((unsigned char)c) || c == '_') {
            std::string s;
            while (i < n && (std::isalnum((unsigned char)input[i]) || input[i] == '_' || input[i] == '.'))
               s += input[i++];
            tokens.push_back({.kind = Token::Kind::Id, .number = 0, .text = s});
         } else {
            switch (c) {
            case '+':
               tokens.push_back({.kind = Token::Kind::Plus, .number = 0, .text = ""});
               break;
            case '-':
               tokens.push_back({.kind = Token::Kind::Minus, .number = 0, .text = ""});
               break;
            case '*':
               tokens.push_back({.kind = Token::Kind::Star, .number = 0, .text = ""});
               break;
            case '/':
               tokens.push_back({.kind = Token::Kind::Slash, .number = 0, .text = ""});
               break;
            case '=':
               tokens.push_back({.kind = Token::Kind::Equal, .number = 0, .text = ""});
               break;
            case '(':
               tokens.push_back({.kind = Token::Kind::LParen, .number = 0, .text = ""});
               break;
            case ')':
               tokens.push_back({.kind = Token::Kind::RParen, .number = 0, .text = ""});
               break;
            case ',':
               tokens.push_back({.kind = Token::Kind::Comma, .number = 0, .text = ""});
               break;
            default:
               throw std::runtime_error(std::string("ERROR: invalid character I do not know what to do!!!: ") + c);
            }
            i++;
         }
      }
      tokens.push_back({.kind = Token::Kind::_EOF, .number = 0, .text = ""});
      return tokens;
   }

   struct GlossaParser {
      std::vector<Token> tokens;
      size_t pos = 0;
      BumpAllocator& arena;

      GlossaParser(std::vector<Token> t, BumpAllocator& a) : tokens(t), arena(a) {}
      const Token& current() const { return tokens[pos]; }
      void advance() { pos++; }

      Expr* parse_expr() { return parse_add(); }
      Expr* parse_add() {
         auto expr = parse_mul();
         for (;;) {
            BinaryOp op;
            if (current().kind == Token::Kind::Plus) {
               op = BinaryOp::Add;
            } else if (current().kind == Token::Kind::Minus) {
               op = BinaryOp::Sub;
            } else {
               break;
            }
            advance();
            if (is_operand_terminator()) {
               throw std::runtime_error("ERROR: missing operand after " + describe_token(current()));
            }
            expr = new_binary(arena, expr, op, parse_mul());
         }
         return expr;
      }

      Expr* parse_mul() {
         auto expr = parse_primary();
         for (;;) {
            BinaryOp op;
            if (current().kind == Token::Kind::Star) {
               op = BinaryOp::Mul;
            } else if (current().kind == Token::Kind::Slash) {
               op = BinaryOp::Div;
            } else {
               break;
            }
            advance();
            if (is_operand_terminator()) {
               throw std::runtime_error("ERROR: missing operand after " + describe_token(current()));
            }
            expr = new_binary(arena, expr, op, parse_primary());
         }
         return expr;
      }

      // A binary operator must always be followed by an operand. Without this check a
      // trailing operator (e.g. "x = 1 +") produces an Expr node with a null child which
      // segfaults at evaluation time.
      bool is_operand_terminator() const {
         switch (current().kind) {
         case Token::Kind::_EOF:
            return true;
         case Token::Kind::RParen:
            return true;
         default:
            return false;
         }
      }

      std::string describe_token(const Token& t) const {
         switch (t.kind) {
         case Token::Kind::_EOF:
            return "end of line";
         case Token::Kind::LParen:
            return "'('";
         case Token::Kind::RParen:
            return "')'";
         case Token::Kind::Plus:
            return "'+'";
         case Token::Kind::Minus:
            return "'-'";
         case Token::Kind::Star:
            return "'*'";
         case Token::Kind::Slash:
            return "'/'";
         case Token::Kind::Equal:
            return "'='";
         case Token::Kind::Comma:
            return "','";
         case Token::Kind::Number:
            return std::to_string(t.number);
         case Token::Kind::Id:
            return t.text;
         case Token::Kind::String:
            return "\"...\"";
         }
         return "?";
      }

      Expr* parse_primary() {
         Token t = current();
         switch (t.kind) {
         case Token::Kind::Number:
            advance();
            return new_number(arena, t.number);
         case Token::Kind::String:
            advance();
            return new_string(arena, t.text);
         case Token::Kind::Id: {
            advance();
            if (current().kind == Token::Kind::LParen) {
               advance();
               std::vector<Expr*> args;
               if (current().kind != Token::Kind::RParen) {
                  for (;;) {
                     args.push_back(parse_expr());
                     if (current().kind == Token::Kind::Comma) {
                        advance();
                     } else {
                        break;
                     }
                  }
               }
               if (current().kind != Token::Kind::RParen) {
                  throw std::runtime_error("ERROR: missing parenthesis");
               }
               advance();
               return new_call(arena, t.text, args);
            }
            return new_variable(arena, t.text);
         }
         case Token::Kind::LParen: {
            advance();
            auto expr = parse_expr();
            if (current().kind != Token::Kind::RParen) {
               throw std::runtime_error("ERROR: missing parenthesis");
            }
            advance();
            return expr;
         }
         case Token::Kind::Minus:
            advance();
            return new_binary(arena, new_number(arena, 0.0), BinaryOp::Sub, parse_primary());
         case Token::Kind::Plus:
            advance();
            return parse_primary();
         default:
            throw std::runtime_error("ERROR: invlaid token I do not know what to do!");
         }
      }
   };

   using Vars = std::unordered_map<std::string, Expr*>;

   // Resolves an expression against vars. inflight tracks the variable names currently being
   // resolved (by name) so that a cyclic definition like a = b / b = a raises a clean
   // "circular dependency" error instead of recursing until the stack overflows.
   inline Value eval(const Expr* expr, const Vars& vars, std::unordered_set<std::string>* inflight = nullptr) {
      switch (expr->kind) {
      case Expr::Kind::Num:
         return numeric(expr->number);
      case Expr::Kind::String:
         return str(expr->text);
      case Expr::Kind::Var: {
         auto it = vars.find(expr->text);
         if (it == vars.end()) {
            throw std::runtime_error("ERROR: invalid var " + expr->text);
         }
         if (inflight != nullptr) {
            if (!inflight->insert(expr->text).second) {
               throw std::runtime_error("ERROR: circular variable dependency involving '" + expr->text + "'");
            }
            Value r = eval(it->second, vars, inflight);
            inflight->erase(expr->text);
            return r;
         }
         return eval(it->second, vars);
      }
      case Expr::Kind::Bin: {
         Value a = eval(expr->lhs, vars, inflight);
         Value b = eval(expr->rhs, vars, inflight);
         if (a.kind == Value::Kind::Number && b.kind == Value::Kind::Number) {
            switch (expr->op) {
            case BinaryOp::Add:
               return numeric(a.number + b.number);
            case BinaryOp::Sub:
               return numeric(a.number - b.number);
            case BinaryOp::Mul:
               return numeric(a.number * b.number);
            case BinaryOp::Div:
               return numeric(a.number / b.number);
            case BinaryOp::Invalid:
               abort();
            }
         }
         if (a.kind == Value::Kind::String && b.kind == Value::Kind::String && expr->op == BinaryOp::Add) {
            return str(a.text + b.text);
         }
         throw std::runtime_error("ERROR: invalid binary operation");
      }
      case Expr::Kind::Callable: {
         std::vector<Value> values;
         for (auto& a : expr->args) {
            values.push_back(eval(a, vars, inflight));
         }
         const std::string& name = expr->text;
         auto arg = [&](size_t i) { return values[i].to_number(); };
         if (name == "sqrt") {
            return numeric(std::sqrt(arg(0)));
         }
         if (name == "cbrt") {
            return numeric(std::cbrt(arg(0)));
         }
         if (name == "sin") {
            return numeric(std::sin(arg(0)));
         }
         if (name == "cos") {
            return numeric(std::cos(arg(0)));
         }
         if (name == "tan") {
            return numeric(std::tan(arg(0)));
         }
         if (name == "asin") {
            return numeric(std::asin(arg(0)));
         }
         if (name == "acos") {
            return numeric(std::acos(arg(0)));
         }
         if (name == "atan") {
            return numeric(std::atan(arg(0)));
         }
         if (name == "sinh") {
            return numeric(std::sinh(arg(0)));
         }
         if (name == "cosh") {
            return numeric(std::cosh(arg(0)));
         }
         if (name == "tanh") {
            return numeric(std::tanh(arg(0)));
         }
         if (name == "asinh") {
            return numeric(std::asinh(arg(0)));
         }
         if (name == "acosh") {
            return numeric(std::acosh(arg(0)));
         }
         if (name == "atanh") {
            return numeric(std::atanh(arg(0)));
         }
         if (name == "exp") {
            return numeric(std::exp(arg(0)));
         }
         if (name == "exp2") {
            return numeric(std::exp2(arg(0)));
         }
         if (name == "ln" || name == "log") {
            return numeric(std::log(arg(0)));
         }
         if (name == "log2") {
            return numeric(std::log2(arg(0)));
         }
         if (name == "log10") {
            return numeric(std::log10(arg(0)));
         }
         if (name == "abs") {
            return numeric(std::fabs(arg(0)));
         }
         if (name == "sign") {
            double x = arg(0);
            return numeric(x > 0 ? 1.0 : (x < 0 ? -1.0 : 0.0));
         }
         if (name == "floor") {
            return numeric(std::floor(arg(0)));
         }
         if (name == "ceil") {
            return numeric(std::ceil(arg(0)));
         }
         if (name == "round") {
            return numeric(std::round(arg(0)));
         }
         if (name == "fract") {
            double ip;
            return numeric(std::modf(arg(0), &ip));
         }
         if (name == "pow") {
            return numeric(std::pow(arg(0), arg(1)));
         }
         if (name == "min") {
            return numeric(std::min(arg(0), arg(1)));
         }
         if (name == "max") {
            return numeric(std::max(arg(0), arg(1)));
         }
         if (name == "atan2") {
            return numeric(std::atan2(arg(0), arg(1)));
         }
         if (name == "hypot") {
            return numeric(std::hypot(arg(0), arg(1)));
         }
         if (name == "clamp") {
            return numeric(std::min(std::max(arg(0), arg(1)), arg(2)));
         }
         if (name == "int") {
            return numeric(static_cast<int>(arg(0)));
         }
         if (name == "print") {
            for (size_t i = 0; i < values.size(); i++) {
               std::cout << (i ? " " : "");
               if (values[i].kind == Value::Kind::String) {
                  std::cout << values[i].text;
               } else {
                  std::cout << values[i].number;
               }
            }
            std::cout << std::endl;
            return values.empty() ? numeric(0) : values[0];
         }
         throw std::runtime_error("ERROR: unknown function: " + name);
      }
      }
      throw std::runtime_error("ERROR: unreachable");
   }

   inline Value execute_line(const std::string& source, Vars& vars, BumpAllocator& arena) {
      auto tokens = lex(source);
      if (tokens.size() >= 2 && tokens[0].kind == Token::Kind::Id && tokens[1].kind == Token::Kind::Equal) {
         const std::string name = tokens[0].text;
         std::vector<Token> rest(tokens.begin() + 2, tokens.end());
         if (rest.empty() || rest.back().kind != Token::Kind::_EOF) {
            rest.push_back({.kind = Token::Kind::_EOF, .number = 0, .text = ""});
         }
         GlossaParser parser(rest, arena);
         auto expr = parser.parse_expr();
         std::unordered_set<std::string> inflight;
         Value value = eval(expr, vars, &inflight);
         vars[name] = expr;
         return value;
      }
      GlossaParser parser(tokens, arena);
      auto expr = parser.parse_expr();
      return eval(expr, vars);
   }

   inline const std::unordered_map<std::string, double>& reserved_constants() {
      return predefined_globals;
   }

   inline size_t find_expr_no_quotes(const std::string& s, char needle) {
      bool in_string = false;
      for (size_t i = 0; i < s.size(); i++) {
         if (in_string && s[i] == '\\' && i + 1 < s.size()) {
            i++;
            continue;
         }
         if (s[i] == '"') {
            in_string = !in_string;
            continue;
         }
         if (!in_string && s[i] == needle) {
            return i;
         }
      }
      return std::string::npos;
   }

   inline std::string evaluate_config(const std::string& source, Vars& vars, BumpAllocator& arena) {
      auto trim = [](const std::string& s) {
         size_t b = s.find_first_not_of(" \t\r\n");
         if (b == std::string::npos) {
            return std::string("");
         }
         size_t e = s.find_last_not_of(" \t\r\n");
         return s.substr(b, e - b + 1);
      };

      for (const auto& [name, val] : reserved_constants()) {
         vars[name] = new_number(arena, val);
      }

      std::istringstream iss(source);
      std::ostringstream result;
      std::string section;
      std::string line;
      std::vector<std::string> scratchVariables;
      size_t lineno = 0;
      while (std::getline(iss, line)) {
         ++lineno;
         size_t first = line.find_first_not_of(" \t");
         if (first == std::string::npos) {
            result << line << "\n";
            continue;
         }

         char c = line[first];
         if (c == '[') {
            for (const auto& name : scratchVariables) {
               vars.erase(name);
               if (!section.empty()) {
                  vars.erase(section + "." + name);
               }
            }
            scratchVariables.clear();
            size_t end = line.find(']', first);
            if (end != std::string::npos) {
               section = line.substr(first + 1, end - first - 1);
            }
            result << line << "\n";
            continue;
         }
         if (c == '#') {
            result << line << "\n";
            continue;
         }

         size_t eq_pos = find_expr_no_quotes(line, '=');
         if (eq_pos == std::string::npos) {
            bool executed = false;
            try {
               auto tokens = lex(trim(line));
               GlossaParser parser(tokens, arena);
               Expr* expr = parser.parse_expr();
               if (expr->kind == Expr::Kind::Callable && parser.current().kind == Token::Kind::_EOF) {
                  eval(expr, vars);
                  executed = true;
               }
            } catch (const std::exception&) {
               executed = false;
            }
            if (!executed) {
               result << line << "\n";
            }
            continue;
         }
         std::string key = trim(line.substr(0, eq_pos));
         bool isLocal = false;
         bool isGlobal = false;
         if (key.size() > LOCALKW.size() && key.compare(0, LOCALKW.size(), LOCALKW) == 0 && std::isspace((unsigned char)key[LOCALKW.size()])) {
            isLocal = true;
            key = trim(key.substr(LOCALKW.size()));
         } else if (key.size() > GLOBALKW.size() && key.compare(0, GLOBALKW.size(), GLOBALKW) == 0 && std::isspace((unsigned char)key[GLOBALKW.size()])) {
            isGlobal = true;
            key = trim(key.substr(GLOBALKW.size()));
         }
         std::string after_eq = line.substr(eq_pos + 1);
         size_t hash_pos = find_expr_no_quotes(after_eq, '#');
         std::string value_raw = hash_pos == std::string::npos ? after_eq : after_eq.substr(0, hash_pos);
         std::string comment_suffix = hash_pos == std::string::npos ? "" : after_eq.substr(hash_pos);
         std::string expr_text = trim(value_raw);
         if (key.empty() || expr_text.empty()) {
            result << line << "\n";
            continue;
         }
         if (reserved_constants().count(key)) {
            std::cerr << "WARNING: '" << key << "' is a reserved global and you have just redefined it!!! " << line << std::endl;
            continue;
         }

         Expr* expr = nullptr;
         Value value;
         bool ok = false;
         std::string errmsg;
         try {
            auto tokens = lex(expr_text);
            GlossaParser parser(tokens, arena);
            expr = parser.parse_expr();
            if (parser.current().kind == Token::Kind::_EOF) {
               std::unordered_set<std::string> inflight;
               value = eval(expr, vars, &inflight);
               ok = true;
            } else {
               errmsg = "incomplete expression";
            }
         } catch (const std::exception& e) {
            errmsg = e.what();
         }
         if (!ok) {
            // Report the failure and pass the line through unchanged so that downstream
            // option handling sees the original (broken) value instead of us silently
            // converting a bad expression into an opaque string.
            std::cerr << "ERROR: glossa: cannot evaluate '" << expr_text << "' for variable '" << key
                      << "' at line " << lineno << ": " << errmsg << "\n";
            result << line << "\n";
            continue;
         }

         // Snapshot: store the evaluated VALUE, not the live expression tree. This keeps
         // later re-evaluations independent of the variables used to compute this one:
         // - `a = a + 1` must not commit a self-referential Expr graph (stack overflow on
         //   any later use of 'a');
         // - locals removed at section end must not dangle inside other expressions.
         Expr* stored = (value.kind == Value::Kind::Number) ? new_number(arena, value.number) : new_string(arena, value.text);
         vars[key] = stored;
         if (!section.empty()) {
            vars[section + "." + key] = stored;
         }
         if (isLocal) {
            scratchVariables.push_back(key);
            continue;
         }
         if (isGlobal) {
            continue;
         }
         if (!ok) {
            result << line << "\n";
            continue;
         }

         std::string out_line = key + " = ";
         if (value.kind == Value::Kind::Number) {
            std::ostringstream oss;
            // Integral values inside the long long range are emitted in plain decimal. The
            // defaultfloat formatting switches to scientific notation for large magnitudes,
            // which downstream integer option parsing (stoll) would silently truncate at
            // the decimal point (e.g. 9007199254740992 -> "9.00719925474099e+15" -> 9).
            constexpr double LL_MIN_SAFE = -9.223372036854775e18;
            constexpr double LL_MAX_SAFE = 9.223372036854775e18;
            if (std::isfinite(value.number) && value.number == std::floor(value.number) &&
                value.number >= LL_MIN_SAFE && value.number <= LL_MAX_SAFE) {
               oss << static_cast<long long>(value.number);
            } else {
               oss << std::setprecision(15) << value.number;
            }
            out_line += oss.str();
         } else {
            out_line += value.text;
         }
         if (!comment_suffix.empty()) {
            out_line += " " + comment_suffix;
         }
         result << out_line << "\n";
      }
      for (const auto& name : scratchVariables) {
         vars.erase(name);
         if (!section.empty()) {
            vars.erase(section + "." + name);
         }
      }
      return result.str();
   }

   inline std::string evaluate_config(const std::string& source, const std::unordered_map<std::string, double> &supplied_globals) {
      constexpr std::size_t N = 1024 * 1024;
      for (auto [key, val] : supplied_globals) {
         predefined_globals[key] = val;
      }
      void* mem = malloc(N);
      if (!mem) {
         throw std::runtime_error("Could not allocate memory for glossa's bump allocator!");
      }
      BumpAllocator arena(mem, N);
      Vars vars;
      std::string result = evaluate_config(source, vars, arena);
      arena.release();
      free(mem);
      return result;
   }
} // namespace glossa
