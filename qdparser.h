/* File: qdargparser.h
A sub 100LOC quick and dirty parser.
Author: Kostis Papadakis, 2025 (kpapadakis@protonmail.com)
Usage:
  constexpr auto delim=' '; //for no delim use '0'
  QdArgParser<delim> cli_args(argc, argv);
  for (auto it=cli_args.begin(); it!=cli_args.end();++it){
    std::cout<<*it<<"\n";
  }
This program is free software; you can redistribute it and/or
modify it as you wish. You can make money with it or scam people.
*/
#include <cstring>
#include <iostream>
#include <iterator>
#include <memory>
#include <string_view>
#include <vector>
#include <algorithm>
template <char _delim, typename Allocator = std::allocator<char>>
struct QdArgParser {
  std::vector<char, Allocator> _data;
  QdArgParser(int argc, char **buffer) {
    if (argc > 0) {
      for (int i = 0; i < argc; ++i) {
        _data.insert(_data.end(), buffer[i],
                     buffer[i] + std::strlen(buffer[i]) + 1);
      }
      std::replace(_data.begin(), _data.end(), '\0', ' ');
      _data.pop_back();
      _data.shrink_to_fit();
    }
  }
  QdArgParser(const char *buffer, std::size_t size) {
    if (size > 0) {
      _data.insert(_data.end(), buffer, buffer + size + 1);
      std::replace(_data.begin(), _data.end(), '\0', ' ');
      _data.pop_back();
      _data.shrink_to_fit();
    }
  }
  template <typename F> auto apply(F &&functor) const noexcept {
    return functor(_data.begin(), _data.end());
  }
  template <typename T, typename F>
  std::vector<T> apply_collect_per_element(F &&functor) noexcept {
    std::vector<T> vals;
    for (auto it = begin(); it != end(); it = it.next()) {
      vals.push_back(functor(*it));
    }
    return vals;
  }

  std::size_t len() noexcept {
    std::size_t cnt = {0};
    for (const auto &token : *this) {
      (void)token;
      cnt++;
    }
    return cnt;
  }

  struct Iterator {
    using iterator_category = std::forward_iterator_tag;
    std::string_view view;
    char *_pos;
    std::size_t advance = 0;
    std::size_t left_over = 0;
    Iterator(char *pos, std::size_t len) : _pos(pos) {
      if (len == 0) {
        return;
      }
      if constexpr (_delim == '0') {
        view = std::string_view(pos, 1);
        advance = 1;
        left_over = len - 1;
        return;
      }
      char *delimiter = std::find(pos, pos + len, _delim);
      const std::size_t token_length = static_cast<std::size_t>(delimiter - pos);
      view = std::string_view(pos, token_length);
      if (delimiter == pos + len) {
        advance = token_length;
        left_over = 0;
      } else {
        advance = token_length + 1;
        left_over = len - advance;
      }
    }
    std::string_view operator*() { return view; }
    const char *operator->() { return view.data(); }
    bool operator==(const Iterator &other) const noexcept { return _pos == other._pos;}
    bool operator!=(const Iterator &other) const noexcept { return !(*this == other); }
    char *data() { return _pos; }
    Iterator next() { return Iterator(_pos + advance, left_over); }
    Iterator &operator++() {*this = next();return *this; }
  };
  Iterator begin() { return Iterator(_data.data(), _data.size()); }
  Iterator end() { return Iterator(_data.data() + _data.size(), 0); }
};
