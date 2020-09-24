/* © 2014-2019, ETH Zurich, Institute of Fluid Dynamics, Stephan Küchlin */ 

/*
 * MyString.h
 *
 *  Created on: Jul 9, 2014
 *      Author: kustepha
 */

#ifndef MYSTRING_H_
#define MYSTRING_H_

#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>

namespace MyString {


  template <typename T>
  std::string to_string_with_precision(const T a_value, const int n = 6)
  {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
  }

  inline void
  toLower(std::string& s)
  {
    std::transform(s.begin(),s.end(),s.begin(),::tolower);
  }

  inline void
  toUpper(std::string& s)
  {
    std::transform(s.begin(),s.end(),s.begin(),::toupper);
  }

  inline std::vector<std::string>&
  split(const std::string& s, char delim, std::vector<std::string>& elems) {
    if (s.size())
      {
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
      }
    return elems;
  }


  inline std::vector<std::string>
  split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
  }

  inline std::string
  cat(const std::vector<std::string>& src, const std::string& sep = " ")
  {
    if (src.size() == 0)
      return "";

    std::string dest( src[0] );

    for ( uint64_t i = 1; i < src.size(); ++i )
      {
        dest += sep;
        dest += src[i];
      }

    return dest;
  }

  inline uint64_t
  find(const std::vector<std::string>& strings_to_search, const std::string& string_to_find)
  {
    uint64_t idx = 0;
    while(idx < strings_to_search.size())
      if (strings_to_search[idx].compare( string_to_find ) == 0)
        return idx;
    return -1;
  }

  template<typename It>
  uint64_t
  find(It begin, It end, const std::string& string_to_find)
  {
    uint64_t idx = 0;
    while(begin != end)
      if ((*begin++).compare( string_to_find ) == 0)
        return idx;
      else
        ++idx;
    return -1;
  }

  namespace {
    inline int32_t
    sub(const int32_t in, const int32_t toSub, const int32_t sub) {
      return (in == toSub) ? sub : in;
    }
  }

  inline std::vector<int32_t>&
  str2ASCII(const std::string& s, char delim, std::vector<int32_t>& ascii) {

    std::for_each(
        s.begin(),s.end(),
        [&ascii,&delim](const char& c)
        {
      ascii.push_back(sub(static_cast<int32_t>(c),delim,0));
        }
    );

    if (ascii.back() != 0) ascii.push_back(0);

    return ascii;

  }

  inline std::vector<int32_t>
  str2ASCII(const std::string& s, char delim) {

    std::vector<int32_t> ascii;

    str2ASCII(s, delim, ascii);

    return ascii;

  }



  inline std::string&
  ASCII2str(const std::vector<int32_t>& ascii, char delim, std::string& s)
  {

    std::for_each(
        ascii.begin(),ascii.end(),
        [&s,&delim](const int32_t& c)
        {
      s.push_back(sub(c,0,delim));
        }
    );

    return s;

  }

  inline std::string
  ASCII2str(const std::vector<int32_t>& ascii, char delim) {

    std::string s;

    ASCII2str(ascii, delim, s);

    return s;

  }




  inline std::string
  makeLastCharDirSep(const std::string& str)
  {
    std::string s(str);
    if (s.empty())
      s = "./";
    if (s.back() != '/')
      s += "/";
    return s;
  }

  inline std::string
  catDirs(const std::string& a, const std::string& b)
  {
    if (
        (b.compare(0,1,".") == 0)
        || (b.compare(0,1,"/") == 0) )
      return b; // already a full directory

    return makeLastCharDirSep(a).append(b);
  }

  template<typename T>
  inline std::string
  to_bin_string(const T& num)
  {
    constexpr uint64_t n_bits_( sizeof(T)*8 );
    union
    {
      uint64_t num_as_uint_;
      T num_cpy_;
    };
    num_cpy_ = num;
    std::string s(n_bits_,'0');
    for(size_t i = 0; i < n_bits_; ++i)
      s[n_bits_-1-i] = '0' + (  ((num_as_uint_>>i) & 1) );
    return s;
  }


} // namespace

#endif /* MYSTRING_H_ */
