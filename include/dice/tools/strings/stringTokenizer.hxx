#ifndef __STRING_TOKENIZER_HXX__
#define __STRING_TOKENIZER_HXX__

#include "../../internal/namespace.header"

class StringTokenizer
{
public:
  template <class T>
  static bool from_string(T &t, const std::string &s)
  {
    std::istringstream iss(s);
    return !(iss >> t).fail();
  }

  template <class T>
  static void to_string(const T t, std::string &s)
  {
    std::ostringstream oss;
    oss << t;
    s = oss.str();
  }

  template <typename Container>
  static Container &split(Container &result,
                          const typename Container::value_type &s,
                          const typename Container::value_type &delimiters)
  {
    result.clear();
    size_t current;
    size_t next = -1;
    do
    {
      next = s.find_first_not_of(delimiters, next + 1);
      if (next == Container::value_type::npos)
        break;
      next -= 1;
      current = next + 1;
      next = s.find_first_of(delimiters, current);
      result.push_back(s.substr(current, next - current));
    } while (next != Container::value_type::npos);
    return result;
  }
};

#include "../../internal/namespace.footer"
#endif
