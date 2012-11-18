#ifndef IGL_INDEXCOMPARISON_H
#define IGL_INDEXCOMPARISON_H
// Comparison struct used by sort
// http://bytes.com/topic/c/answers/132045-sort-get-index
namespace igl{

  // For use with functions like std::sort
  template<class T> struct IndexLessThan
  {
    IndexLessThan(const T arr) : arr(arr) {}
    bool operator()(const size_t a, const size_t b) const
    {
      return arr[a] < arr[b];
    }
    const T arr;
  };

  // For use with functions like std::unique
  template<class T> struct IndexEquals
  {
    IndexEquals(const T arr) : arr(arr) {}
    bool operator()(const size_t a, const size_t b) const
    {
      return arr[a] == arr[b];
    }
    const T arr;
  };
}

#endif
