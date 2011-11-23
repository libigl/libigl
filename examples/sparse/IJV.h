// Templates: 
template <typename I,typename V>
struct IJV
{
  // Subscripts
  I i,j;
  // Value
  V v;
  IJV(int i, int j, double v)
  {
    this->i = i;
    this->j = j;
    this->v = v;
  }
  // Compare based on i then j
  bool operator<(const IJV & B) const
  {
    if(this->i < B.i)
    {
      return true;
    }else if(this->i > B.i)
    {
      return false;
    }else
    {
      if(this->j < B.j)
      {
        return true;
      }else if(this->j > B.j)
      {
        return false;
      }else
      {
        // i and j are equal
        return false;
      }
    }
    return false;
  }
};
