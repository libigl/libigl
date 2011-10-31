#include "readOBJ.h"

#include <cstdio>
#include <vector>

using namespace igl;
using namespace std;

// Template:
//   T  type that can be safely cast to float
// Inputs:
//   vv  vector of vectors of type T
template <typename T>
void print_vector_of_vectors_as_floats(const std::vector<std::vector<T > > & vv)
{
  for(int i = 0;i<vv.size();i++)
  {
    for(int j = 0;j<vv[i].size();j++)
    {
      printf("%g ",(float)(vv[i][j]));
    }
    printf("\n");
  }
}

int main(int argc, char * argv[])
{
  if(argc <= 1)
  {
    printf("USAGE:\n  ./example [path_1] [path_2] ... [path_n]\n");
    return 1;
  }
  std::vector<std::vector<double> > V,TC,N;
  std::vector<std::vector<int> > F,FTC,FN;
  // loop over arguments
  for(int i = 1; i < argc; i++)
  {
    if(i != 1)
    {
      printf("-----------------------------------------------------------\n");
    }
    readOBJ(argv[i],V,TC,N,F,FTC,FN);
    cout<<"V=[";  print_vector_of_vectors_as_floats(V);  cout<<"];"<<endl;
    cout<<"TC=["; print_vector_of_vectors_as_floats(TC); cout<<"];"<<endl;
    cout<<"N=[";  print_vector_of_vectors_as_floats(N);  cout<<"];"<<endl;
    cout<<"F=[";  print_vector_of_vectors_as_floats(F);  cout<<"];"<<endl;
    cout<<"FTC=[";print_vector_of_vectors_as_floats(FTC);cout<<"];"<<endl;
    cout<<"FN=["; print_vector_of_vectors_as_floats(FN); cout<<"];"<<endl;
  }
  return 0;
}
