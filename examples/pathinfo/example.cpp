#define VERBOSE
#include <igl/pathinfo.h>
using namespace igl;
#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char * argv[])
{
  ifstream fin("input.txt");
  if (fin.is_open() == false)
  {
    // error
    return 1;
  }

  string line;
  while( getline(fin, line) )
  {
    string dirname,basename,extension,filename;
    pathinfo(line,dirname,basename,extension,filename);
    printf("%s -> %s,%s,%s,%s\n",
      line.c_str(),
      dirname.c_str(),
      basename.c_str(),
      extension.c_str(),
      filename.c_str());
  }
}
