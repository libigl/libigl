#define VERBOSE
#include <igl/pathinfo.h>
#include <igl/C_STR.h>
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

  const char * format = "%-25s | %-10s %-10s %-10s %-10s\n";
  printf(format,
    "input",
    "dirname",
    "basename",
    "extension",
    "filename");
  string line;
  while( getline(fin, line) )
  {
    string dirname,basename,extension,filename;
    pathinfo(line,dirname,basename,extension,filename);
    printf(format,
      C_STR("\""<<line<<"\""),
      C_STR("\""<<dirname<<"\""),
      C_STR("\""<<basename<<"\""),
      C_STR("\""<<extension<<"\""),
      C_STR("\""<<filename<<"\""));
  }
}
