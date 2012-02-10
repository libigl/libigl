/**
 * Simple example showing how to use stdin_to_temp
 *
 * Compile with:
 *   g++ -O3 -o stdin_to_temp_demo stdin_to_temp_demo.cpp
 *
 * Run examples:
 *   cat file1 | ./stdin_to_temp_demo
 *   cat file1 | ./stdin_to_temp_demo | cat >file2
 *   cat file1 | ./stdin_to_temp_demo dummy1 dummy2 | cat >file2
 *   ./stdin_to_temp_demo <file1 | cat >file2
 *   ./stdin_to_temp_demo <file1 >file2
 *
 */

#include <igl/stdin_to_temp.h>
using namespace igl;

#include <cstdio>
using namespace std;

int main(int argc,char * argv[])
{
  // Process arguements and print to stderr
  for(int i = 1;i<argc;i++)
  {
    fprintf(stderr,"argv[%d] = %s\n",i,argv[i]);
  }

  FILE * temp_file;
  bool success = stdin_to_temp(&temp_file);
  if(!success)
  {
    fprintf(stderr,"Fatal Error: could not convert stdin to temp file\n");
    // try to close temp file
    fclose(temp_file);
    return 1;
  }

  // Do something interesting with the temporary file. 
  // Read the file and write it to stdout
  char c;
  // Read file one character at a time and write to stdout
  while(fread(&c,sizeof(char),1,temp_file)==1)
  {
    fwrite(&c,sizeof(char),1,stdout);
  }
  // close file
  fclose(temp_file);
}
