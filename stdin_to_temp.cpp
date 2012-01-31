#include "stdin_to_temp.h"

#include <iostream>

IGL_INLINE bool igl::stdin_to_temp(FILE ** temp_file)
{
  // get a temporary file
  *temp_file = tmpfile();
  if(*temp_file == NULL)
  {
    fprintf(stderr,"IOError: temp file could not be created.\n");
    return false;
  }
  char c;
  // c++'s cin handles the stdind input in a reasonable way
  while (std::cin.good())
  {
    c = std::cin.get();
    if(std::cin.good())
    {
      if(1 != fwrite(&c,sizeof(char),1,*temp_file))
      {
        fprintf(stderr,"IOError: error writing to tempfile.\n");
        return false;
      }
    }
  }
  // rewind file getting it ready to read from
  rewind(*temp_file);
  return true;
}
