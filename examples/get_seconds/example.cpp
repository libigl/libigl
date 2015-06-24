#include <igl/get_seconds.h>
#include <cstdio>
#include <cmath>

int main(int argc, char * argv[])
{
  using namespace igl;
  using namespace std;
  double start = get_seconds();
  printf("start: %lgs\n",start);
  double lap = start;
  double now;
  do
  {
    now = get_seconds();
    if((now-lap)>=1.0)
    {
      printf("%lgs - %lgs = %lgs\n",now,start,now-start);
      lap = now-((now-start)-floor(now-start));
    }
  }while((now - start)<10.0);
  printf("end: %lgs\n",now);
  return 0;
}
