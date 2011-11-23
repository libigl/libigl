#!/bin/bash
//usr/bin/tail -n +2 $0 | g++ -o main -x c++ - && ./main && rm main && exit
#include "IJV.h"

#include <cstdio>
#include <vector>
#include <algorithm>
using namespace std;

int print_ijv(const IJV<int,double> & ijv)
{
  printf("%d %d %g\n",
    ijv.i,
    ijv.j,
    ijv.v);
}

int main(int argc, char * argv[])
{
  vector<IJV<int,double> > Aijv;
  Aijv.push_back(IJV<int,double>(1,2,10));
  Aijv.push_back(IJV<int,double>(4,2,10));
  Aijv.push_back(IJV<int,double>(4,3,10));
  Aijv.push_back(IJV<int,double>(9,2,10));
  Aijv.push_back(IJV<int,double>(1,2,10));
  Aijv.push_back(IJV<int,double>(1,3,10));
  Aijv.push_back(IJV<int,double>(1,1,10));
  Aijv.push_back(IJV<int,double>(3,2,10));

  printf("Original:\n");
  for_each(Aijv.begin(),Aijv.end(),print_ijv);
  sort(Aijv.begin(),Aijv.end());
  printf("Sorted:\n");
  for_each(Aijv.begin(),Aijv.end(),print_ijv);

  return 0;
}
