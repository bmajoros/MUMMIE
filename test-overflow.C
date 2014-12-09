/*
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
*/
#include <iostream>
#include <fenv.h>
using namespace std;
int main(int argc,char **argv) {
  feenableexcept(FE_OVERFLOW);
  double d=100;
  for(int i=0 ; i<10 ; ++i) {
    float x=0;
    float y=x/x;
    cout<<"y="<<y<<endl;
    d=d*d;
    cout<<d<<endl;
    if(fetestexcept(FE_OVERFLOW)) cout<<"OVERFLOW!"<<endl;
  }
}

