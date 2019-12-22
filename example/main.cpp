#include <stdio.h>
#include <initializer_list>
#include <type_traits>
#include <iostream>
#include <stdint.h>
#include <iomanip>
#include <bitset>
#include <random>
#include <tgmath.h>
#define PRINT_VAR(v) printf(#v "= %i\n",v)


template<typename ...Ts> void print(Ts... ts){
    (void) std::initializer_list<int>{ (std::cout << ts, 0)...}; }
template<typename ...Ts> void log(Ts... ts){ print(ts..., '\n'); }

#include "apa.h"


int main(){
    apa::s<128> a = 10;
    apa::s<64> b = 172914.78;

    // width of c is going to be 64 + 128 = 192
    auto c = a * b;
    // well, no division no printing....
    // tests against gmp are passing though!
}
