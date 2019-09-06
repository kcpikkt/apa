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

#include "bigint.h"

template<typename ...Ts> void print(Ts... ts){
    (void) std::initializer_list<int>{ (std::cout << ts, 0)...}; }
template<typename ...Ts> void log(Ts... ts){ print(ts..., '\n'); }

int main(){
    [[maybe_unused]]auto d0 = apa::floating_point_t_decompose(3.0f);
}
