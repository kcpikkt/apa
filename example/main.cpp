#include <stdio.h>
#include <initializer_list>
#include <iostream>
#include <stdint.h>
#include <iomanip>
#include <bitset>
#include <random>
#define PRINT_VAR(v) printf(#v "= %i\n",v)

#include "bigint.h"

template<typename ...Ts> void print(Ts... ts){
    (void) std::initializer_list<int>{ (std::cout << ts, 0)...}; }
template<typename ...Ts> void log(Ts... ts){ print(ts..., '\n'); }

template<size_t N> void print_bigint(const bigint::Signed<N>& integer){
    uint64_t num = 0;
    for(uint64_t i=0; i<integer.bit_sz; i++){
        num += integer.bit_at(i) * pow(2,i);
    }
    print(num, '\n');
}

// #define LHS 11
#define B 16

const long LHS = 251;
const long RHS = 1801;
//               4294967296
int main(){
    std::random_device rd;
    std::mt19937_64 mt(rd());

    bigint::s<B> lhs = (LHS);
    bigint::s<32> rhs = (RHS);
    auto result = lhs * rhs;

    log(std::bitset<2*B>(LHS*RHS), " expected");

    log(result.binary_string(), " result");
    log(std::bitset<B>(LHS), " * ", std::bitset<B>(RHS));
    log(lhs.binary_string(), " * ", rhs.binary_string());
    log("s = ", result.bit_sz);
}
