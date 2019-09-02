#define CATCH_CONFIG_MAIN
#include "include/catch2/catch.hpp"
#include "include/gmp-6.1.2/gmp.h"

#include <stdint.h>

#include <iostream>
#include <random>
#include <bitset>

#define private public // :)
#include "bigint.h"

template<typename ...Ts> void print(Ts... ts){
    (void) std::initializer_list<int>{ (std::cout << ts, 0)...}; }
template<typename ...Ts> void log(Ts... ts){ print(ts..., '\n'); }

#define TIMES(N) for(int i=0; i < N; i++)


namespace{
    using impl_t = BIGINT_IMPL_TYPE;
    [[maybe_unused]] constexpr size_t impl_t_byte_sz = sizeof(impl_t);
    [[maybe_unused]] constexpr size_t impl_t_bit_sz = sizeof(impl_t)*8;
}

template<size_t SZ, typename T>
bool equal(bigint::s<SZ> bint, T tint) {
    size_t i = 0;
    for(; i < bint.segments_count && i < sizeof(T) / sizeof(impl_t); i++){
        impl_t mask = static_cast<impl_t> (tint >> i * sizeof(impl_t) * 8);

        if(! (bint.get_segment(i) == mask)) return false;
    }
    for(; i < bint.segments_count; i++){
        if(! (bint.get_segment(i) == 0)) return false;
    }
    return true;
}

std::random_device rd;
std::mt19937    mt32(rd());
std::mt19937_64 mt64(rd());

TEST_CASE( "Constructors" ) {

    SECTION( "Zero initialized when no parameters provided" ) {
        bigint::s<128> bint;
        REQUIRE(equal(bint, 0));
    }

    SECTION( "64 bit int to 128 bit bigint" ) {
        uint64_t tint = mt64(); 

        bigint::s<128> bint(tint);

        REQUIRE(equal(bint, tint));
    }

    SECTION( "64 bit int to 64 bit bigint" ) {
        uint64_t tint = mt64();

        bigint::s<64> bint(tint);

        REQUIRE(equal(bint, tint));
    }

    SECTION( "64 bit int to 32 bit bigint loseless" ) {
        uint64_t tint = mt32();

        bigint::s<32> bint(tint);

        REQUIRE(bint.was_truncated() == false);
        REQUIRE(equal(bint, tint));
    }

    SECTION( "64 bit int to 32 bit bigint with truncation" ) {
        uint64_t tint = (size_t)(mt32() + 1) << 32;

        bigint::s<32> bint(tint);

        REQUIRE(bint.was_truncated() == true);
        REQUIRE(equal(bint, tint));
    }
}

TEST_CASE( "Comparison" ) {

    SECTION( "Two bigints assinged the same data are equal" ) {
        TIMES(100) {
            uint64_t datain[16];

            for(auto& d : datain) d = mt64();

            bigint::s<1024> bint1;
            bigint::s<1024> bint2;
            REQUIRE(bint1.import(datain, 16));
            REQUIRE(bint2.import(datain, 16));

            REQUIRE(bint1 == bint2);
        }
    }
    SECTION( "Two bigints assinged different data are not equal" ) {
        TIMES(100) {
            uint64_t datain1[16];
            uint64_t datain2[16];

            for(auto& d : datain1) d = mt64();
            for(auto& d : datain2) d = mt64();

            bigint::s<1024> bint1;
            bigint::s<1024> bint2;
            REQUIRE(bint1.import(datain1, 16));
            REQUIRE(bint2.import(datain2, 16));

            REQUIRE(bint1 != bint2);
        }
    }

    SECTION( "Two bigints assinged slightly different data are not equal" ) {
        TIMES(100) {
            uint64_t datain[16];

            for(auto& d : datain) d = mt64();

            bigint::s<1024> bint1;
            bigint::s<1024> bint2;

            REQUIRE(bint1.import(datain, 16));

            datain[8] += 1;
            REQUIRE(bint2.import(datain, 16));

            REQUIRE(bint1 != bint2);
        }
    }
}

TEST_CASE( "Addition" ) {

    SECTION( "64 bit range random operator+(bigint, bigint)" ) {
        TIMES(1000) {
            uint64_t tint1 = mt64() >> 1;
            uint64_t tint2 = mt64() >> 1;
            uint64_t tint_result = tint1 + tint2;

            bigint::s<128> bint1(tint1);
            bigint::s<128> bint2(tint2);
            bigint::s<128> bint_result = bint1 + bint2;

            REQUIRE(equal(bint_result, tint_result));
        }
    }
    SECTION( "64 bit range random operator+=(bigint, bigint)" ) {
        TIMES(1000) {
            uint64_t tint1 = mt64() >> 1;
            uint64_t tint2 = mt64() >> 1;
            uint64_t tint_result = tint1 + tint2;

            bigint::s<128> bint1(tint1);
            bigint::s<128> bint2(tint2);
            bint1 += bint2;

            REQUIRE(equal(bint1, tint_result));
        }
    }
}

TEST_CASE( "Subtraction" ) {
}

TEST_CASE( "Multiplication" ) {

    SECTION( "thousand 64 bit range random operator*(bigint, bigint) tests" ) {
        TIMES(1000) {
            uint64_t tint1 = mt32() >> 1;
            uint64_t tint2 = mt32() >> 1;
            uint64_t tint_result = tint1 * tint2;

            bigint::s<64> bint1(tint1);
            bigint::s<64> bint2(tint2);
            auto bint_result = bint1 * bint2;

            REQUIRE(equal(bint_result, tint_result));
        }
    }
    SECTION( "1024 bit range random operator*(bigint, bigint) with gmp" ) {
        TIMES(1000) {
            uint64_t datain1[16];
            uint64_t datain2[16];

            for(auto& d : datain1) d = mt64();
            for(auto& d : datain2) d = mt64();

            bigint::s<1024> bint1;
            bigint::s<1024> bint2;
            REQUIRE(bint1.import(datain1, 16));
            REQUIRE(bint2.import(datain2, 16));

            mpz_t gmpint1;
            mpz_t gmpint2;
            mpz_init(gmpint1);
            mpz_init(gmpint2);
            mpz_import(gmpint1, 16, -1, sizeof(uint64_t), 0, 0, datain1);
            mpz_import(gmpint2, 16, -1, sizeof(uint64_t), 0, 0, datain2);

            bigint::s<2048> bint_result = bint1 * bint2;

            mpz_t gmpint_result;
            mpz_init(gmpint_result);
            mpz_mul(gmpint_result, gmpint1, gmpint2);
            uint64_t gmpint_result_data[32];
            size_t count = 0;
            mpz_export(gmpint_result_data, &count, -1, sizeof(uint64_t), 0, 0,
                       gmpint_result);

            bigint::s<2048> gmpint_result_bint;
            REQUIRE(gmpint_result_bint.import(gmpint_result_data, 32));

            REQUIRE(gmpint_result_bint == bint_result);
        }
    }
    SECTION( "1048576 bit range random operator*(bigint, bigint) with gmp" ) {
        TIMES(0) {
            uint64_t * datain1 = new uint64_t[16384];
            uint64_t * datain2 = new uint64_t[16384];

            for(int i=0; i<16384; i++) *(datain1 + i) = mt64();
            for(int i=0; i<16384; i++) *(datain2 + i) = mt64();

            bigint::s<1048576ull> * bint1 = new bigint::s<1048576ull>();
            bigint::s<1048576ull> * bint2 = new bigint::s<1048576ull>();
            REQUIRE(bint1->import(datain1, 16384));
            REQUIRE(bint2->import(datain2, 16384));

            mpz_t gmpint1;
            mpz_t gmpint2;
            mpz_init(gmpint1);
            mpz_init(gmpint2);
            mpz_import(gmpint1, 16384, -1, sizeof(uint64_t), 0, 0, datain1);
            mpz_import(gmpint2, 16384, -1, sizeof(uint64_t), 0, 0, datain2);

            bigint::s<2097152ull> * bint_result = new bigint::s<2097152ull>();
            // (*bint_result) = (*bint1) * (*bint2);

            mpz_t gmpint_result;
            mpz_init(gmpint_result);
            mpz_mul(gmpint_result, gmpint1, gmpint2);
            uint64_t * gmpint_result_data = new uint64_t[32768];
            size_t count = 0;
            // mpz_export(gmpint_result_data, &count, -1, sizeof(uint64_t), 0, 0,
            //            gmpint_result);

            log(count, "count");

            // bigint::s<2097152> gmpint_result_bint = bint1 * bint2;
            // REQUIRE(gmpint_result_bint.import(gmpint_result_data, 32768));

            // REQUIRE(gmpint_result_bint == bint_result);
        }
    }
}
