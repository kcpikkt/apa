#define CATCH_CONFIG_MAIN
#include "include/catch2/catch.hpp"
#include "include/gmp-6.1.2/gmp.h"

#include <stdint.h>

#include <iostream>
#include <array>
#include <random>
#include <bitset>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wkeyword-macro"
#define private public // :)
#pragma GCC diagnostic pop

#define APA_IMPL_TYPE uint16_t // internal word size
#include "bigint.h"

template<typename ...Ts> void print(Ts... ts){
    (void) std::initializer_list<int>{ (std::cout << ts, 0)...}; }
template<typename ...Ts> void log(Ts... ts){ print(ts..., '\n'); }

#define TIMES(N) for(uint32_t i=0; i < N; i++)
[[maybe_unused]] constexpr uint32_t large_test_repeats = 1000;
[[maybe_unused]] constexpr uint32_t medium_test_repeats = 100;
[[maybe_unused]] constexpr uint32_t small_test_repeats = 20;

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define S_TYPE(N) "apa::s<" STR(N) ">"

namespace{
    using impl_t = APA_IMPL_TYPE;
    [[maybe_unused]] constexpr size_t impl_t_size = sizeof(impl_t);
    [[maybe_unused]] constexpr size_t impl_t_bits = sizeof(impl_t)*8;
}

template<size_t SZ, typename T>
bool equal(apa::s<SZ> bint, T tint) {
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

template<typename T, size_t ARR_SZ>
void init_random_data(std::array<T, ARR_SZ>& data) {
    for(size_t i=0; i<data.size(); i++) 
        data[i] = static_cast<T>(mt64());
}

template<size_t SZ, typename T, size_t ARR_SZ>
apa::s<SZ> create_signed(std::array<T, ARR_SZ>& data, bool sign = 0) {
    apa::s<SZ> ret;
    REQUIRE(ret.import(data.data(), data.size()));
    ret.set_sign_bool(sign);
    return ret;
}

template<typename T, size_t ARR_SZ>
void create_mpz(mpz_t * mpz, std::array<T, ARR_SZ>& data, bool sign = 0) {
    mpz_init(*mpz);
    mpz_import(*mpz, data.size(), -1, sizeof(T), 0, 0, data.data());
    if(sign) mpz_neg(*mpz, *mpz);
}

template<size_t SZ>
apa::s<SZ> mpz_to_signed(mpz_t * mpz) {
    constexpr size_t ARR_SZ = (SZ / (sizeof(uint64_t) * 8)) << 1;
    apa::s<SZ> ret;
    std::array<uint64_t, ARR_SZ> data = {0};
    size_t count = 0;
    mpz_export(data.data(), &count, -1, sizeof(uint64_t), 0, 0, *mpz);
    REQUIRE(ret.import(data.data(), data.size()));
    ret.set_sign( mpz_sgn(*mpz));
    return ret;
}

TEST_CASE( "Constructors" ) {

    SECTION( "Zero initialized when no parameters provided" ) {
        apa::s<128> bint;
        REQUIRE(equal(bint, 0));
    }

    SECTION( "64 bit int to 128 bit apa" ) {
        uint64_t tint = mt64(); 

        apa::s<128> bint(tint);

        // REQUIRE(equal(bint, tint));
    }

    SECTION( "64 bit int to 64 bit apa" ) {
        uint64_t tint = mt64();

        apa::s<64> bint(tint);

        // REQUIRE(equal(bint, tint));
    }

    SECTION( "64 bit int to 32 bit apa loseless" ) {
        uint64_t tint = mt32();

        apa::s<32> bint(tint);

        // REQUIRE(bint.trucated() == false);
        // REQUIRE(equal(bint, tint));
    }

    SECTION( "64 bit int to 32 bit apa with truncation" ) {
        uint64_t tint = (size_t)(mt32() + 1) << 32;

        apa::s<32> bint(tint);

        // REQUIRE(bint.trucated() == true);
        // REQUIRE(equal(bint, tint));
    }
}

TEST_CASE( "Basics" ) {

    SECTION("Setting sign") {
        apa::s<1024> bint;
        bint.set_sign_bool(0);
        REQUIRE(bint.sign() ==  1);

        bint.negate();
        REQUIRE(bint.sign() == -1);

        bint.set_sign_bool(1);
        REQUIRE(bint.sign() == -1);

        bint.negate();
        REQUIRE(bint.sign() ==  1);

        bint.set_sign(1);
        REQUIRE(bint.sign() ==  1);

        bint.set_sign(-1);
        REQUIRE(bint.sign() == -1);
    }

}

TEST_CASE( "Comparison" ) {

    SECTION( "Equality" ) {
        constexpr size_t TEST_BIT_SZ = 1024;
        constexpr size_t DATA_ARR_SZ = TEST_BIT_SZ / (sizeof(uint64_t) * 8);

        SECTION( "operator==(" S_TYPE(TEST_BIT_SZ) ", " S_TYPE(TEST_BIT_SZ) "), " \
                 "assigned the same date - equal" )
        {
            TIMES(medium_test_repeats) {
                bool sign = mt32() & 1;

                std::array<uint64_t, DATA_ARR_SZ> datain;
                init_random_data(datain);

                auto bint1 = create_signed<TEST_BIT_SZ>(datain, sign);
                auto bint2 = create_signed<TEST_BIT_SZ>(datain, sign);

                REQUIRE(bint1 == bint2);
            }
        }


        SECTION( "operator!=(" S_TYPE(TEST_BIT_SZ) ", " S_TYPE(TEST_BIT_SZ) "), " \
                 "assigned different data - not equal" )
        {
            TIMES(medium_test_repeats) {
                bool sign = mt32() & 1;

                std::array<uint64_t, DATA_ARR_SZ> datain1, datain2;
                init_random_data(datain1); init_random_data(datain2);

                auto bint1 = create_signed<TEST_BIT_SZ>(datain1, sign);
                auto bint2 = create_signed<TEST_BIT_SZ>(datain2, sign);

                REQUIRE(bint1 != bint2);
            }
        }


        SECTION( "operator!=(" S_TYPE(TEST_BIT_SZ) ", " S_TYPE(TEST_BIT_SZ) "), " \
                 "assigned slightly different data - not equal" )
        {
            TIMES(medium_test_repeats) {
                std::array<uint64_t, DATA_ARR_SZ> datain;
                init_random_data(datain);

                auto bint1 = create_signed<TEST_BIT_SZ>(datain);
                datain[mt32()%16] += 1;
                auto bint2 = create_signed<TEST_BIT_SZ>(datain);

                REQUIRE(bint1 != bint2);
            }
        }

        SECTION( "operator!=(" S_TYPE(TEST_BIT_SZ) ", " S_TYPE(TEST_BIT_SZ) "), " \
                 "additive inverses - not equal" )
        {
            TIMES(medium_test_repeats) {
                std::array<uint64_t, DATA_ARR_SZ> datain;
                init_random_data(datain);

                auto bint1 = create_signed<TEST_BIT_SZ>(datain, true);
                auto bint2 = create_signed<TEST_BIT_SZ>(datain, false);

                REQUIRE(bint1 != bint2);
            }
        }
    }

    SECTION("Inequality") {
        constexpr size_t TEST_BIT_SZ = 1024;
        constexpr size_t DATA_ARR_SZ = TEST_BIT_SZ / (sizeof(uint64_t) * 8);

        SECTION( "operator>(" S_TYPE(TEST_BIT_SZ) ", " S_TYPE(TEST_BIT_SZ) ") against gmp" ) {
            TIMES(medium_test_repeats) {
                std::array<uint64_t, DATA_ARR_SZ> datain1, datain2;
                init_random_data(datain1); init_random_data(datain2);

                bool sign1 = mt32() & 1, sign2 = mt32() & 1;

                auto bint1 = create_signed<1024>(datain1, sign1);
                auto bint2 = create_signed<1024>(datain2, sign2);

                bool bint_result = bint1 > bint2;

                mpz_t gmpint1, gmpint2;
                create_mpz(&gmpint1, datain1, sign1);
                create_mpz(&gmpint2, datain2, sign2);

                int gmpint_result = mpz_cmp(gmpint1, gmpint2);

                if(gmpint_result > 0) {
                    REQUIRE( bint_result);
                } else if(gmpint_result < 0) {
                    REQUIRE(!bint_result);
                }
            }
        }

        SECTION( "operator>(" S_TYPE(TEST_BIT_SZ) ", " S_TYPE(TEST_BIT_SZ) ") against gmp" ) {
            TIMES(medium_test_repeats) {
                std::array<uint64_t, DATA_ARR_SZ> datain1, datain2;
                init_random_data(datain1); init_random_data(datain2);

                bool sign1 = mt32() & 1, sign2 = mt32() & 1;

                auto bint1 = create_signed<1024>(datain1, sign1);
                auto bint2 = create_signed<1024>(datain2, sign2);

                bool bint_result = bint1 < bint2;

                mpz_t gmpint1, gmpint2;
                create_mpz(&gmpint1, datain1, sign1);
                create_mpz(&gmpint2, datain2, sign2);

                int gmpint_result = mpz_cmp(gmpint1, gmpint2);

                if(gmpint_result > 0) {
                    REQUIRE(!bint_result);
                } else if(gmpint_result < 0) {
                    REQUIRE( bint_result);
                }
            }
        }
    }
}


template<size_t B1, size_t B2>
void addition_test()
{
    constexpr size_t DATA_ARR_SZ1 = B1 / (sizeof(uint64_t) * 8);
    constexpr size_t DATA_ARR_SZ2 = B2 / (sizeof(uint64_t) * 8);

    std::array<uint64_t, DATA_ARR_SZ1> datain1;
    std::array<uint64_t, DATA_ARR_SZ2> datain2;
    init_random_data(datain1); init_random_data(datain2);

    bool sign1 = mt32() & 1, sign2 = mt32() & 1;

    auto bint1 = create_signed<B1>(datain1, sign1);
    auto bint2 = create_signed<B2>(datain2, sign2);

    mpz_t gmpint1, gmpint2;
    create_mpz(&gmpint1, datain1, sign1);
    create_mpz(&gmpint2, datain2, sign2);

    // auto bint_result = apa::s<1024>::operator+(bint1, bint2, bint1::_A);
    auto bint_result = bint1 + bint2;

    mpz_t gmpint_result;
    mpz_init(gmpint_result);
    mpz_add(gmpint_result, gmpint1, gmpint2);

    REQUIRE(mpz_to_signed<std::max(B1, B2) + 1>(&gmpint_result) == bint_result);
}

template<size_t B1, size_t B2>
void subtraction_test()
{
    constexpr size_t DATA_ARR_SZ1 = B1 / (sizeof(uint64_t) * 8);
    constexpr size_t DATA_ARR_SZ2 = B2 / (sizeof(uint64_t) * 8);

    std::array<uint64_t, DATA_ARR_SZ1> datain1;
    std::array<uint64_t, DATA_ARR_SZ2> datain2;
    init_random_data(datain1); init_random_data(datain2);

    bool sign1 = mt32() & 1, sign2 = mt32() & 1;

    auto bint1 = create_signed<B1>(datain1, sign1);
    auto bint2 = create_signed<B2>(datain2, sign2);

    mpz_t gmpint1, gmpint2;
    create_mpz(&gmpint1, datain1, sign1);
    create_mpz(&gmpint2, datain2, sign2);

    auto bint_result = bint1 - bint2;

    mpz_t gmpint_result;
    mpz_init(gmpint_result);
    mpz_sub(gmpint_result, gmpint1, gmpint2);

    REQUIRE(mpz_to_signed<std::max(B1, B2) + 1>(&gmpint_result) == bint_result);
}

template<size_t B1, size_t B2>
void multiplication_test()
{
    constexpr size_t DATA_ARR_SZ1 = B1 / (sizeof(uint64_t) * 8);
    constexpr size_t DATA_ARR_SZ2 = B2 / (sizeof(uint64_t) * 8);

    std::array<uint64_t, DATA_ARR_SZ1> datain1;
    std::array<uint64_t, DATA_ARR_SZ2> datain2;
    init_random_data(datain1); init_random_data(datain2);

    bool sign1 = mt32() & 1, sign2 = mt32() & 1;

    std::string str;

    auto bint1 = create_signed<B1>(datain1, sign1);
    auto bint2 = create_signed<B2>(datain2, sign2);

    // str = bint1.binary_string();
    // TIMES(64) print(str[i]); log("ok?");
    // log(std::bitset<64>(datain1[DATA_ARR_SZ1 -1]));

    // str = bint2.binary_string();
    // TIMES(64) print(str[i]); log("ok?");
    // log(std::bitset<64>(datain1[DATA_ARR_SZ2 -1]));

    mpz_t gmpint1, gmpint2;
    create_mpz(&gmpint1, datain1, sign1);
    create_mpz(&gmpint2, datain2, sign2);

    auto bint_result = bint1 * bint2;

    mpz_t gmpint_result;
    mpz_init(gmpint_result);
    mpz_mul(gmpint_result, gmpint1, gmpint2);

    // str = bint_result.binary_string();
    // TIMES(100) print(str[i]); log("ok?");

    // str = mpz_to_signed<B1 + B2>(&gmpint_result).binary_string();
    // TIMES(100) print(str[i]); log("ok?");
    // log();

    auto prt = [](apa::s<B1 + B2> bint){
                   auto str = bint.binary_string();
                   TIMES(64) print(str[i]); log();
               };
    prt(mpz_to_signed<B1 + B2>(&gmpint_result) );
    prt(bint_result);

    REQUIRE(mpz_to_signed<B1 + B2>(&gmpint_result) == bint_result);
}

template<size_t B>
void bitshift_left_test()
{
    constexpr size_t TEST_BIT_SZ = B;
    constexpr size_t DATA_ARR_SZ = TEST_BIT_SZ / (sizeof(uint64_t) * 8);

    std::array<uint64_t, DATA_ARR_SZ> datain1, datain2;
    init_random_data(datain1); init_random_data(datain2);

    auto bint1 = create_signed<TEST_BIT_SZ>(datain1, false);
    auto shift = mt64() % TEST_BIT_SZ;

    mpz_t gmpint1;
    create_mpz(&gmpint1, datain1, false);

    auto bint_result = bint1 << shift;

    mpz_t gmpint_result;
    mpz_init(gmpint_result);
    mpz_mul_2exp(gmpint_result, gmpint1, shift);

    REQUIRE(mpz_to_signed<TEST_BIT_SZ>(&gmpint_result) == bint_result);
}

template<size_t B>
void bitshift_right_test()
{
    constexpr size_t DATA_ARR_SZ = B / (sizeof(uint64_t) * 8);

    std::array<uint64_t, DATA_ARR_SZ> datain1, datain2;
    init_random_data(datain1); init_random_data(datain2);

    auto bint1 = create_signed<B>(datain1, false);
    auto shift = mt64() % B;

    mpz_t gmpint1;
    create_mpz(&gmpint1, datain1, false);

    auto bint_result = bint1 >> shift;

    mpz_t gmpint_result;
    mpz_init(gmpint_result);
    mpz_fdiv_q_2exp(gmpint_result, gmpint1, shift);

    REQUIRE(mpz_to_signed<B>(&gmpint_result) == bint_result);
}

TEST_CASE( "Arithmentic" ) {

    SECTION( "Addition" ) {

        SECTION( "operator+(apa::s<SZ1>, apa::s<SZ2>) against gmp" ) {

            TIMES(medium_test_repeats) addition_test<1024,1024>();
            TIMES(medium_test_repeats) addition_test<1024,512>();
            TIMES(medium_test_repeats) addition_test<4096,1024>();
        }
    }

    SECTION( "Subtraction" ) {

        SECTION( "operator-(apa::s<SZ1>, apa::s<SZ2>) against gmp" ) {

            TIMES(medium_test_repeats) subtraction_test<1024,1024>();
            TIMES(medium_test_repeats) subtraction_test<1024,512>();
            TIMES(medium_test_repeats) subtraction_test<4096,1024>();
        }
    }

    SECTION( "Multiplication" ) {

        SECTION( "operator*(apa::s<SZ1>, apa::s<SZ2>) against gmp" ) {

            TIMES(medium_test_repeats) multiplication_test<1024,1024>();
            // TIMES(medium_test_repeats) multiplication_test<1024,512>();
            // TIMES(medium_test_repeats) multiplication_test<4096,1024>();
        }
    }

    SECTION( "Bitshifts" ) {

        SECTION( "operator<<(apa::s<SZ1>, size_t) against gmp" ) {

            TIMES(medium_test_repeats) bitshift_left_test<1024>();
            TIMES(medium_test_repeats) bitshift_left_test<4096>();
            TIMES(medium_test_repeats) bitshift_left_test<512>();
        }

        SECTION( "operator>>(apa::s<SZ1>, size_t) against gmp" ) {

            TIMES(medium_test_repeats) bitshift_right_test<1024>();
            TIMES(medium_test_repeats) bitshift_right_test<4096>();
            TIMES(medium_test_repeats) bitshift_right_test<512>();
        }
    }
}
