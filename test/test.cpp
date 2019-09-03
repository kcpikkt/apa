#define CATCH_CONFIG_MAIN
#include "include/catch2/catch.hpp"
#include "include/gmp-6.1.2/gmp.h"

#include <stdint.h>

#include <iostream>
#include <array>
#include <random>
#include <bitset>

#define private public // :)
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
#define S_TYPE(N) "bigint::s<" STR(N) ">"

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

template<typename T, size_t ARR_SZ>
void init_random_data(std::array<T, ARR_SZ>& data) {
    for(size_t i=0; i<data.size(); i++) 
        data[i] = static_cast<T>(mt64());
}

template<size_t SZ, typename T, size_t ARR_SZ>
bigint::s<SZ> create_signed(std::array<T, ARR_SZ>& data, bool sign = 0) {
    bigint::s<SZ> ret;
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
bigint::s<SZ> mpz_to_signed(mpz_t * mpz) {
    constexpr size_t ARR_SZ = (SZ / (sizeof(uint64_t) * 8)) << 1;
    bigint::s<SZ> ret;
    std::array<uint64_t, ARR_SZ> data = {0};
    size_t count = 0;
    mpz_export(data.data(), &count, -1, sizeof(uint64_t), 0, 0, *mpz);
    REQUIRE(ret.import(data.data(), data.size()));
    ret.set_sign( mpz_sgn(*mpz));
    return ret;
}

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

        REQUIRE(bint.trucated() == false);
        REQUIRE(equal(bint, tint));
    }

    SECTION( "64 bit int to 32 bit bigint with truncation" ) {
        uint64_t tint = (size_t)(mt32() + 1) << 32;

        bigint::s<32> bint(tint);

        REQUIRE(bint.trucated() == true);
        REQUIRE(equal(bint, tint));
    }
}

TEST_CASE( "Basics" ) {

    SECTION("Setting sign") {
        bigint::s<1024> bint;
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

                // bool sign1 = mt32() & 1, sign2 = mt32() & 1;
                bool sign1 = false, sign2 = false;

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

                // bool sign1 = mt32() & 1, sign2 = mt32() & 1;
                bool sign1 = false, sign2 = false;

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

TEST_CASE( "Addition" ) {
    constexpr size_t TEST_BIT_SZ = 1024;
    constexpr size_t DATA_ARR_SZ = TEST_BIT_SZ / (sizeof(uint64_t) * 8);

    SECTION( "operator+(bigint::s<1024>, bigint::s<1024>) against gmp" ) {
        TIMES(large_test_repeats) {
            std::array<uint64_t, DATA_ARR_SZ> datain1, datain2;
            init_random_data(datain1); init_random_data(datain2);

            bool sign1 = mt32() & 1, sign2 = mt32() & 1;

            auto bint1 = create_signed<TEST_BIT_SZ>(datain1, sign1);
            auto bint2 = create_signed<TEST_BIT_SZ>(datain2, sign2);

            mpz_t gmpint1, gmpint2;
            create_mpz(&gmpint1, datain1, sign1);
            create_mpz(&gmpint2, datain2, sign2);

            auto bint_result = bint1 + bint2;

            mpz_t gmpint_result;
            mpz_init(gmpint_result);
            mpz_add(gmpint_result, gmpint1, gmpint2);

            REQUIRE(mpz_to_signed<TEST_BIT_SZ + 1>(&gmpint_result) == bint_result);
        }
    }
}

TEST_CASE( "Subtraction" ) {

    SECTION( "operator-(bigint::s<1024>, bigint::s<1024>) against gmp" ) {
        TIMES(large_test_repeats) {
            constexpr size_t TEST_BIT_SZ = 1024;
            constexpr size_t DATA_ARR_SZ = TEST_BIT_SZ / (sizeof(uint64_t) * 8);

            std::array<uint64_t, DATA_ARR_SZ> datain1, datain2;
            init_random_data(datain1); init_random_data(datain2);

            bool sign1 = mt32() & 1, sign2 = mt32() & 1;

            auto bint1 = create_signed<TEST_BIT_SZ>(datain1, sign1);
            auto bint2 = create_signed<TEST_BIT_SZ>(datain2, sign2);

            mpz_t gmpint1, gmpint2;
            create_mpz(&gmpint1, datain1, sign1);
            create_mpz(&gmpint2, datain2, sign2);

            auto bint_result = bint1 - bint2;

            mpz_t gmpint_result;
            mpz_init(gmpint_result);
            mpz_sub(gmpint_result, gmpint1, gmpint2);

            REQUIRE(mpz_to_signed<TEST_BIT_SZ + 1>(&gmpint_result) == bint_result);
        }
    }
}

TEST_CASE( "Multiplication" ) {

    SECTION( "operator*(bigint::s<1024>, bigint::s<1024>) against gmp" ) {

        TIMES(large_test_repeats) {
            constexpr size_t TEST_BIT_SZ = 1024;
            constexpr size_t DATA_ARR_SZ = TEST_BIT_SZ / (sizeof(uint64_t) * 8);

            std::array<uint64_t, DATA_ARR_SZ> datain1, datain2;
            init_random_data(datain1); init_random_data(datain2);

            bool sign1 = mt32() & 1, sign2 = mt32() & 1;
            auto bint1 = create_signed<TEST_BIT_SZ>(datain1, sign1);
            auto bint2 = create_signed<TEST_BIT_SZ>(datain2, sign2);

            mpz_t gmpint1, gmpint2;
            create_mpz(&gmpint1, datain1, sign1);
            create_mpz(&gmpint2, datain2, sign2);

            auto bint_result = bint1 * bint2;

            mpz_t gmpint_result;
            mpz_init(gmpint_result);
            mpz_mul(gmpint_result, gmpint1, gmpint2);

            auto test = (mpz_to_signed<2 * TEST_BIT_SZ>(&gmpint_result) == bint_result);
            if(!test) {

                log(bint_result.binary_string());
                log(mpz_to_signed<2 * TEST_BIT_SZ>(&gmpint_result).binary_string());
                log();
                REQUIRE(test);
            }
        }
    }

    /*
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
            std::cout << "helo" << std::endl;
            // (*bint_result) = (*bint1) * (*bint2);
            std::cout << "helo2" << std::endl;
            std::cout << ((*bint1) * (*bint2)).binary_string() << std::endl;

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
    */
}
