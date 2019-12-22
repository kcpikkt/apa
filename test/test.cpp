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

// NOTE(kacper): just in case, this is technically non-standard macro
#ifndef __COUNTER__
    #define __COUNTER__ 27u
#endif

// HACK(kacper): compile time counter for 'randomness'
#define HASHLINE(X) \
((__LINE__ >> ((__COUNTER__ + X) % 32)) ^ (__LINE__ * 185791827599696181ULL))

#include "apa.h"

template<typename ...Ts> void print(Ts... ts){
    (void) std::initializer_list<int>{ (std::cout << ts, 0)...}; }
template<typename ...Ts> void log(Ts... ts){ print(ts..., '\n'); }

#define TIMES(N) for(uint32_t i=0; i < N; i++)
[[maybe_unused]] constexpr uint32_t large_test_repeats = 1000;
[[maybe_unused]] constexpr uint32_t medium_test_repeats = 100;
[[maybe_unused]] constexpr uint32_t small_test_repeats = 20;

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define S_TYPE(N) "_signed<" STR(N) ", A>"

// namespace{
//     using impl_t = APA_IMPL_TYPE;
//     [[maybe_unused]] constexpr size_t impl_t_size = sizeof(impl_t);
//     [[maybe_unused]] constexpr size_t impl_t_bits = sizeof(impl_t)*8;
// }

// template<size_t SZ, typename T>
// bool equal(apa::s<SZ> bint, T tint) {
//     size_t i = 0;
//     for(; i < bint.segments_count && i < sizeof(T) / sizeof(impl_t); i++){
//         impl_t mask = static_cast<impl_t> (tint >> i * sizeof(impl_t) * 8);

//         if(! (bint.get_segment(i) == mask)) return false;
//     }
//     for(; i < bint.segments_count; i++){
//         if(! (bint.get_segment(i) == 0)) return false;
//     }
//     return true;
// }

std::random_device rd;
std::mt19937    mt32(rd());
std::mt19937_64 mt64(rd());

template<typename T, size_t ARR_SZ>
void init_random_data(std::array<T, ARR_SZ>& data) {
    for(size_t i=0; i<data.size(); i++) 
        data[i] = static_cast<T>(mt64());
}

template<typename T, size_t ARR_SZ>
void create_mpz(mpz_t * mpz, std::array<T, ARR_SZ>& data, bool sign = 0) {
    mpz_init(*mpz);
    mpz_import(*mpz, data.size(), -1, sizeof(T), 0, 0, data.data());
    if(sign) mpz_neg(*mpz, *mpz);
}

template<size_t SZ, typename A>
apa::_signed<SZ, A> mpz_to_signed(mpz_t * mpz) {
    apa::_signed<SZ, A> ret;
    using byte = unsigned char;
    size_t count = (mpz_sizeinbase(*mpz, 2) / (sizeof(byte) * 8)) << 1;

    byte * data = (byte *) malloc(count);

        size_t rcount = 0;
        mpz_export(data, &rcount, -1, sizeof(byte), 0, 0, *mpz);
        REQUIRE( ret.import(data, rcount) );
        ret.set_sign( mpz_sgn(*mpz) );

    free(data);

    return ret;
}

// Compile time random numbers for use in templates
// (Stolen from Jason Turner's C++ weekly ep 44)
struct PCG {
private:
    constexpr static uint64_t seed() {
        uint64_t shifted = 0;
        for(const auto c : __TIME__) {
            shifted <<= 8;
            shifted |= c;
        }
        return shifted;
    }

    constexpr uint64_t pcg_random(size_t iter) {
        uint64_t old = iter * 6364136223846793006ULL + (seed() | 1);
        uint64_t xorshifted = ((old >> 18u) ^ old) >> 27u;
        uint64_t rot = old >> 59u;
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }

public:
    constexpr uint64_t operator()(uint64_t iter) {
        return pcg_random(iter);
    }
    constexpr uint64_t range(uint64_t min, uint64_t max, uint64_t iter){
        return (pcg_random(iter) % (max - min)) + min;
    }
};

struct Test {
    template<typename Attr>
    void random_signed_test(size_t times);

    template<size_t SZ, typename Attr, typename T, size_t ARR_SZ>
    apa::_signed<SZ, Attr> create_signed
        (std::array<T, ARR_SZ>& data, bool sign = 0);

    template<size_t B, typename A>
    void signed_basic_test();

    template<size_t B1, size_t B2, typename A>
    void signed_addition_test();

    template<size_t B1, size_t B2, typename A>
    void signed_subtraction_test();

    template<size_t B1, size_t B2, typename A>
    void signed_multiplication_test();

    template<size_t B, typename A>
    void signed_bitshift_right_test();

    template<size_t B, typename A>
    void signed_bitshift_left_test();

    template<size_t B, typename A>
    void signed_sign_test();

    template<typename Attr>
    void signed_full_random_test(size_t times);

};

TEST_CASE("_signed<?, NumAttr< uint8_t, double >>")
{
    using Attr = apa::NumAttr<uint8_t, double>;
    Test().signed_full_random_test< Attr >(medium_test_repeats);
}

TEST_CASE("_signed<?, NumAttr< uint16_t, double >>")
{
    using Attr = apa::NumAttr<uint16_t, double>;
    Test().signed_full_random_test< Attr >(medium_test_repeats);
}

// TEST_CASE("_signed<?, NumAttr< uint32_t, double >>")
// {
//     using Attr = apa::NumAttr<uint32_t, long double>;
//     Test().signed_full_random_test< Attr >(medium_test_repeats);
// }

// TEST_CASE("_signed<?, NumAttr< uint64_t, double >>")
// {
//     using Attr = apa::NumAttr<uint64_t, double>;
//     Test().signed_full_random_test< Attr >(medium_test_repeats);
// }




template<typename Attr>
void Test::signed_full_random_test(size_t times) {
    PCG pcg;
    constexpr auto pow2 =
        [](uint64_t p) -> uint64_t { return ((uint64_t)1 << p); };

    SECTION( "Basics" ) {
        TIMES(times) signed_sign_test<
            pcg.range(0, pow2(10), HASHLINE(0)),
            Attr>();
    }
    SECTION( "Addition against gmp" ) {

        TIMES(times) signed_addition_test<
            pcg.range(0, pow2(10), HASHLINE(0)),
            pcg.range(0, pow2(10), HASHLINE(1)),
            Attr>();

        TIMES(times) signed_addition_test<
            pcg.range(pow2(10), pow2(12), HASHLINE(0)),
            pcg.range(pow2(10), pow2(12), HASHLINE(1)),
            Attr>();

        TIMES(times) signed_addition_test<
            pcg.range(pow2(12), pow2(14), HASHLINE(0)),
            pcg.range(pow2(12), pow2(14), HASHLINE(1)),
            Attr>();
    }

    SECTION( "Subtraction against gmp" ) {

        TIMES(times) signed_subtraction_test<
            pcg.range(0, pow2(10), HASHLINE(0)),
            pcg.range(0, pow2(10), HASHLINE(1)),
            Attr>();

        TIMES(times) signed_subtraction_test<
            pcg.range(pow2(10), pow2(12), HASHLINE(0)),
            pcg.range(pow2(10), pow2(12), HASHLINE(1)),
            Attr>();

        TIMES(times) signed_subtraction_test<
            pcg.range(pow2(12), pow2(14), HASHLINE(0)),
            pcg.range(pow2(12), pow2(14), HASHLINE(1)),
            Attr>();
    }


    SECTION( "Multiplication against gmp" ) {

        TIMES(times) signed_multiplication_test<
            pcg.range(0, pow2(10), HASHLINE(0)),
            pcg.range(0, pow2(10), HASHLINE(1)),
            Attr>();

        TIMES(times) signed_multiplication_test<
            pcg.range(pow2(10), pow2(12), HASHLINE(0)),
            pcg.range(pow2(10), pow2(12), HASHLINE(1)),
            Attr>();
    }

    SECTION( "Bitshift right against gmp" ) {

        TIMES(times) signed_bitshift_right_test<
            pcg.range(0, pow2(10), HASHLINE(0)),
            Attr>();

        TIMES(times) signed_bitshift_right_test<
            pcg.range(pow2(10), pow2(12), HASHLINE(0)),
            Attr>();

        TIMES(times) signed_bitshift_right_test<
            pcg.range(pow2(12), pow2(14), HASHLINE(0)),
            Attr>();
    }

    SECTION( "Bitshift left against gmp" ) {

        TIMES(times) signed_bitshift_left_test<
            pcg.range(0, pow2(10), HASHLINE(0)),
            Attr>();

        TIMES(times) signed_bitshift_left_test<
            pcg.range(pow2(10), pow2(12), HASHLINE(0)),
            Attr>();

        TIMES(times) signed_bitshift_left_test<
            pcg.range(pow2(12), pow2(14), HASHLINE(0)),
            Attr>();
    }
}

template<size_t SZ, typename Attr, typename T, size_t ARR_SZ>
apa::_signed<SZ, Attr> Test::create_signed
    (std::array<T, ARR_SZ>& data, bool sign)
{
    apa::_signed<SZ, Attr> ret;
    REQUIRE(ret.import(data.data(), data.size()));
    ret.set_sign_bool(sign);
    return ret;
}

/*
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

*/

template<size_t B1, size_t B2, typename A>
void Test::signed_addition_test()
{
    constexpr size_t DATA_ARR_SZ1 = B1 / (sizeof(uint64_t) * 8);
    constexpr size_t DATA_ARR_SZ2 = B2 / (sizeof(uint64_t) * 8);

    std::array<uint64_t, DATA_ARR_SZ1> datain1;
    std::array<uint64_t, DATA_ARR_SZ2> datain2;
    init_random_data(datain1); init_random_data(datain2);

    bool sign1 = mt32() & 1, sign2 = mt32() & 1;

    auto bint1 = create_signed<B1, A>(datain1, sign1);
    auto bint2 = create_signed<B2, A>(datain2, sign2);

    mpz_t gmpint1, gmpint2;
    create_mpz(&gmpint1, datain1, sign1);
    create_mpz(&gmpint2, datain2, sign2);

    auto bint_result = bint1 + bint2;

    mpz_t gmpint_result;
    mpz_init(gmpint_result);
    mpz_add(gmpint_result, gmpint1, gmpint2);

    REQUIRE(mpz_to_signed<std::max(B1, B2) + 1, A>(&gmpint_result) == bint_result);
}

template<size_t B1, size_t B2, typename A>
void Test::signed_subtraction_test()
{
    constexpr size_t DATA_ARR_SZ1 = B1 / (sizeof(uint64_t) * 8);
    constexpr size_t DATA_ARR_SZ2 = B2 / (sizeof(uint64_t) * 8);

    std::array<uint64_t, DATA_ARR_SZ1> datain1;
    std::array<uint64_t, DATA_ARR_SZ2> datain2;
    init_random_data(datain1); init_random_data(datain2);

    bool sign1 = mt32() & 1, sign2 = mt32() & 1;

    auto bint1 = create_signed<B1, A>(datain1, sign1);
    auto bint2 = create_signed<B2, A>(datain2, sign2);

    mpz_t gmpint1, gmpint2;
    create_mpz(&gmpint1, datain1, sign1);
    create_mpz(&gmpint2, datain2, sign2);

    auto bint_result = bint1 - bint2;

    mpz_t gmpint_result;
    mpz_init(gmpint_result);
    mpz_sub(gmpint_result, gmpint1, gmpint2);

    REQUIRE(mpz_to_signed<std::max(B1, B2) + 1, A>(&gmpint_result) == bint_result);
}

template<size_t B1, size_t B2, typename A>
void Test::signed_multiplication_test()
{
    constexpr size_t DATA_ARR_SZ1 = B1 / (sizeof(uint64_t) * 8);
    constexpr size_t DATA_ARR_SZ2 = B2 / (sizeof(uint64_t) * 8);

    std::array<uint64_t, DATA_ARR_SZ1> datain1;
    std::array<uint64_t, DATA_ARR_SZ2> datain2;
    init_random_data(datain1); init_random_data(datain2);

    bool sign1 = mt32() & 1, sign2 = mt32() & 1;

    std::string str;

    auto bint1 = create_signed<B1, A>(datain1, sign1);
    auto bint2 = create_signed<B2, A>(datain2, sign2);


    mpz_t gmpint1, gmpint2;
    create_mpz(&gmpint1, datain1, sign1);
    create_mpz(&gmpint2, datain2, sign2);

    auto bint_result = bint1 * bint2;

    mpz_t gmpint_result;
    mpz_init(gmpint_result);
    mpz_mul(gmpint_result, gmpint1, gmpint2);

    REQUIRE(mpz_to_signed<B1 + B2, A>(&gmpint_result) == bint_result);
}

template<size_t B, typename A>
void Test::signed_bitshift_right_test()
{
    constexpr size_t DATA_ARR_SZ = B / (sizeof(uint64_t) * 8);

    std::array<uint64_t, DATA_ARR_SZ> datain1, datain2;
    init_random_data(datain1); init_random_data(datain2);

    auto bint1 = create_signed<B, A>(datain1, false);
    auto shift = mt64() % B;

    mpz_t gmpint1;
    create_mpz(&gmpint1, datain1, false);

    auto bint_result = bint1 >> shift;

    mpz_t gmpint_result;
    mpz_init(gmpint_result);
    mpz_fdiv_q_2exp(gmpint_result, gmpint1, shift);

    REQUIRE(mpz_to_signed<B, A>(&gmpint_result) == bint_result);
}

template<size_t B, typename A>
void Test::signed_bitshift_left_test()
{
    constexpr size_t DATA_ARR_SZ = B / (sizeof(uint64_t) * 8);

    std::array<uint64_t, DATA_ARR_SZ> datain1, datain2;
    init_random_data(datain1); init_random_data(datain2);

    auto bint1 = create_signed<B, A>(datain1, false);
    auto shift = mt64() % B;

    mpz_t gmpint1;
    create_mpz(&gmpint1, datain1, false);

    auto bint_result = bint1 << shift;

    mpz_t gmpint_result;
    mpz_init(gmpint_result);
    mpz_mul_2exp(gmpint_result, gmpint1, shift);

    REQUIRE(mpz_to_signed<B, A>(&gmpint_result) == bint_result);
}

template<size_t B, typename A>
void Test::signed_sign_test()
{
    apa::_signed<B, A> bint;
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
