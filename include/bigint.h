#pragma once
//TODO: reduce stl dependancy
#include <algorithm>
#include <functional>
#include <complex>
#include <vector>
#include <iostream>
#include <type_traits>
#include <assert.h>
#include <array>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#ifndef BIGINT_IMPL_TYPE
    #define BIGINT_IMPL_TYPE uint8_t
#else
    static_assert(std::is_integral<BIGINT_IMPL_TYPE>::value);
    static_assert(std::is_unsigned<BIGINT_IMPL_TYPE>::value);
#endif

namespace bigint{
//DEBUG
template<typename ...Ts> void print(Ts... ts){
    (void) std::initializer_list<int>{ (std::cout << ts, 0)...}; }
template<typename ...Ts> void log(Ts... ts){ print(ts..., '\n'); }
//DEBUG
namespace{
    using impl_t = BIGINT_IMPL_TYPE;
    constexpr size_t impl_t_byte_sz = sizeof(impl_t);
    constexpr size_t impl_t_bit_sz = sizeof(impl_t)*8;
}

constexpr static size_t MSB(size_t n){
    n |= (n >> 1); n |= (n >> 2);  n |= (n >> 4);
    n |= (n >> 8); n |= (n >> 16); n |= (n >> 32);
    n += 1; 
    if(n == 0){ return ((size_t)0 | ((size_t)1 << 63)); }
    else{       return (n >> 1);     }
}

template<size_t _SZ>
class Signed{

    static_assert(std::is_unsigned<impl_t>::value);

public:
    enum { 
        NEGATIVE  = 1,
        TRUNCATED = 2
    };

    Signed();

    template<typename T> 
    Signed(T val);

    template<size_t SZ> 
    Signed(const Signed<SZ>& other);

    template<typename T>
    bool import(T* data, size_t count); // TODO: endianness options and stuff

    constexpr static size_t get_segments_count(){
        if((_SZ-1) & ~_SZ) { 
            return _SZ/(sizeof(impl_t)*8); }
        else { 
            size_t ceil_log_2 = (MSB(_SZ) << 1);
            //static_assert( ceil_log_2 == 0);
            return ceil_log_2/(sizeof(impl_t)*8);
        }
    }

    constexpr static size_t segments_count = get_segments_count();
    constexpr static size_t bit_sz = _SZ;
    constexpr static size_t real_bit_sz = segments_count * impl_t_bit_sz;

    //debug
    void print_segs() const {
        for(impl_t s : _segments) std::cout << (uint32_t)s << " ";
        std::puts("\n");
    }

    inline impl_t  get_segment(size_t index) const;
    inline bool bit_at(size_t index) const;

    inline uint8_t  get_flags()     const { return flags; }
    inline bool     is_negative()   const { return  (flags & NEGATIVE);  }
    inline bool     is_positive()   const { return !(flags & NEGATIVE);  }
    inline int8_t   sign()          const { return (flags & NEGATIVE) ? -1 : 1; }
    inline int8_t   sign_bool()     const { return (flags & NEGATIVE) ?  1 : 0; }
    inline void     set_sign_bool(bool b) { flags &= ~NEGATIVE; flags |= NEGATIVE * b; }
    inline void     set_sign(int s) { flags &= ~NEGATIVE; flags |= NEGATIVE * (s < 0); }
    inline bool     trucated() const { return (flags & TRUNCATED); }
    inline void     negate()         { flags ^= NEGATIVE; }
    inline bool     is_zero()       const { 
        for(auto e : _segments) if(e != 0) return false;
        return true; 
    }

    size_t ctz() const; // count trailing zeros
    size_t clz() const; // count leading zeros

    std::string binary_string() const {
        // size_t real_sz = _SZ - clz();
        std::string str; 
        str.resize(_SZ);
        // str.at(0) = is_negative() ? '-' : '+';
        // for(size_t i=0; i < real_sz; i++)
        //     str[real_sz - i - 1] = bit_at(i) + 48;

        for(size_t i=0; i < _SZ; i++)
            str[_SZ - i - 1] = bit_at(i) + 48;
        return str;
    }

    std::string decimal_string() const {
        size_t real_sz = _SZ - clz();
        //str.at(0) = is_negative() ? '-' : '+';
        size_t num = 0;
        for(size_t i=0; i < real_sz; i++)
            num += bit_at(i) * pow(2,i);
        std::string str = std::to_string(num);
        return str;
    }

// Assignment
    template<size_t SZ>
    Signed<_SZ>& operator=(const Signed<SZ>& other);


// Type Conversions
    operator bool() const { return !is_zero(); }

// Arithmetic Operators
    //add unsigned
    template<size_t SZ1, size_t SZ2> 
    friend Signed<std::max(SZ1, SZ2)+1> add_u(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

    //operator+
    template<size_t SZ, typename T>
    friend inline Signed<std::max(SZ,sizeof(T)*8)+1> operator+(const Signed<SZ>& lhs, T rhs);

    template<size_t SZ, typename T>
    friend inline Signed<std::max(SZ,sizeof(T)*8)+1> operator+(T lhs, const Signed<SZ>& rhs);

    template<size_t SZ1, size_t SZ2>
    friend inline Signed<std::max(SZ1,SZ2)+1> operator+(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

    //operator+=
    // template<size_t SZ, typename T>
    // friend inline Signed<std::max(SZ,sizeof(T)*8)+1> operator+=(const Signed<SZ>& lhs, T rhs);

    // template<size_t SZ, typename T>
    // friend inline Signed<std::max(SZ,sizeof(T)*8)+1> operator+=(T lhs, const Signed<SZ>& rhs);

    template<size_t SZ1, size_t SZ2>
    friend inline Signed<SZ1>& operator+=(Signed<SZ1>& lhs, const Signed<SZ2>& rhs);
 


    //add unsigned
    template<size_t SZ1, size_t SZ2> 
    friend Signed<std::max(SZ1,SZ2)+1> sub_u(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

    //operator-
    template<size_t SZ, typename T>
    friend inline Signed<std::max(SZ,sizeof(T)*8)+1> operator-(T lhs, const Signed<SZ>& rhs);

    template<size_t SZ, typename T>
    friend inline Signed<std::max(SZ,sizeof(T)*8)+1> operator-(const Signed<SZ>& lhs, T rhs);

    template<size_t SZ1, size_t SZ2>
    friend inline Signed<std::max(SZ1,SZ2)+1> operator-(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

    //multiply unsigned
    template<size_t SZ1, size_t SZ2> 
    friend Signed<SZ1+SZ2> mul_u(const Signed<SZ1>& lhs, Signed<SZ2> rhs);

    //operator*
    //template<size_t SZ, typename T>
    //friend inline Signed<SZ+sizeof(T)*8> operator*(T lhs, const Signed<SZ>& rhs);
    //
    //template<size_t SZ, typename T>
    //friend inline Signed<SZ+sizeof(T)*8> operator*(const Signed<SZ>& lhs, T rhs);
    //
    template<size_t SZ1, size_t SZ2>
    friend inline Signed<SZ1+SZ2> operator*(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

// Uniary Operators
    //operator~
    template<size_t SZ>
    friend inline Signed<SZ> operator~(Signed<SZ> lhs);

    //operator prefix++;
    template<size_t SZ>
    friend inline Signed<SZ>& operator++(Signed<SZ>& lhs);
    template<size_t SZ>
    friend inline Signed<SZ> operator++(Signed<SZ>& lhs, int);

    //operator prefix--
    template<size_t SZ>
    friend inline Signed<SZ>& operator--(Signed<SZ> lhs);

    //operator postfix--
    template<size_t SZ>
    friend inline Signed<SZ> operator--(Signed<SZ> lhs, int);

// Binary Operators
    //operator<<
    template<size_t SZ>
    friend inline Signed<SZ> operator<<(const Signed<SZ>& lhs, size_t shift);

    //operator<<=
    template<size_t SZ>
    friend inline Signed<SZ>& operator<<=(Signed<SZ>& lhs, size_t shift);

    //operator>>
    template<size_t SZ>
    friend inline Signed<SZ> operator>>(const Signed<SZ>& lhs, size_t shift);

    //operator>>=
    template<size_t SZ>
    friend inline Signed<SZ>& operator>>=(Signed<SZ>& lhs, size_t shift);

    //operator&
    template<size_t SZ, typename T>
    friend inline Signed<std::min(SZ, sizeof(T)*8)> operator&(const Signed<SZ>& lhs, const T rhs);

    template<size_t SZ1, size_t SZ2>
    friend inline Signed<std::min(SZ1, SZ2)> operator&(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

    //operator|
    template<size_t SZ, typename T>
    friend inline Signed<std::max(SZ, sizeof(T)*8)> operator|(const Signed<SZ>& lhs, const T rhs);

    template<size_t SZ1, size_t SZ2>
    friend inline Signed<std::max(SZ1, SZ2)> operator|(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

    //operator^
    template<size_t SZ, typename T>
    friend inline Signed<std::max(SZ, sizeof(T)*8)> operator^(const Signed<SZ>& lhs, const T rhs);

    template<size_t SZ1, size_t SZ2>
    friend inline Signed<std::max(SZ1, SZ2)> operator^(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

// Relational Operators
    // is lhs greater
    template<size_t SZ1, size_t SZ2>
    friend inline int8_t comp_u      (const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

    //operator>
    template<size_t SZ1, size_t SZ2>
    friend inline bool operator>   (const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

    //operator<
    template<size_t SZ1, size_t SZ2>
    friend inline bool operator<   (const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

    //operator==
    template<size_t SZ1, size_t SZ2>
    friend inline bool operator==  (const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

    //operator!=
    template<size_t SZ1, size_t SZ2>
    friend inline bool operator!=  (const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

    //operator<=
    template<size_t SZ1, size_t SZ2>
    friend inline bool operator<=  (const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

    //operator>=
    template<size_t SZ1, size_t SZ2>
    friend inline bool operator>=  (const Signed<SZ1>& lhs, const Signed<SZ2>& rhs);

private:
    std::array<impl_t, segments_count> _segments = {};
    uint8_t flags = 0;
    double multiplication_error_bound;
}; // class Signed

    template<size_t _SZ>
    using s = Signed<_SZ>;
} // namespace bigint

// =================================================================================
namespace bigint{

// Utility
template<size_t _SZ>
inline impl_t Signed<_SZ>::get_segment(size_t index) const {
    return (index < segments_count) ? _segments[index] : 0;
}

template<size_t _SZ>
inline bool Signed<_SZ>::bit_at(size_t index) const {
    size_t segment_i = floor(index/impl_t_bit_sz);
    size_t bit_i = index - ( segment_i * impl_t_bit_sz);
    return (get_segment(segment_i) & (1 << bit_i));
}


template<size_t _SZ>
size_t Signed<_SZ>::ctz() const {          // count trailing zeros
    size_t ret = 0, i = 0; 
    for(; i < segments_count; i++) {
        if(_segments[i] == 0) ret += impl_t_bit_sz;
        else break;
    }
    if(i != segments_count){
        impl_t mask, j = 0;
        do{
            j++;
            mask = (impl_t)pow(2, j)-1;
        } while( !(_segments[i] & mask) );
        ret += j - 1;
    }
    return ret;
}

template<size_t _SZ>
size_t Signed<_SZ>::clz() const {          // count leading zeros
    size_t ret = 0, i = segments_count; 
    for(; i > 0; i--) {
        if(_segments[i-1] == 0) ret += impl_t_bit_sz;
        else break;
    }
    if(i != 0){
        impl_t mask; size_t j = impl_t_bit_sz;
        do{
            j--;
            mask = ~((impl_t)pow(2, j)-1);
        } while( !(_segments[i-1] & mask) );
        ret += impl_t_bit_sz - j-1;
    }
    return ret;
}


// Constructors
//
template<size_t _SZ>
Signed<_SZ>::Signed() {}

template<size_t _SZ>
template<typename T>
Signed<_SZ>::Signed(T val){
    static_assert(std::is_integral<T>::value);
    flags &= ~TRUNCATED;
    if(val < 0) { flags |= NEGATIVE; }

    size_t uval = val;
    if constexpr(!std::is_unsigned<T>::value) { // just to suppress warnings
        uval = llabs(val);
    }
    size_t i = 0;
    for(; i < segments_count && i < sizeof(uval) / sizeof(impl_t); i++) {
        _segments.at(i) = (impl_t)(uval >> i * sizeof(impl_t) * 8);
    }

    for(; i < segments_count; i++) {
        _segments.at(i) = (impl_t)0;
    }

    if(segments_count * sizeof(impl_t) < sizeof(T)){
        for(size_t i= segments_count; i < sizeof(T)/sizeof(impl_t); i++){
            if((impl_t)(uval >> sizeof(impl_t)*8*i) != 0){ flags |= TRUNCATED; break; }
        }
    }
}

template<size_t _SZ>
template<size_t SZ>
Signed<_SZ>::Signed(const Signed<SZ>& other){
    flags = other.get_flags();
    flags &= ~TRUNCATED;
    for(size_t i=0; i< segments_count; i++){ _segments[i] = other.get_segment(i); }
    if constexpr (Signed<SZ>::segments_count > segments_count){
        for(size_t i = segments_count; i < Signed<SZ>::segments_count; i++){
            if(other.get_segment(i) != 0){ flags |= TRUNCATED; break; }
        }
    }
}

template<size_t _SZ>
template<typename T>
bool Signed<_SZ>::import(T * data, size_t count){ 
    // if(sizeof(T) * count > sizeof(impl_t) * segments_count) return false;

    for(auto s : _segments) s = 0;

    if(sizeof(impl_t) > sizeof(T)){
        for(size_t i = 0; i < count; i++) {
            for(size_t j = 0; j < sizeof(impl_t) / sizeof(T); j++) {
                T dataval = data[i * sizeof(impl_t) / sizeof(T) + j];
                _segments[i * sizeof(impl_t) / sizeof(T)] += dataval << j * sizeof(T) * 8;
            }
        }
    } else {
        for(size_t i = 0; i < count; i++) {
            T dataval = data[i];
            for(size_t j = 0; j < sizeof(T) / sizeof(impl_t); j++) {

                if(i * sizeof(T) / sizeof(impl_t) + j > segments_count) break;

                _segments[i * sizeof(T) / sizeof(impl_t) + j] =
                    dataval >> j * sizeof(impl_t) * 8;
            }
        }
    }
    return true;
}

//Assignment
template<size_t _SZ>
template<size_t SZ>
Signed<_SZ>& Signed<_SZ>::operator=(const Signed<SZ>& other){
    flags = other.get_flags();
    flags &= ~TRUNCATED;
    for(size_t i=0; i < segments_count; i++){ _segments[i] = other.get_segment(i); }
    if constexpr (Signed<SZ>::segments_count > segments_count){
        for(size_t i = segments_count; i < Signed<SZ>::segments_count; i++){
            if(other.get_segment(i) != 0){ flags |= TRUNCATED; break; }
        }
    }
    return *this;
}


//subtract unsigned
template<size_t SZ1, size_t SZ2> 
Signed<std::max(SZ1, SZ2)+1> add_u(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    constexpr size_t RET_SZ = std::max(SZ1, SZ2)+1;
    Signed<RET_SZ> ret = lhs;

    std::function<bool(size_t, impl_t)>
    add_to_segment = [&](size_t index, impl_t val) {
                         impl_t temp = ret._segments[index];
                         ret._segments[index] += val;
                         if(ret._segments[index] < temp) {
                             if(index == ret.segments_count - 1) //overflow
                                 return false;
                             else
                                 return add_to_segment(index + 1, 1);
                         }
                         return true;
                     };

    //TODO: signal overflow
    for(size_t i=0; i < rhs.segments_count; i++)
        add_to_segment(i, rhs._segments[i]);

    return ret;
}

//subtract unsigned
template<size_t SZ1, size_t SZ2> 
Signed<std::max(SZ1,SZ2)+1> sub_u(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){

    assert(comp_u(lhs, rhs) > 0);

    constexpr size_t RET_SZ = std::max(SZ1, SZ2)+1;

    Signed<RET_SZ> ret = lhs;

    std::function<bool(size_t, impl_t)>
        sub_from_segment = [&](size_t index, impl_t val) {
                             impl_t temp = ret._segments[index];
                             ret._segments[index] -= val;
                             if(ret._segments[index] > temp) {
                                 if(index == ret.segments_count - 1) //underflow
                                     return false;
                                 else
                                     return sub_from_segment(index + 1, 1);
                             }
                             return true;
                         };

    //TODO: signal underflow (shouldn't happen - lhs >= rhs)
    for(size_t i=0; i < rhs.segments_count; i++)
        sub_from_segment(i, rhs._segments[i]);

    return ret;
}
//     constexpr size_t RET_SZ = std::max(SZ1, SZ2)+1;
//     Signed<RET_SZ> ret;
//     impl_t carry = 0;
//     impl_t limit = -1;

//     for(size_t i=0; i<ret.segments_count; i++){
//         ret._segments[i] = lhs.get_segment(i) - rhs.get_segment(i) - carry;
//         if  (  ( ret.get_segment(i) > lhs.get_segment(i)) ||
//             ( rhs.get_segment(i) == limit && carry == 1) )
//             { carry = 1; }
//         else{ carry = 0; }
//     }
//     if(carry){ 
//         for(size_t i=0; i<ret.segments_count; i++){ ret._segments[i] = ~ret.get_segment(i); }
//         ret = add_u(ret,Signed<8>(1) );
//         ret.negate();
//     }
//     return ret;
// }

//operator+
template<size_t SZ, typename T>
inline Signed<std::max(SZ,sizeof(T)*8)+1> operator+(const Signed<SZ>& lhs, T rhs){
    return operator+(lhs, Signed<sizeof(T)*8>(rhs));
}

template<size_t SZ, typename T>
inline Signed<std::max(SZ,sizeof(T)*8)+1> operator+(T lhs, const Signed<SZ>& rhs){
    return operator+(Signed<sizeof(T)*8>(lhs), rhs);
}

template<size_t SZ1, size_t SZ2>
inline Signed<std::max(SZ1,SZ2)+1> operator+(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    Signed<std::max(SZ1,SZ2)+1> ret;

    if(lhs.sign() != rhs.sign()) {
        if(comp_u(lhs, rhs) > 0) {
            ret = sub_u(lhs,rhs);
            ret.set_sign( lhs.sign());
        } else {
            ret = sub_u(rhs,lhs);
            ret.set_sign(-lhs.sign());
        }
    } else {
        ret = add_u(lhs, rhs);
        ret.set_sign(lhs.sign());
    }

    return ret;
}

//operator+=
//TODO: SFINAE is integral
//template<size_t SZ, typename T>
//inline Signed<std::max(SZ,sizeof(T)*8)+1>& operator+=(const Signed<SZ>& lhs, T rhs){
//    lhs = lhs + rhs;
//    return lhs;
//}

// template<size_t SZ, typename T>
// inline Signed<std::max(SZ,sizeof(T)*8)+1>& operator+=(T lhs, const Signed<SZ>& rhs){

// }

template<size_t SZ1, size_t SZ2>
inline Signed<SZ1>& operator+=(Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    lhs = lhs + rhs;
    return lhs;
}

//operator-
template<size_t SZ, typename T>
inline Signed<std::max(SZ,sizeof(T)*8)+1> operator-(T lhs, const Signed<SZ>& rhs){
    return operator-(Signed<sizeof(T)*8>(lhs), rhs);
}

template<size_t SZ, typename T>
inline Signed<std::max(SZ,sizeof(T)*8)+1> operator-(const Signed<SZ>& lhs, T rhs){
    return operator-(lhs, Signed<sizeof(T)*8>(rhs));
}

template<size_t SZ1, size_t SZ2>
inline Signed<std::max(SZ1,SZ2)+1> operator-(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    Signed<std::max(SZ1,SZ2)+1> ret;

    if(lhs.sign() != rhs.sign()) {
        ret = add_u(lhs, rhs);
        ret.set_sign(lhs.sign());
    } else {
        if(comp_u(lhs, rhs) > 0) {
            ret = sub_u(lhs, rhs);
            ret.set_sign( lhs.sign());
        } else {
            ret = sub_u(rhs, lhs);
            ret.set_sign(-lhs.sign());
        }
    }
    return ret;
}

// FFT
template<typename T>
static bool is_pow2(T v) {
    static_assert(std::is_unsigned<T>::value);
    static_assert(std::is_integral<T>::value);
    return v && !(v & (v - 1));
}

template<typename T>
static constexpr T bit_reverse(T in, uint8_t count = sizeof(T) * 8) {
    static_assert(std::is_unsigned<T>::value);
    static_assert(std::is_integral<T>::value);
    T out = in & 1;

    for(uint8_t i = 0; i < count - 1; i++) {
        in >>= 1;
        out <<= 1;
        out |= in & 1;
    }
    return out;
}
static_assert(bit_reverse(0u) == 0); // constexpr check

template<typename I>
void fast_fourier_transform(I first, I last, bool inverse = false){
    size_t size = last - first;

    assert(is_pow2(size));

    for(size_t i = 0; i < size; i++) {
        size_t temp = bit_reverse(i, std::log2(MSB(size - 1)) + 1);

        if(temp < i)
            std::swap( *(first + temp), *(first + i));
    }

    for(int64_t i = std::log2(size) - 1;  i >= 0; i--) {
        for(uint64_t j = 0; j < (1 << i); j++) {
            auto part_sz = size / (1 << i);
            for(uint64_t k = 0; k < size / (1 << (i + 1)); k++) {
                auto w = std::exp(std::complex<double>
                                    (0, (inverse ? -1 : 1) * 2.0 * M_PI * k / part_sz));

                auto index = j * part_sz + k;

                auto bottom = first[index];
                auto top =    first[index + part_sz/2];

                first[index]             = bottom + w * top;
                first[index + part_sz/2] = bottom - w * top;
            }
        }
    }
}

template<size_t SZ1, size_t SZ2> 
Signed<SZ1+SZ2> mul_u(const Signed<SZ1>& lhs, Signed<SZ2> rhs){

    constexpr size_t RET_SZ = 2 * std::max(SZ1, SZ2);
    constexpr size_t pow2_sz = MSB(Signed<2 * std::max(SZ1, SZ2)>::segments_count);

    Signed<RET_SZ> ret;
    std::array<std::complex<double>, pow2_sz> X, Y, Z;

    for(size_t i = 0; i < pow2_sz; i++) {
        X[i] = lhs.get_segment(i);
        Y[i] = rhs.get_segment(i);
    }

    fast_fourier_transform(X.begin(), X.end(), false);
    fast_fourier_transform(Y.begin(), Y.end(), false);

    for(size_t i=0; i<pow2_sz; i++) Z[i] = X[i] * Y[i];

    fast_fourier_transform(Z.begin(), Z.end(), true);

    Signed<RET_SZ + 64> temp; // during 'reassemly' of the number little bit more space is needed
    for(size_t i = 0; i < pow2_sz; i++) {
        temp = std::llround( Z[i].real() );
        temp <<= i * impl_t_bit_sz;
        temp >>= std::log2(pow2_sz);
        ret += temp;
    }
    return ret;
}

template<size_t SZ1, size_t SZ2>
inline Signed<SZ1+SZ2> operator*(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    Signed<SZ1+SZ2> ret = mul_u(lhs,rhs);
    ret.set_sign_bool(lhs.sign() xor rhs.sign());
    return ret;
}

template<size_t SZ, typename T>
inline Signed<SZ+sizeof(T)*8> operator*(T lhs, const Signed<SZ>& rhs){
    return operator*(Signed<sizeof(T)*8>(lhs), rhs);
}

template<size_t SZ, typename T>
inline Signed<SZ+sizeof(T)*8> operator*(const Signed<SZ>& lhs, T rhs){
    return operator*(lhs, Signed<sizeof(T)*8>(rhs));
}


// template<size_t SZ1, size_t SZ2> 
// Signed<SZ1+SZ2> div_u(const Signed<SZ1>& lhs, Signed<SZ2> rhs){
//    Signed<SZ1+SZ2> ret;
//    auto diff = lhs.get_segments_count() - rhs.get_segments_count();
//    for(size_t i=0; i<diff; i++){

//    }
//    return ret;
// }

// template<size_t SZ1, size_t SZ2>
// inline Signed<SZ1+SZ2> operator/(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
// }



//template<size_t SZ1, size_t SZ2> // take two unsinged
//static Signed<SZ1+SZ2> div_u(const Signed<SZ1>& lhs, Signed<SZ2> rhs){
//    if(lhs == rhs) return Signed<8>::ONE;
//    if(lhs <  rhs) return Signed<8>::ZERO; // use comp_u
//}


// Uniary Operators
//operator~
template<size_t _SZ>
inline Signed<_SZ> operator~(Signed<_SZ> lhs){
    //TODO: truncation?
    for(auto& s : lhs._segments) s = ~s;
    return lhs;
}

//operator prefix++;
template<size_t _SZ>
inline Signed<_SZ>& operator++(Signed<_SZ>& lhs){
    lhs = (lhs+1);
    return lhs;
}

//operator postfix++;
template<size_t _SZ>
inline Signed<_SZ> operator++(Signed<_SZ>& lhs, int){
    Signed<_SZ> ret(lhs);
    ++lhs;
    return ret;
}

//operator prefix--
template<size_t _SZ>
inline Signed<_SZ>& operator--(Signed<_SZ> lhs){
    lhs = (lhs-1);
    return lhs;
}

//operator postfix--
template<size_t _SZ>
inline Signed<_SZ> operator--(Signed<_SZ> lhs, int){
    Signed<_SZ> ret(lhs);
    --lhs;
    return ret;
}

// FIXME: not sure about these
// TODO: SFAINE for different 
////operator uniary+
//template<size_t SZ>
//inline Signed<SZ> operator+(Signed<SZ> lhs){
//    lhs.set_sign_bool(false);
//    return lhs;
//}
////operator uniary-
//template<size_t SZ>
//inline Signed<SZ> operator-(Signed<SZ> lhs){
//    lhs.set_sign_bool(true);
//    return lhs;
//}


// Binary Operators
//#pragma GCC diagnostic ignored "-Wshift-overflow"
//operator<<
template<size_t SZ>
inline Signed<SZ> operator<<(const Signed<SZ>& lhs, size_t shift){
    Signed<SZ> ret = lhs;
    for(size_t i=0; i<lhs.segments_count; i++){
        size_t msseg_d = shift / (sizeof(impl_t)*8);       // most / least significant segment distnace
        size_t lsseg_d = shift / (sizeof(impl_t)*8) + 1;   // given in indices
        size_t msseg_i = i - msseg_d;
        size_t lsseg_i = i - lsseg_d;
        impl_t msseg = lhs.get_segment(msseg_i) << (shift - msseg_d*impl_t_bit_sz);
        impl_t lsseg = lhs.get_segment(lsseg_i) >> (impl_t_bit_sz - (shift - msseg_d*impl_t_bit_sz));
        ret._segments[i] = msseg | lsseg;
    }
    return ret;
}


//operator>>=
template<size_t SZ>
inline Signed<SZ>& operator<<=(Signed<SZ>& lhs, size_t shift){
    lhs = lhs << shift;
    return lhs;
}

//operator>>
template<size_t SZ>
inline Signed<SZ> operator>>(const Signed<SZ>& lhs, size_t shift){
    Signed<SZ> ret = lhs;
    for(size_t i=0; i<lhs.segments_count; i++){
        size_t msseg_d = shift / (sizeof(impl_t)*8) + 1;   // given in indices
        size_t lsseg_d = shift / (sizeof(impl_t)*8);       // most /least significant segment distnace
        size_t msseg_i = i + msseg_d;
        size_t lsseg_i = i + lsseg_d;
        impl_t msseg = lhs.get_segment(msseg_i) << (impl_t_bit_sz - (shift - lsseg_d*impl_t_bit_sz));
        impl_t lsseg = lhs.get_segment(lsseg_i) >> (shift - lsseg_d*impl_t_bit_sz);
        ret._segments[i] = msseg | lsseg;
    }
    return ret;
}

//operator>>=
template<size_t SZ>
inline Signed<SZ>& operator>>=(Signed<SZ>& lhs, size_t shift){
    lhs = lhs >> shift;
    return lhs;
}


//#pragma GCC diagnostic pop
//operator&
template<size_t SZ, typename T>
//std::enable_if_t<std::is_integral_v<T>>
inline Signed<std::min(SZ, sizeof(T)*8)> operator&(const Signed<SZ>& lhs, const T rhs){
    return (lhs & Signed<(sizeof(T)*8)>(rhs));
}
template<size_t SZ1, size_t SZ2>
inline Signed<std::min(SZ1, SZ2)> operator&(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    Signed<std::min(SZ1, SZ2)> ret;
    for(size_t i=0; i < std::min(lhs.segments_count, rhs.segments_count); i++)
        ret._segments[i] = lhs.get_segment(i) & rhs.get_segment(i);
    return ret;
}

//operator|
template<size_t SZ, typename T>
//std::enable_if_t<std::is_integral_v<T>>
inline Signed<std::max(SZ, sizeof(T)*8)> operator|(const Signed<SZ>& lhs, const T rhs){
    return (lhs & Signed<(sizeof(T)*8)>(rhs));
}
template<size_t SZ1, size_t SZ2>
inline Signed<std::max(SZ1, SZ2)> operator|(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    Signed<std::max(SZ1, SZ2)> ret;
    for(size_t i=0; i < std::min(lhs.segments_count, rhs.segments_count); i++)
        ret._segments[i] = lhs.get_segment(i) | rhs.get_segment(i);
    return ret;
}

//operator^
template<size_t SZ, typename T>
//std::enable_if_t<std::is_integral_v<T>>
inline Signed<std::max(SZ, sizeof(T)*8)> operator^(const Signed<SZ>& lhs, const T rhs){
    return (lhs & Signed<(sizeof(T)*8)>(rhs));
}
template<size_t SZ1, size_t SZ2>
inline Signed<std::max(SZ1, SZ2)> operator^(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    Signed<std::max(SZ1, SZ2)> ret;
    for(size_t i=0; i < std::min(lhs.segments_count, rhs.segments_count); i++)
        ret._segments[i] = lhs.get_segment(i) ^ rhs.get_segment(i);
    return ret;
}

// Relational Operators
template<size_t SZ1, size_t SZ2>
inline int8_t comp_u(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    for(size_t i = std::max(lhs.segments_count, rhs.segments_count); i > 0; i--){
        impl_t lhs_seg = lhs.get_segment(i-1);
        impl_t rhs_seg = rhs.get_segment(i-1);

        if(lhs_seg == rhs_seg) continue;
        else {
            if(lhs_seg > rhs_seg) return 1;
            return -1;
        }
    }
    return 0;
}

//operator>
template<size_t SZ1, size_t SZ2>
inline bool operator>(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    if(lhs.sign_bool() xor rhs.sign_bool())
        return lhs.sign() > 0;
    else
        return (lhs.sign() * comp_u(lhs, rhs)) > 0;
}

//operator<
template<size_t SZ1, size_t SZ2>
inline bool operator<(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    if(lhs.sign_bool() xor rhs.sign_bool())
        return lhs.sign() < 0;
    else
        return (lhs.sign() * comp_u(lhs, rhs)) < 0;
}

//operator==
template<size_t SZ1, size_t SZ2>
inline bool operator==(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    if(lhs.sign() == rhs.sign())
        return comp_u(lhs, rhs) == 0;
    else
        return lhs.is_zero() && rhs.is_zero();
}

//operator!=
template<size_t SZ1, size_t SZ2>
inline bool operator!=(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    return !(lhs == rhs);
}

//operator<=
template<size_t SZ1, size_t SZ2>
inline bool operator<=(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    return (lhs < rhs || lhs == rhs);
}

//operator>=
template<size_t SZ1, size_t SZ2>
inline bool operator>=(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    return (lhs > rhs || lhs == rhs);
}
} //namespace bigint

