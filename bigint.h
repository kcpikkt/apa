#ifndef BIGINT_H
#define BIGINT_H
#include <algorithm>
#include <type_traits>
#include <assert.h>
#include <array>
#include <stdint.h>
#include <math.h>

namespace{
    typedef uint8_t impl_t;
    constexpr size_t impl_t_byte_sz = sizeof(impl_t);
    constexpr size_t impl_t_bit_sz = sizeof(impl_t)*8;
}

template<size_t _SZ>
struct bigint{

    static_assert(std::is_unsigned<impl_t>::value);

    //template<size_t SZ1, size_t SZ2>
    //friend bigint<std::max(SZ1, SZ2)+1> add_u(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs);
    //template<size_t SZ1, size_t SZ2> 
    //static bigint<std::max(SZ1, SZ2)+1> add_u(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs);

    constexpr static size_t get_segments_count(){
        if((_SZ-1) & ~_SZ) { 
            return _SZ/(sizeof(impl_t)*8); }
        else { 
            size_t ceil_log_2 = (MSB(_SZ) << 1);
            //static_assert( ceil_log_2 == 0);
            return ceil_log_2/(sizeof(impl_t)*8);
        }
    }
    constexpr static size_t MSB(size_t n){
        n |= (n >> 1); n |= (n >> 2);  n |= (n >> 4);
        n |= (n >> 8); n |= (n >> 16); n |= (n >> 32);
        n += 1; 
        if(n == 0){ return ((size_t)0 | ((size_t)1 << 63)); }
        else{       return (n >> 1);     }
    }

    constexpr static size_t segments_count = get_segments_count();
    constexpr static size_t bit_sz = _SZ;
    constexpr static size_t real_bit_sz = segments_count * impl_t_bit_sz;
    constexpr static size_t binary_width = _SZ + 2;
    //constexpr static size_t decimal_width = need compile-time log10;
    //impl_t _segments[segments_count] = {};
    std::array<impl_t, segments_count> _segments = {};

    uint8_t flags = 0;
    enum { 
        NEGATIVE = 1, 
        TRUNCATED = 2
    };

    bigint(){}
    template<typename T> bigint(T val){
        static_assert(std::is_integral<T>::value);
        flags &= ~TRUNCATED;
        if(val < 0) { flags |= NEGATIVE; }
        size_t uval = llabs(val);
        for(size_t i=0; i < segments_count; i++){
            _segments.at(i) = (impl_t)(uval >> sizeof(impl_t)*8*i);
        }
        if(segments_count * sizeof(impl_t) < sizeof(T)){
            for(size_t i= segments_count; i < sizeof(T)/sizeof(impl_t); i++){
                if((impl_t)(uval >> sizeof(impl_t)*8*i) != 0){ flags |= TRUNCATED; break; }
            }
        }
    }
    template<size_t SZ> bigint(const bigint<SZ>& other){
        flags &= ~TRUNCATED;
        flags = other.flags;
        for(size_t i=0; i< segments_count; i++){ _segments[i] = other.get_segment(i); }
        if constexpr (bigint<SZ>::segments_count > segments_count){
            for(size_t i = segments_count; i < bigint<SZ>::segments_count; i++){
                if(other.get_segment(i) != 0){ flags |= TRUNCATED; break; }
            }
        }
    }
    inline impl_t get_segment(size_t index) const {
        return (index < segments_count) ? _segments[index] : 0;
    }
    inline bool bit_at(size_t index) const {
        size_t segment_i = floor(index/impl_t_bit_sz);
        size_t bit_i = index - ( segment_i * impl_t_bit_sz);
        return (get_segment(segment_i) & (1 << bit_i));
    }
    inline bool     is_negative()   const { return  (flags & NEGATIVE);  }
    inline bool     is_positive()   const { return !(flags & NEGATIVE);  }
    inline int8_t   sign()          const { return (flags & NEGATIVE) ? -1 : 1; }
    inline void     set_sign(bool s)      { flags &= ~NEGATIVE; flags |= NEGATIVE * s; }
    inline bool     was_truncated() const { return (flags & TRUNCATED); }
    inline void     toggle_sign()         { flags ^= NEGATIVE; }
    inline bool     is_zero()       const { 
        for(auto e : _segments) if(e != 0) return false;
        return true; 
    }
    size_t clz() const {          // count leading zeros
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
    size_t ctz() const {          // count trailing zeros
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

    operator bool() const { return !is_zero(); }

    void to_string_binary(char * str) const {
        str[0] = is_negative() ? '-' : '+';
        for(size_t i=0; i < _SZ; i++){
            size_t segment_i = floor(i/impl_t_bit_sz);
            size_t bit_i = i - ( segment_i * impl_t_bit_sz);
            str[i+1] = 
            _segments[segments_count - segment_i-1] & (1 << (impl_t_bit_sz - bit_i-1)) 
            ? '1' : '0';
        }
        str[binary_width-1] = '\0';
    }
};

// Arithmetic Operators
// TODO: /, %
// IN PLACE: +, -, *
//

template<size_t SZ1, size_t SZ2> 
static bigint<std::max(SZ1, SZ2)+1> add_u(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    constexpr size_t ret_sz = std::max(SZ1, SZ2)+1;
    bigint<ret_sz> ret;
    impl_t carry = 0;
    impl_t limit = -1;

    for(size_t i=0; i < ret.segments_count; i++){
        //TODO: change truncation side
        //printf("%i %i | %i\n", lhs.get_segment(i), rhs.get_segment(i), carry);
        ret._segments[i] = lhs.get_segment(i) + rhs.get_segment(i);
        if(ret._segments[i] < lhs.get_segment(i)){ 
            ret._segments[i] += carry; carry = 1; } 
        else { 
            if(carry){
                if(ret._segments[i] == limit)  { ret._segments[i] = 0; }
                else                           { ret._segments[i] += carry; carry = 0;}
            }
            else{ carry = 0; }
        }
    }
    return ret;
};

template<size_t SZ1, size_t SZ2> 
static bigint<std::max(SZ1,SZ2)+1> sub_u(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    constexpr size_t ret_sz = std::max(SZ1, SZ2)+1;
    bigint<ret_sz> ret;
    impl_t carry = 0;
    impl_t limit = -1;

    for(size_t i=0; i<ret.segments_count; i++){
        ret._segments[i] = lhs.get_segment(i) - rhs.get_segment(i) - carry;
        if  (  ( ret.get_segment(i) > lhs.get_segment(i)) ||
               ( rhs.get_segment(i) == limit && carry == 1) )
            { carry = 1; }
        else{ carry = 0; }
    }
    if(carry){ 
        for(size_t i=0; i<ret.segments_count; i++){ ret._segments[i] = ~ret.get_segment(i); }
        ret = add_u(ret,bigint<8>(1) );
        ret.toggle_sign();
    }
    return ret;
}
template<size_t SZ, typename T>
inline bigint<std::max(SZ,sizeof(T)*8)+1> operator+(const bigint<SZ>& lhs, T rhs){
    return operator+(lhs, bigint<sizeof(T)*8>(rhs));
}

template<size_t SZ, typename T>
inline bigint<std::max(SZ,sizeof(T)*8)+1> operator+(T lhs, const bigint<SZ>& rhs){
    return operator+(bigint<sizeof(T)*8>(lhs), rhs);
}
template<size_t SZ1, size_t SZ2>
inline bigint<std::max(SZ1,SZ2)+1> operator+(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    if(lhs.sign() == rhs.sign())    { return add_u<SZ1,SZ2>(lhs, rhs);}
    else                            { 
        if(rhs.is_negative())         return sub_u<SZ1,SZ2>(lhs, rhs);
        else                          return sub_u<SZ1,SZ2>(rhs, lhs);
    }
}
template<size_t SZ, typename T>
inline bigint<std::max(SZ,sizeof(T)*8)+1> operator-(T lhs, const bigint<SZ>& rhs){
    return operator-(bigint<sizeof(T)*8>(lhs), rhs);
}
template<size_t SZ, typename T>
inline bigint<std::max(SZ,sizeof(T)*8)+1> operator-(const bigint<SZ>& lhs, T rhs){
    return operator-(lhs, bigint<sizeof(T)*8>(rhs));
}
template<size_t SZ1, size_t SZ2>
inline bigint<std::max(SZ1,SZ2)+1> operator-(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    if(lhs.sign() == rhs.sign())    { return sub_u<SZ1,SZ2>(lhs, rhs);}
    else                            { return add_u<SZ1,SZ2>(lhs, rhs);}
}

template<size_t SZ1, size_t SZ2> 
static bigint<SZ1+SZ2> mul_u(const bigint<SZ1>& lhs, bigint<SZ2> rhs){
    bigint<SZ1+SZ2> ret = lhs; 
    ret.set_sign(0); rhs.set_sign(0);
    if(!rhs.is_zero() && !lhs.is_zero()){
        auto two_factors = rhs.clz();
        rhs = rhs >> two_factors;
        ret = ret << two_factors;
        bigint<SZ1+SZ2> temp;
        for(size_t i=0; i < (rhs.bit_sz - rhs.ctz()); i++){
            if(rhs.bit_at(i)){
                temp = temp + (ret << i);
            }
        }
        ret = temp;
    }else{
        ret = 0;
    }
    return ret;
}
template<size_t SZ, typename T>
inline bigint<SZ+sizeof(T)*8> operator*(T lhs, const bigint<SZ>& rhs){
    return operator*(bigint<sizeof(T)*8>(lhs), rhs);
}

template<size_t SZ, typename T>
inline bigint<SZ+sizeof(T)*8> operator*(const bigint<SZ>& lhs, T rhs){
    return operator*(lhs, bigint<sizeof(T)*8>(rhs));
}
template<size_t SZ1, size_t SZ2>
inline bigint<SZ1+SZ2> operator*(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    bigint<SZ1+SZ2> ret = mul_u(lhs,rhs);
    //FIXME
    //ret.set_sign(lhs.sign() xor rhs.sign());
    return ret;
}



// Uniary Operators
// TODO: ++, --, +, -, 
// IN PLACE:


// Binary Operators
// TODO: ~, |, ^
// IN PLACE: <<, >>, &,
#pragma GCC diagnostic ignored "-Wshift-overflow"
template<size_t SZ>
inline bigint<SZ> operator<<(const bigint<SZ>& lhs, size_t shift){
    bigint<SZ> ret = lhs;
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
template<size_t SZ>
inline bigint<SZ> operator>>(const bigint<SZ>& lhs, size_t shift){
    bigint<SZ> ret = lhs;
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
#pragma GCC diagnostic pop
template<size_t SZ, typename T>
//std::enable_if_t<std::is_integral_v<T>>
inline bigint<std::min(SZ, sizeof(T)*8)> operator&(const bigint<SZ>& lhs, const T rhs){
    return (lhs & bigint<(sizeof(T)*8)>(rhs));
}

template<size_t SZ1, size_t SZ2>
inline bigint<std::min(SZ1, SZ2)> operator&(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    bigint<std::min(SZ1, SZ2)> ret;
    for(size_t i=0; i < std::min(lhs.segments_count, rhs.segments_count); i++)
        ret._segments[i] = lhs.get_segment(i) & rhs.get_segment(i);
    return ret;
}
// Relational Operators
// TODO:
// IN PLACE: ==, !=, <, >, >=, <=
template<size_t SZ1, size_t SZ2>
inline bool comp_u(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){ //is lhs greater
    for(size_t i = std::max(lhs.segments_count, rhs.segments_count); i > 0; i--){
        auto lhs_seg = lhs.get_segment(i-1);
        auto rhs_seg = rhs.get_segment(i-1);
        if((bool)lhs_seg xor (bool)rhs_seg) return (lhs_seg == 0);
        else                                if (lhs_seg > rhs_seg) return true;
    }
    return false;
}
template<size_t SZ1, size_t SZ2>
inline bool operator>(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    if(lhs.sign() == rhs.sign())
        return (lhs.is_positive() ? comp_u(lhs, rhs) : comp_u(rhs, lhs));
    else{
        if(lhs.is_zero() && rhs.is_zero()) return false;
        else return lhs.is_positive();
    }
}
template<size_t SZ1, size_t SZ2>
inline bool operator<(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    if(lhs.sign() == rhs.sign())
        return (lhs.is_positive() ? comp_u(rhs, lhs) : comp_u(lhs, rhs));
    else{
        if(lhs.is_zero() && rhs.is_zero()) return false;
        else return rhs.is_positive();
    }
}

template<size_t SZ1, size_t SZ2>
inline bool operator==(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    if( lhs.sign() != rhs.sign()) 
        return lhs.is_zero() && rhs.is_zero();

    for(size_t i = 0; i < std::max(lhs.segments_count, rhs.segments_count); i++)
        if(lhs.get_segment(i-1) != rhs.get_segment(i-1)) return false;
    return true;
}
template<size_t SZ1, size_t SZ2>
inline bool operator!=(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    return !(lhs == rhs);
}

template<size_t SZ1, size_t SZ2>
inline bool operator<=(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    return (lhs < rhs || lhs == rhs);
}
template<size_t SZ1, size_t SZ2>
inline bool operator>=(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    return (lhs > rhs || lhs == rhs);
}



#endif

