#ifndef BIGINT_H
#define BIGINT_H
#include <algorithm>
#include <type_traits>
#include <assert.h>
#include <stdint.h>
#include <math.h>

namespace{
    typedef uint8_t impl_t;
}

template<size_t _SZ>
struct bigint{

    static_assert(std::is_unsigned<impl_t>::value);
    static_assert(_SZ >= 8);

    constexpr static size_t get_segments_count(){
        if(_SZ & ~_SZ == 0) { 
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
    constexpr static size_t binary_width = _SZ + 2;
    //constexpr static size_t decimal_width = need compile-time log10;
    impl_t _segments[segments_count] = {};

    uint8_t flags = 0;
    enum { 
        NEGATIVE = 1, 
        TRUNCATED = 2
    };

    bigint(){}
    template<typename T> bigint(T val){
        static_assert(std::is_integral<T>::value);
        if(val < 0) { flags |= NEGATIVE; }
        size_t uval = abs(val);
        for(size_t i=0; i < sizeof(T)/sizeof(impl_t); i++){
            _segments[i] = (impl_t)(uval >> sizeof(impl_t)*8*i);
        }
    }
    template<size_t SZ> bigint(const bigint<SZ>& other){
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
    inline bool     is_negative()   const { return (flags & NEGATIVE);  }
    inline int8_t   sign()          const { return (flags & NEGATIVE) ? -1 : 1; }
    inline bool     was_truncated() const { return (flags & TRUNCATED); }
    inline void     toggle_sign()         { flags ^= NEGATIVE; }

    void to_string_binary(char * str) const {
        str[0] = is_negative() ? '-' : '+';
        for(size_t i=0; i < _SZ; i++){
            size_t segment_i = floor(i/(sizeof(impl_t)*8));
            size_t bit_i = i - ( segment_i * sizeof(impl_t) * 8 );
            if (_segments[segment_i] & (1 << bit_i)) { str[i+1] = '1';}
            else                                     { str[i+1] = '0';}
        }
        str[binary_width-1] = '\0';
    }
};

template<size_t SZ1, size_t SZ2> 
static bigint<std::max(SZ1, SZ2)+1> add_u(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    constexpr size_t ret_sz = std::max(SZ1, SZ2)+1;
    bigint<ret_sz> ret;
    uint8_t carry = 0;
    impl_t  limit = -1;

    for(size_t i=0; i < ret.segments_count; i++){
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
    bigint<std::max(SZ1,SZ2)+1> ret;
    uint8_t carry = 0;
    impl_t limit = -1;

    for(int i=0; i<ret.segments_count; i++){
        ret._segments[i] = lhs.get_segment(i) - rhs.get_segment(i) - carry;
        if  (  (rhs.get_segment(i) > lhs.get_segment(i)) ||
               (rhs.get_segment(i) == limit && carry == 1) )
        {
            carry = 1;
        }
    }
    if(carry){ 
        ret.toggle_sign();
        for(int i=0; i<std::max(SZ1, SZ2); i++){ ret._segments[i] = ~ret.get_segment(i); }
        //ret = ret + 1;
    }
    return ret;
}

template<size_t SZ, typename T>
inline bigint<std::max(SZ,sizeof(T)*8)+1> operator+(const bigint<SZ>& lhs, T rhs){
    return operator+(lhs, bigint<sizeof(T)*8>(rhs));
}
template<size_t SZ, typename T>
inline bigint<std::max(SZ,sizeof(T)*8)+1> operator-(const bigint<SZ>& lhs, T rhs){
    return operator-(lhs, bigint<sizeof(T)*8>(rhs));
}
template<size_t SZ, typename T>
inline bigint<std::max(SZ,sizeof(T)*8)+1> operator+(T lhs, const bigint<SZ>& rhs){
    return operator+(bigint<sizeof(T)*8>(lhs), rhs);
}
template<size_t SZ, typename T>
inline bigint<std::max(SZ,sizeof(T)*8)+1> operator-(T lhs, const bigint<SZ>& rhs){
    return operator-(bigint<sizeof(T)*8>(lhs), rhs);
}

template<size_t SZ1, size_t SZ2>
inline bigint<std::max(SZ1,SZ2)+1> operator+(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    if(lhs.sign() == rhs.sign())    { return add_u<SZ1,SZ2>(lhs, rhs);}
    else                            { return sub_u<SZ1,SZ2>(lhs, rhs);}
}
template<size_t SZ1, size_t SZ2>
inline bigint<std::max(SZ1,SZ2)+1> operator-(const bigint<SZ1>& lhs, const bigint<SZ2>& rhs){
    if(lhs.sign() == rhs.sign())    { return sub_u<SZ1,SZ2>(lhs, rhs);}
    else                            { return add_u<SZ1,SZ2>(lhs, rhs);}
}
#endif

