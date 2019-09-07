#pragma once
//TODO: reduce stl dependancy
#include <functional>
#include <complex>
#include <iostream>
#include <type_traits>
#include <assert.h>
#include <stdint.h>
#include <tgmath.h>

#ifndef APA_SEG_TYPE
    #define APA_SEG_TYPE uint16_t
#endif

#ifndef APA_FLOAT_TYPE
    #define APA_FLOAT_TYPE double
#endif

namespace apa {

template<typename F>
struct float_info {
// NOTE(kacper):
//       Is there platform agnostic way of acquiring floating point types
//       exponent and mantissa bit counts for a given float type?
//       Since I don't know how to do that I just assumed IEC 559/IEEE 754
//       standard floating points, hence the assert here - so at least
//       it won't compile in case of non compatibility
private:
    constexpr static size_t _sgn_bits() {
        return 1;
    }

    constexpr static size_t _exp_bits() {
        switch( sizeof(F) ) {
            case 4:  return 8;
            case 8:  return 11;
            case 10: return 15;
            case 16: return 15;
        }
    }
    constexpr static size_t _mnt_bits() {
        switch( sizeof(F) ){
            case 4:  return 23;
            case 8:  return 52;
            case 10: return 64;
            case 16: return 112;
        }
    }
public:
    constexpr static size_t size = sizeof(F) * 8;
    constexpr static size_t bits = sizeof(F) * 8;

    constexpr static size_t sgn_bits = _sgn_bits(); // sign
    constexpr static size_t exp_bits = _exp_bits(); // exponent
    constexpr static size_t mnt_bits = _mnt_bits(); // mantissa

    static_assert(bits == sgn_bits + exp_bits + mnt_bits);

    constexpr static size_t mnt_full_words =
        mnt_bits / (sizeof(unsigned char) * 8);
};
//     }
// }

// }

// constexpr size_t float_t_bits =
//     floating_point_info<float_t>::sgn_bits +
//     float_t_exp_bits() +
//     float_t_mnt_bits();

inline bool is_little_endian() {
    const int x { 0x01 };
    const void * addr = static_cast<const void *>(&x);
    return static_cast<const unsigned char *>(addr);
}
// constexpr size_t word_bits = sizeof(uint8_t) * 8;

// constexpr size_t mantissa_full_words =
//     float_t_mnt_bits() / word_bits;

// constexpr size_t mantissa_residue_bits =
    // float_t_mnt_bits() - mantissa_full_words * word_bits;

template<typename F>
struct float_decomp {
    // little-endian here!
    uint8_t mnt[float_info<F>::mnt_full_words + 1]; 
    int exp, sgn;
};

// union float_t_cast{
//     float_t_cast(float_t _f) : f(_f) {}
//     float_t f;
//     uint8_t b[sizeof(f)/sizeof(uint8_t)];
//     static_assert(sizeof(b) == sizeof(float_t));
// };


template<typename F>
float_decomp<F> float_decompose(F f) {
    static_assert(std::is_floating_point<F>::value);

    using info = float_info<F>;

    auto n_mask = [] (size_t n) { return ((size_t)1 << n) - 1; };
    float_decomp<F> d;

    F norm = frexp(f, &d.exp);

    using uword = unsigned char;
    uword * fbyte = reinterpret_cast<uword *>(&norm);

    // TODO(kacper): use C++20 std::endian
    if( is_little_endian() ){
        size_t i = 0;
        for(; i < info::mnt_full_words; i++)
          d.mnt[i] = fbyte[i];


        d.mnt[i] = fbyte[i] &
            n_mask(info::mnt_bits - info::mnt_full_words * sizeof(uword));
    } else {
        // TODO(kacper): implement for big-endian
        assert(false);
    }
    return d;
}

constexpr static size_t MSB(size_t n){
    n |= (n >> 1); n |= (n >> 2);  n |= (n >> 4);
    n |= (n >> 8); n |= (n >> 16); n |= (n >> 32);
    n += 1; 
    if(n == 0){ return ((size_t)0 | ((size_t)1 << 63)); }
    else{       return (n >> 1);     }
}


template<typename _seg_t, typename _float_t> struct NumAttr {
    using seg_t   = _seg_t;
    using float_t = _float_t;

    static_assert(std::is_integral<seg_t>::value);
    static_assert(std::is_unsigned<seg_t>::value);
    static_assert(std::is_floating_point<float_t>::value);
    static_assert(std::numeric_limits<float_t>::is_iec559);
    // NOTE(kacper): (the same thing as before)
    //       Is there platform agnostic way of acquiring floating point types
    //       exponent and mantissa bit counts for a given float type?
    //       Since I don't know how to do that I just assumed IEC 559/IEEE 754
    //       standard floating points, hence the assert here - so at least
    //       it won't compile in case of non compatibility

    constexpr static size_t seg_t_size   = sizeof(seg_t);
    constexpr static size_t seg_t_bits   = sizeof(seg_t) * 8;

    constexpr static size_t float_t_size = sizeof(float_t);
    constexpr static size_t float_t_bits = sizeof(float_t) * 8;

    using float_t_info = float_info<float_t>;
};

using DefaultNumAttr = NumAttr<APA_SEG_TYPE, APA_FLOAT_TYPE>;

template<size_t _SZ, typename _A = DefaultNumAttr>
class _signed {

    using seg_t    = typename _A::seg_t;
    using float_t  = typename _A::float_t;

    constexpr static size_t seg_t_size   = _A::seg_t_size;
    constexpr static size_t seg_t_bits   = _A::seg_t_bits;

    constexpr static size_t float_t_size = _A::float_t_size;
    constexpr static size_t float_t_bits = _A::float_t_bits;

    using float_t_info = typename _A::float_t_info;

public:
    enum { 
        NEGATIVE  = 1,
        TRUNCATED = 2
    };

    _signed();

    // Constructor for integers
    template<typename I,
            typename =
            typename std::enable_if<std::is_integral<I>::value>::type>
    _signed(I val);

    // Constructor for floats
    template<typename F,
            typename = 
            typename std::enable_if<std::is_floating_point<F>::value>::type,
            typename = void>
    //      ^ HACK(kacper): so I can have unambigous definition outside class
    _signed(F val);


    template<size_t SZ> 
    _signed(const _signed<SZ, _A>& other);

    template<typename T>
    bool import(T* data, size_t count); // TODO: endianness options and stuff

    constexpr static size_t get_segments_count(){
        if((_SZ-1) & ~_SZ) { 
            return _SZ/(sizeof(seg_t)*8); }
        else {
            size_t ceil_log_2 = (MSB(_SZ) << 1);
            //static_assert( ceil_log_2 == 0);
            return ceil_log_2/(sizeof(seg_t)*8);
        }
    }

    constexpr static size_t segments_count = get_segments_count();
    constexpr static size_t bit_sz = _SZ;
    constexpr static size_t real_bit_sz = segments_count * seg_t_bits;

    //debug
    void print_segs() const {
        for(seg_t s : _segments) std::cout << (uint32_t)s << " ";
        std::puts("\n");
    }

    inline seg_t  get_segment(size_t index) const;
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
    _signed<_SZ, _A>& operator=
    (const _signed<SZ, _A>& other);


// Type Conversions
    explicit operator bool() const { return !is_zero(); }

// Arithmetic Operators
    //add unsigned
    template<size_t SZ1, size_t SZ2, typename A> 
    friend _signed<std::max(SZ1, SZ2)+1, A> add_u
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    //operator+
    template<size_t SZ, typename A, typename T, typename /* SFINAE */>
    friend inline _signed<std::max(SZ,sizeof(T)*8)+1> operator+
      (const _signed<SZ, A>& lhs, T rhs);

    template<size_t SZ, typename A, typename T, typename /* SFINAE */>
    friend inline _signed<std::max(SZ,sizeof(T)*8)+1, A> operator+
      (T lhs, const _signed<SZ, A>& rhs);

    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<std::max(SZ1,SZ2)+1, A> operator+
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    //operator+=
    // template<size_t SZ, typename T>
    // friend inline _signed<std::max(SZ,sizeof(T)*8)+1> operator+=
    // (const _signed<SZ>& lhs, T rhs);

    // template<size_t SZ, typename T>
    // friend inline _signed<std::max(SZ,sizeof(T)*8)+1> operator+=
    // (T lhs, const _signed<SZ>& rhs);

    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<SZ1, A>& operator+=
      (_signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    //add unsigned
    template<size_t SZ1, size_t SZ2, typename A>
    friend _signed<std::max(SZ1,SZ2)+1, A> sub_u
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    //operator-
    template<size_t SZ, typename A, typename T>
    friend inline _signed<std::max(SZ,sizeof(T)*8)+1, A> operator-
      (T lhs, const _signed<SZ, A>& rhs);

    template<size_t SZ, typename A, typename T>
    friend inline _signed<std::max(SZ,sizeof(T)*8)+1, A> operator-
      (const _signed<SZ, A>& lhs, T rhs);

    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<std::max(SZ1,SZ2)+1, A> operator-
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    //multiply unsigned
    template<size_t SZ1, size_t SZ2, typename A>
    friend _signed<SZ1+SZ2, A> mul_u
      (const _signed<SZ1, A>& lhs, _signed<SZ2, A> rhs);

    //operator*
    //template<size_t SZ, typename T>
    //friend inline _signed<SZ+sizeof(T)*8> operator*
    // (T lhs, const _signed<SZ>& rhs);
    //
    //template<size_t SZ, typename T>
    //friend inline _signed<SZ+sizeof(T)*8> operator*
    //(const _signed<SZ>& lhs, T rhs);
    //
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<SZ1+SZ2, A> operator*
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

// Uniary Operators
    //operator~
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A> operator~
      (_signed<SZ, A> lhs);

    //operator prefix++;
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A>& operator++
      (_signed<SZ, A>& lhs);

    template<size_t SZ, typename A>
    friend inline _signed<SZ, A> operator++
      (_signed<SZ, A>& lhs, int);

    //operator prefix--
    template<size_t SZ, typename A>
    friend inline _signed<SZ>& operator--
      (_signed<SZ> lhs);

    //operator postfix--
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A> operator--
      (_signed<SZ, A> lhs, int);

// Binary Operators
    //operator<<
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A> operator<<
      (const _signed<SZ, A>& lhs, size_t shift);

    //operator<<=
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A>& operator<<=
      (_signed<SZ, A>& lhs, size_t shift);

    //operator>>
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A> operator>>
      (const _signed<SZ, A>& lhs, size_t shift);

    //operator>>=
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A>& operator>>=
      (_signed<SZ, A>& lhs, size_t shift);

    //operator&
    template<size_t SZ, typename A, typename T>
    friend inline _signed<std::min(SZ, sizeof(T)*8), A> operator&
      (const _signed<SZ, A>& lhs, const T rhs);

    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<std::min(SZ1, SZ2), A> operator&
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    //operator|
    template<size_t SZ, typename A, typename T>
    friend inline _signed<std::max(SZ, sizeof(T)*8), A> operator|
      (const _signed<SZ, A>& lhs, const T rhs);

    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<std::max(SZ1, SZ2), A> operator|
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    //operator^
    template<size_t SZ, typename A, typename T>
    friend inline _signed<std::max(SZ, sizeof(T)*8), A> operator^
      (const _signed<SZ, A>& lhs, const T rhs);

    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<std::max(SZ1, SZ2), A> operator^
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

// Relational Operators
    // is lhs greater
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline int8_t comp_u
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    //operator>
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline bool operator>
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    //operator<
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline bool operator<
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    //operator==
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline bool operator==
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    //operator!=
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline bool operator!=
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    //operator<=
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline bool operator<=
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    //operator>=
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline bool operator>=
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

private:
    std::array<seg_t, segments_count> _segments = {};
    uint8_t flags = 0;
    double multiplication_error_bound;
}; // class _signed<size_t, typename>

  template<size_t _SZ>
  using s = _signed<_SZ>;

// ================================================= IMPLEMENTATION

// ======================================= UTILITY METHODS
template<size_t _SZ, typename A>
inline typename _signed<_SZ, A>::seg_t
_signed<_SZ, A>::get_segment(size_t index) const
{
    return (index < segments_count) ? _segments[index] : 0;
}

template<size_t _SZ, typename A>
inline bool _signed<_SZ, A>::bit_at(size_t index) const
{
    size_t segment_i = floor(index/seg_t_bits);
    size_t bit_i = index - ( segment_i * seg_t_bits);
    return (get_segment(segment_i) & (1 << bit_i));
}

// count trailing zeros
template<size_t _SZ, typename A>
size_t _signed<_SZ, A>::ctz() const
{
    size_t ret = 0, i = 0; 
    for(; i < segments_count; i++) {
        if(_segments[i] == 0) ret += seg_t_bits;
        else break;
    }
    if(i != segments_count){
        seg_t mask, j = 0;
        do{
            j++;
            mask = (seg_t)pow(2, j)-1;
        } while( !(_segments[i] & mask) );
        ret += j - 1;
    }
    return ret;
}

// count leading zeros
template<size_t _SZ, typename A>
size_t _signed<_SZ, A>::clz() const
{
    size_t ret = 0, i = segments_count; 
    for(; i > 0; i--) {
        if(_segments[i-1] == 0) ret += seg_t_bits;
        else break;
    }
    if(i != 0){
        seg_t mask; size_t j = seg_t_bits;
        do{
            j--;
            mask = ~((seg_t)pow(2, j)-1);
        } while( !(_segments[i-1] & mask) );
        ret += seg_t_bits - j-1;
    }
    return ret;
}


// ======================================= CONSTRUCTORS

template<size_t _SZ, typename _A>
_signed<_SZ, _A>::_signed()
{}


template<size_t _SZ, typename _A>
template<typename I, typename /* SFINAE */>
_signed<_SZ, _A>::_signed(I val)
{
    flags &= ~TRUNCATED;
    if(val < 0) { flags |= NEGATIVE; }

    size_t uval = val;
    if constexpr(!std::is_unsigned<I>::value) { // just to suppress warnings
        uval = llabs(val);
    }
    size_t i = 0;
    for(; i < segments_count && i < sizeof(uval) / sizeof(seg_t); i++) {
        _segments.at(i) = (seg_t)(uval >> i * sizeof(seg_t) * 8);
    }

    for(; i < segments_count; i++) {
        _segments.at(i) = (seg_t)0;
    }

    if(segments_count * sizeof(seg_t) < sizeof(I)){
        for(size_t i= segments_count; i < sizeof(I)/sizeof(seg_t); i++){
            if((seg_t)(uval >> sizeof(seg_t)*8*i) != 0){ flags |= TRUNCATED; break; }
        }
    }
}


template<size_t _SZ, typename _A>
template<typename F, typename /* SFINAE */, typename /* HACK */>
_signed<_SZ, _A>::_signed([[maybe_unused]]F val)
{
    using info = float_info<F>;
    float_decomp<F> d = float_decompose<F>(val);
    import(d.mnt, info::mnt_full_words + 1);
    for(size_t i=0; i<(info::mnt_full_words + 1); i++)
        print(std::bitset<8>(d.mnt[i]), " ");

    log();
    for(size_t i=0; i<4; i++)
        print(std::bitset<_A::seg_t_bits>(_segments[i]), " ");
}

template<size_t _SZ, typename _A>
template<size_t SZ>
_signed<_SZ, _A>::_signed(const _signed<SZ, _A>& other)
{
    flags = other.get_flags();
    flags &= ~TRUNCATED;
    for(size_t i=0; i< segments_count; i++){ _segments[i] = other.get_segment(i); }
    if constexpr (_signed<SZ, _A>::segments_count > segments_count){
            for(size_t i = segments_count; i < _signed<SZ, _A>::segments_count; i++){
            if(other.get_segment(i) != 0){ flags |= TRUNCATED; break; }
        }
    }
}

// import raw binary little-endian data
template<size_t _SZ, typename _A>
template<typename T>
bool _signed<_SZ, _A>::import
(T * data, size_t count)
{
    // if(sizeof(T) * count > sizeof(seg_t) * segments_count) return false;

    for(auto s : _segments) s = 0;

    if(sizeof(seg_t) > sizeof(T)){
        for(size_t i = 0; i < count; i++) {
            for(size_t j = 0; j < sizeof(seg_t) / sizeof(T); j++) {
                T dataval = data[i * sizeof(seg_t) / sizeof(T) + j];
                _segments[i * sizeof(seg_t) / sizeof(T)] += dataval << j * sizeof(T) * 8;
            }
        }
    } else {
        for(size_t i = 0; i < count; i++) {
            T dataval = data[i];
            for(size_t j = 0; j < sizeof(T) / sizeof(seg_t); j++) {

                if(i * sizeof(T) / sizeof(seg_t) + j > segments_count) break;

                _segments[i * sizeof(T) / sizeof(seg_t) + j] =
                    dataval >> j * sizeof(seg_t) * 8;
            }
        }
    }
    return true;
}

//Assignment
template<size_t _SZ, typename _A>
template<size_t SZ>
_signed<_SZ, _A>& _signed<_SZ, _A>::operator=
(const _signed<SZ, _A>& other)
{
    flags = other.get_flags();
    flags &= ~TRUNCATED;
    for(size_t i=0; i < segments_count; i++){ _segments[i] = other.get_segment(i); }
    if constexpr (_signed<SZ, _A>::segments_count > segments_count){
        for(size_t i = segments_count; i < _signed<SZ, _A>::segments_count; i++){
            if(other.get_segment(i) != 0){ flags |= TRUNCATED; break; }
        }
    }
    return *this;
}


//subtract unsigned
template<size_t SZ1, size_t SZ2, typename A>
_signed<std::max(SZ1, SZ2)+1, A> add_u
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    using seg_t = typename A::seg_t;

    constexpr size_t RET_SZ = std::max(SZ1, SZ2)+1;
    _signed<RET_SZ, A> ret = lhs;

    std::function<bool(size_t, seg_t)>
    add_to_segment = [&](size_t index, seg_t val) {
                        seg_t temp = ret._segments[index];
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
template<size_t SZ1, size_t SZ2, typename A>
_signed<std::max(SZ1,SZ2)+1, A> sub_u
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    assert(comp_u(lhs, rhs) > 0);

    using seg_t = typename A::seg_t;

    constexpr size_t RET_SZ = std::max(SZ1, SZ2)+1;

    _signed<RET_SZ, A> ret = lhs;

    std::function<bool(size_t, seg_t)>
        sub_from_segment = [&](size_t index, seg_t val) {
                            seg_t temp = ret._segments[index];
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
//     _signed<RET_SZ> ret;
//     seg_t carry = 0;
//     seg_t limit = -1;

//     for(size_t i=0; i<ret.segments_count; i++){
//         ret._segments[i] = lhs.get_segment(i) - rhs.get_segment(i) - carry;
//         if  (  ( ret.get_segment(i) > lhs.get_segment(i)) ||
//             ( rhs.get_segment(i) == limit && carry == 1) )
//             { carry = 1; }
//         else{ carry = 0; }
//     }
//     if(carry){ 
//         for(size_t i=0; i<ret.segments_count; i++){ ret._segments[i] = ~ret.get_segment(i); }
//         ret = add_u(ret,_signed<8>(1) );
//         ret.negate();
//     }
//     return ret;
// }

//operator+
template<size_t SZ, typename A,typename T,
  typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
inline _signed<std::max(SZ,sizeof(T)*8)+1, A> operator+
(const _signed<SZ, A>& lhs, T rhs)
{
    return operator+(lhs, _signed<sizeof(T)*8, A>(rhs));
}

template<size_t SZ, typename A, typename T,
  typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
inline _signed<std::max(SZ,sizeof(T)*8)+1, A> operator+
(T lhs, const _signed<SZ, A>& rhs)
{
    return operator+(_signed<sizeof(T)*8, A>(lhs), rhs);
}

template<size_t SZ1, size_t SZ2, typename A>
inline _signed<std::max(SZ1,SZ2)+1, A> operator+
  (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    _signed<std::max(SZ1,SZ2)+1, A> ret;

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
//inline _signed<std::max(SZ,sizeof(T)*8)+1>& operator+=(const _signed<SZ>& lhs, T rhs){
//    lhs = lhs + rhs;
//    return lhs;
//}

// template<size_t SZ, typename T>
// inline _signed<std::max(SZ,sizeof(T)*8)+1>& operator+=(T lhs, const _signed<SZ>& rhs){

// }

template<size_t SZ1, size_t SZ2, typename A>
inline _signed<SZ1, A>& operator+=
(_signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs){
    lhs = lhs + rhs;
    return lhs;
}

//operator-
template<size_t SZ, typename A, typename T>
inline _signed<std::max(SZ, sizeof(T) * 8) + 1, A> operator-
(T lhs, const _signed<SZ, A>& rhs)
{
    return operator-(_signed<sizeof(T)*8, A>(lhs), rhs);
}

template<size_t SZ, typename A, typename T>
inline _signed<std::max(SZ, sizeof(T) * 8) + 1, A> operator-
(const _signed<SZ, A>& lhs, T rhs){
    return operator-(lhs, _signed<sizeof(T)*8, A>(rhs));
}

template<size_t SZ1, size_t SZ2, typename A>
inline _signed<std::max(SZ1,SZ2)+1, A> operator-
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    _signed<std::max(SZ1,SZ2)+1, A> ret;

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
                auto w = std::exp(std::complex<float_t>
                                    (0, (inverse ? 1 : -1) * 2.0 * M_PI * k / part_sz));

                auto index = j * part_sz + k;

                auto bottom = first[index];
                auto top    = first[index + part_sz/2];

                first[index]             = bottom + w * top;
                first[index + part_sz/2] = bottom - w * top;
            }
        }
    }
}

template<size_t SZ1, size_t SZ2, typename A>
_signed<SZ1+SZ2, A> mul_u
(const _signed<SZ1, A>& lhs, _signed<SZ2, A> rhs)
{
    constexpr size_t seg_t_bits = A::seg_t_bits;

    constexpr size_t RET_SZ = SZ1 + SZ2;
    // constexpr size_t pow2_sz = MSB(_signed<2 * std::max(SZ1, SZ2)>::segments_count);
    constexpr size_t pow2_sz = MSB(_signed<RET_SZ, A>::segments_count);

    _signed<RET_SZ, A> ret;
    std::array<std::complex<float_t>, pow2_sz> X, Y, Z;

    for(size_t i = 0; i < pow2_sz; i++) {
        X[i] = lhs.get_segment(i);
        Y[i] = rhs.get_segment(i);
    }

    fast_fourier_transform(X.begin(), X.end(), false);
    fast_fourier_transform(Y.begin(), Y.end(), false);

    for(size_t i=0; i<pow2_sz; i++) Z[i] = X[i] * Y[i];

    fast_fourier_transform(Z.begin(), Z.end(), true);

    // during 'reassemly' of the number little bit more space is needed
    _signed<RET_SZ + 512, A> temp = 0; 

    //TODO: double to apa integer
    for(size_t i = 0; i < pow2_sz; i++) {
        temp = std::llround( Z[i].real() );
        temp <<= i * seg_t_bits;
        temp >>= std::log2(pow2_sz);
        ret += temp;
    }
    return ret;
}

template<size_t SZ1, size_t SZ2, typename A>
inline _signed<SZ1+SZ2, A> operator*
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    _signed<SZ1+SZ2> ret = mul_u(lhs,rhs);
    ret.set_sign_bool(lhs.sign() xor rhs.sign());
    return ret;
}

template<size_t SZ, typename A, typename T>
inline _signed<SZ+sizeof(T)*8, A> operator*
(T lhs, const _signed<SZ, A>& rhs)
{
    return operator*(_signed<sizeof(T)*8, A>(lhs), rhs);
}

template<size_t SZ, typename A, typename T>
inline _signed<SZ+sizeof(T)*8, A> operator*
(const _signed<SZ, A>& lhs, T rhs)
{
    return operator*(lhs, _signed<sizeof(T)*8, A>(rhs));
}


// template<size_t SZ1, size_t SZ2> 
// _signed<SZ1+SZ2> div_u(const _signed<SZ1, A>& lhs, _signed<SZ2, A> rhs){
//    _signed<SZ1+SZ2> ret;
//    auto diff = lhs.get_segments_count() - rhs.get_segments_count();
//    for(size_t i=0; i<diff; i++){

//    }
//    return ret;
// }

// template<size_t SZ1, size_t SZ2>
// inline _signed<SZ1+SZ2> operator/(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs){
// }



//template<size_t SZ1, size_t SZ2> // take two unsinged
//static _signed<SZ1+SZ2> div_u(const _signed<SZ1, A>& lhs, _signed<SZ2, A> rhs){
//    if(lhs == rhs) return _signed<8>::ONE;
//    if(lhs <  rhs) return _signed<8>::ZERO; // use comp_u
//}


// Uniary Operators
//operator~
template<size_t _SZ, typename A>
inline _signed<_SZ, A> operator~
(_signed<_SZ, A> lhs)
{
    //TODO: truncation?
    for(auto& s : lhs._segments) s = ~s;
    return lhs;
}

//operator prefix++;
template<size_t _SZ, typename A>
inline _signed<_SZ, A>& operator++
(_signed<_SZ, A>& lhs)
{
    lhs = (lhs+1);
    return lhs;
}

//operator postfix++;
template<size_t _SZ, typename A>
inline _signed<_SZ, A> operator++
(_signed<_SZ, A>& lhs, int)
{
    _signed<_SZ, A> ret(lhs);
    ++lhs;
    return ret;
}

//operator prefix--
template<size_t _SZ, typename A>
inline _signed<_SZ, A>& operator--
(_signed<_SZ, A> lhs)
{
    lhs = (lhs-1);
    return lhs;
}

//operator postfix--
template<size_t _SZ, typename A>
inline _signed<_SZ, A> operator--
(_signed<_SZ, A> lhs, int)
{
    _signed<_SZ, A> ret(lhs);
    --lhs;
    return ret;
}

// FIXME: not sure about these
// TODO: SFAINE for different 
////operator uniary+
//template<size_t SZ>
//inline _signed<SZ> operator+(_signed<SZ> lhs){
//    lhs.set_sign_bool(false);
//    return lhs;
//}
////operator uniary-
//template<size_t SZ>
//inline _signed<SZ> operator-(_signed<SZ> lhs){
//    lhs.set_sign_bool(true);
//    return lhs;
//}


// Binary Operators
//operator<<
// TODO: reduce redundant copy
template<size_t SZ, typename A>
inline _signed<SZ, A> operator<<
(const _signed<SZ, A>& lhs, size_t shift)
{
    using seg_t = typename A::seg_t;
    constexpr size_t seg_t_bits = A::seg_t_bits;

    _signed<SZ, A> ret = lhs;
    if(shift == 0) return ret;

    size_t seg_dist = shift / seg_t_bits;
    size_t seghi_shift = shift - seg_dist * seg_t_bits;
    size_t seglo_shift = seg_t_bits - seghi_shift;

    for(size_t i=0; i < lhs.segments_count; i++){
        // NOTE: this conditional is here due to bitshift operator restriction
        //       to range [0, sizeof(T)*8) on Intel cpu's
        seg_t seglo = seglo_shift < seg_t_bits
            ? lhs.get_segment(i - seg_dist - 1) >> seglo_shift
            : 0;

        seg_t seghi = lhs.get_segment(i - seg_dist) << seghi_shift;

        ret._segments[i] = seglo | seghi;
    }
    return ret;
}


//operator<<=
template<size_t SZ, typename A>
inline _signed<SZ, A>& operator<<=
(_signed<SZ, A>& lhs, size_t shift)
{
    lhs = lhs << shift;
    return lhs;
}

//operator>>
// TODO: reduce redundant copy
template<size_t SZ, typename A>
inline _signed<SZ, A> operator>>
(const _signed<SZ, A>& lhs, size_t shift)
{
    using seg_t = typename A::seg_t;
    constexpr size_t seg_t_bits = A::seg_t_bits;

    _signed<SZ, A> ret = lhs;
    if(shift == 0) return ret;

    size_t seg_dist = shift / seg_t_bits;
    size_t seglo_shift = shift - seg_dist * seg_t_bits;
    size_t seghi_shift = seg_t_bits - seglo_shift;

    for(size_t i=0; i<lhs.segments_count; i++){
        seg_t seglo = lhs.get_segment(i + seg_dist) >> seglo_shift;

        // NOTE: this conditional is here due to bitshift operator restriction
        //       to range [0, sizeof(T)*8) on Intel cpu's
        seg_t seghi = seghi_shift < seg_t_bits
            ? lhs.get_segment(i + seg_dist + 1) << seghi_shift
            : 0;
        ret._segments[i] = seglo | seghi;
    }
    return ret;
}

//operator>>=
template<size_t SZ, typename A>
inline _signed<SZ, A>& operator>>=
(_signed<SZ, A>& lhs, size_t shift)
{
    lhs = lhs >> shift;
    return lhs;
}

//operator&
template<size_t SZ, typename A, typename T>
inline _signed<std::min(SZ, sizeof(T)*8), A> operator&
(const _signed<SZ, A>& lhs, const T rhs)
{
    return (lhs & _signed<(sizeof(T)*8), A>(rhs));
}

template<size_t SZ1, size_t SZ2, typename A>
inline _signed<std::min(SZ1, SZ2), A> operator&
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    _signed<std::min(SZ1, SZ2), A> ret;
    for(size_t i=0; i < std::min(lhs.segments_count, rhs.segments_count); i++)
        ret._segments[i] = lhs.get_segment(i) & rhs.get_segment(i);
    return ret;
}

//operator|
template<size_t SZ, typename T, typename A>
inline _signed<std::max(SZ, sizeof(T)*8), A> operator|
(const _signed<SZ, A>& lhs, const T rhs)
{
    return (lhs & _signed<(sizeof(T)*8)>(rhs));
}

template<size_t SZ1, size_t SZ2, typename A>
inline _signed<std::max(SZ1, SZ2)> operator|
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    _signed<std::max(SZ1, SZ2), A> ret;
    for(size_t i=0; i < std::min(lhs.segments_count, rhs.segments_count); i++)
        ret._segments[i] = lhs.get_segment(i) | rhs.get_segment(i);
    return ret;
}

//operator^
template<size_t SZ, typename T, typename A>
inline _signed<std::max(SZ, sizeof(T)*8)> operator^
(const _signed<SZ, A>& lhs, const T rhs)
{
    return (lhs & _signed<(sizeof(T)*8), A>(rhs));
}

template<size_t SZ1, size_t SZ2, typename A>
inline _signed<std::max(SZ1, SZ2)> operator^
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    _signed<std::max(SZ1, SZ2)> ret;
    for(size_t i=0; i < std::min(lhs.segments_count, rhs.segments_count); i++)
        ret._segments[i] = lhs.get_segment(i) ^ rhs.get_segment(i);
    return ret;
}

// Relational Operators
template<size_t SZ1, size_t SZ2, typename A>
inline int8_t comp_u
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    using seg_t = typename A::seg_t;

    for(size_t i = std::max(lhs.segments_count, rhs.segments_count); i > 0; i--){
        seg_t lhs_seg = lhs.get_segment(i-1);
        seg_t rhs_seg = rhs.get_segment(i-1);

        if(lhs_seg == rhs_seg) continue;
        else {
            if(lhs_seg > rhs_seg) return 1;
            return -1;
        }
    }
    return 0;
}

//operator>
template<size_t SZ1, size_t SZ2, typename A>
inline bool operator>
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    if(lhs.sign_bool() xor rhs.sign_bool())
        return lhs.sign() > 0;
    else
        return (lhs.sign() * comp_u(lhs, rhs)) > 0;
}

//operator<
template<size_t SZ1, size_t SZ2, typename A>
inline bool operator<
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    if(lhs.sign_bool() xor rhs.sign_bool())
        return lhs.sign() < 0;
    else
        return (lhs.sign() * comp_u(lhs, rhs)) < 0;
}

//operator==
template<size_t SZ1, size_t SZ2, typename A>
inline bool operator==
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    if(lhs.sign() == rhs.sign())
        return comp_u(lhs, rhs) == 0;
    else
        return lhs.is_zero() && rhs.is_zero();
}

//operator!=
template<size_t SZ1, size_t SZ2, typename A>
inline bool operator!=
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    return !(lhs == rhs);
}

//operator<=
template<size_t SZ1, size_t SZ2, typename A>
inline bool operator<=
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    return (lhs < rhs || lhs == rhs);
}

//operator>=
template<size_t SZ1, size_t SZ2, typename A>
inline bool operator>=
(const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    return (lhs > rhs || lhs == rhs);
}

} //namespace apa

