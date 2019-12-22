#pragma once
#include <functional>
#include <complex>
#include <iostream>
#include <type_traits>
#include <cassert>
#include <cstdint>
#include <ctgmath>

#ifndef APA_SEG_TYPE
    #define APA_SEG_TYPE uint16_t
#endif

#ifndef APA_FP_TYPE
    #define APA_FP_TYPE double
#endif

#ifndef APA_FP_CACHE_NUMERICAL_ERRORS
    #define APA_FP_CACHE_NUMERICAL_ERRORS false
#endif
namespace apa {

    template<typename ...Ts> void print(Ts... ts){
        (void) std::initializer_list<int>{ (std::cout << ts, 0)...}; }
    template<typename ...Ts> void log(Ts... ts){ print(ts..., '\n'); }

template<typename F>
struct float_info {
// NOTE(kacper):
//       Is there platform agnostic way of acquiring floating point types
//       exponent and mantissa bit counts for a given float type?
//       Since I don't know how to do that I just assumed IEC 559/IEEE 754
//       standard floating points sizes.
private:
    constexpr static size_t _sgn_bits()
    {
        return 1;
    }

    constexpr static size_t _exp_bits()
    {
        switch( sizeof(F) ) {
            case 4:  return 8;
            case 8:  return 11;
            case 10: return 15;
            case 16: return 15;
        }
    }

    constexpr static size_t _mnt_bits()
    {
        switch( sizeof(F) ){
            case 4:  return 23;
            case 8:  return 52;
            case 10: return 64;
            case 16: return 112;
        }
    }

public:
    using byte = unsigned char;
    constexpr static size_t size = sizeof(F) * 8;
    constexpr static size_t bits = sizeof(F) * 8;
    constexpr static size_t sgn_bits = _sgn_bits(); // sign
    constexpr static size_t exp_bits = _exp_bits(); // exponent
    constexpr static size_t mnt_bits = _mnt_bits(); // mantissa
    constexpr static size_t mnt_full_words = mnt_bits / (sizeof(byte) * 8);

    static_assert(bits == sgn_bits + exp_bits + mnt_bits);
};


template<typename T>
constexpr bool is_pow2(T v)
{
    static_assert(std::is_unsigned<T>::value);
    static_assert(std::is_integral<T>::value);
    return v && !(v & (v - 1));
}


inline bool is_little_endian()
{
    const int x { 0x01 };
    const void * addr = static_cast<const void *>(&x);
    return static_cast<const unsigned char *>(addr);
}


template<typename F>
struct float_decomp
{
    // little-endian here!
    uint8_t mnt[float_info<F>::mnt_full_words + 1]; 
    int exp, sgn;
    F norm;
};


template<typename F>
float_decomp<F> float_decompose(F f)
{
    static_assert(std::is_floating_point<F>::value);

    using info = float_info<F>;
    using byte = unsigned char;

    auto n_mask = [] (size_t n) { return ((size_t)1 << n) - 1; };
    float_decomp<F> d;

    d.norm = frexp(f, &d.exp);

    byte * fbyte = reinterpret_cast<byte *>(&d.norm);

    // TODO(kacper): use C++20 std::endian
    if( is_little_endian() ){
        size_t i = 0;
        for(; i < info::mnt_full_words; i++)
          d.mnt[i] = fbyte[i];

        d.mnt[i] = fbyte[i] &
            n_mask(info::mnt_bits - info::mnt_full_words * sizeof(byte) * 8);
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


template<typename _seg_t, typename _fp_t> struct NumAttr {
    using seg_t   = _seg_t;
    using fp_t = _fp_t;

    static_assert(std::is_integral<seg_t>::value);
    static_assert(std::is_unsigned<seg_t>::value);
    static_assert(std::is_floating_point<fp_t>::value);
    static_assert(std::numeric_limits<fp_t>::is_iec559);
    // NOTE(kacper): (the same thing as before)
    //       Is there platform agnostic way of acquiring floating point types
    //       exponent and mantissa bit counts for a given float type?
    //       Since I don't know how to do that I just assumed IEC 559/IEEE 754
    //       standard floating points, hence the assert here - so at least
    //       it won't compile in case of non compatibility

    constexpr static size_t seg_t_size   = sizeof(seg_t);
    constexpr static size_t seg_t_bits   = sizeof(seg_t) * 8;

    constexpr static size_t fp_t_size = sizeof(fp_t);
    constexpr static size_t fp_t_bits = sizeof(fp_t) * 8;

    using fp_t_info = float_info<fp_t>;
};

using DefaultNumAttr = NumAttr<APA_SEG_TYPE, APA_FP_TYPE>;

// NOTE(kacper):
//       zero size (bit count) is now possible... should it be?
template<size_t _SZ, typename _A = DefaultNumAttr>
class _signed {

public:
    using seg_t = typename _A::seg_t;
    using fp_t  = typename _A::fp_t;

    constexpr static size_t seg_t_size = _A::seg_t_size;
    constexpr static size_t seg_t_bits = _A::seg_t_bits;

    constexpr static size_t fp_t_size = _A::fp_t_size;
    constexpr static size_t fp_t_bits = _A::fp_t_bits;

    using fp_t_info = typename _A::fp_t_info;

    constexpr static size_t bits = _SZ;

private:
    constexpr static size_t get_segments_count(){
        size_t full_segs = bits / seg_t_bits;
        if(full_segs * seg_t_bits < bits)
            return full_segs + 1;
        return full_segs;
    }

public:
    constexpr static size_t segments_count = get_segments_count();
    constexpr static size_t real_bits = segments_count * seg_t_bits;

    constexpr static bool pow2_bits = is_pow2(bits);
    constexpr static bool pow2_segs = is_pow2(segments_count);

    enum {
        NEGATIVE  = 1,
        TRUNCATED = 2,
        OVERFLOW  = 4,
        UNDERFLOW = 8
    };

// ======================================= CONSTRUCTORS

    _signed();

    template<typename I,
            typename = typename
            std::enable_if<std::is_integral<I>::value>::type>
    _signed(I val); // for integers

    template<typename F,
            typename = typename
            std::enable_if<std::is_floating_point<F>::value>::type,
            typename = void>
    //      ^ HACK(kacper): so I can have unambigous definition outside class
    _signed(F val); // for floating points

    template<size_t SZ>
    _signed(const _signed<SZ, _A>& other);

    template<typename T>
    bool import(T* data, size_t count);

// ======================================= UTILITY METHODS

    inline void set_bit       (size_t index);
    inline void unset_bit     (size_t index);
    inline void toggle_bit    (size_t index);
    inline void negate        ();
    inline void set_sign      (int s);
    inline void set_sign_bool (bool b);

    inline seg_t   get_segment (size_t index) const;
    inline bool    bit_at      (size_t index) const;
    inline uint8_t get_flags   ()             const;
    inline int8_t  sign        ()             const;
    inline int8_t  sign_bool   ()             const;
    inline bool    trucated    ()             const;
    inline bool    is_zero     ()             const;

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

// ======================================= TYPE CONVERSIONS

    explicit operator bool() const { return !is_zero(); }

// ======================================= ASSIGNMENT OPERATORS

    template<size_t SZ>
    _signed<_SZ, _A>& operator=
        (const _signed<SZ, _A>& other);

// ======================================= ARITHMETIC OPERATORS

    // internal add unsigned helper
    template<size_t SZ1, size_t SZ2, typename A> 
    friend _signed<std::max(SZ1, SZ2)+1, A> add_u
        (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // +
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<std::max(SZ1,SZ2)+1, A> operator+
      (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // +=
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<SZ1, A>& operator+=
        (_signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // internal subtraction helper
    template<size_t SZ1, size_t SZ2, typename A>
    friend _signed<std::max(SZ1,SZ2)+1, A> sub_u
        (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // -
    template<size_t SZ, typename A, typename T>
    friend inline _signed<std::max(SZ,sizeof(T)*8)+1, A> operator-
        (T lhs, const _signed<SZ, A>& rhs);

    // -=
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<SZ1, A>& operator-=
        (_signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // internal multiply unsigned helper
    template<size_t SZ1, size_t SZ2, typename A>
    friend _signed<SZ1+SZ2, A> mul_u
        (const _signed<SZ1, A>& lhs, _signed<SZ2, A> rhs);

    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<SZ1+SZ2, A> operator*
        (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // ~
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A> operator~
        (_signed<SZ, A> lhs);

    // prefix++
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A>& operator++
        (_signed<SZ, A>& lhs);

    // postfix++
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A> operator++
        (_signed<SZ, A>& lhs, int);

    // prefix--
    template<size_t SZ, typename A>
    friend inline _signed<SZ>& operator--
        (_signed<SZ> lhs);

    // postfix--
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A> operator--
        (_signed<SZ, A> lhs, int);

    // <<
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A> operator<<
        (const _signed<SZ, A>& lhs, size_t shift);

    // <<=
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A>& operator<<=
        (_signed<SZ, A>& lhs, size_t shift);

    // >>
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A> operator>>
        (const _signed<SZ, A>& lhs, size_t shift);

    // >>=
    template<size_t SZ, typename A>
    friend inline _signed<SZ, A>& operator>>=
        (_signed<SZ, A>& lhs, size_t shift);

    // &
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<std::min(SZ1, SZ2), A> operator&
        (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // |
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<std::max(SZ1, SZ2), A> operator|
        (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // ^
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline _signed<std::max(SZ1, SZ2), A> operator^
        (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // internal comparison helper
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline int8_t comp_u
        (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // >
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline bool operator>
        (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // <
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline bool operator<
        (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // ==
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline bool operator==
        (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // !=
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline bool operator!=
        (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // <=
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline bool operator<=
        (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

    // >=
    template<size_t SZ1, size_t SZ2, typename A>
    friend inline bool operator>=
        (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs);

private:
    std::array<seg_t, segments_count> _segments = {};
    uint8_t flags = 0;

}; // class _signed<size_t, typename>

template<size_t _SZ>
using s = _signed<_SZ>;

// ================================================= IMPLEMENTATION

// ======================================= UTILITY METHODS
template<size_t _SZ, typename A>
inline typename _signed<_SZ, A>::seg_t _signed<_SZ, A>::get_segment
    (size_t index) const
{
    return (index < segments_count) ? _segments[index] : 0;
}


template<size_t _SZ, typename A>
inline bool _signed<_SZ, A>::bit_at
    (size_t index) const
{
    size_t segment_i = floor(index/seg_t_bits);
    size_t bit_i = index - ( segment_i * seg_t_bits);
    return (get_segment(segment_i) & (1 << bit_i));
}

template<size_t _SZ, typename _A>
inline void _signed<_SZ, _A>::set_bit
    (size_t index)
{
    using seg_t = typename _A::seg_t;
    constexpr size_t seg_t_bits = _A::seg_t_bits;

    if(index >= _signed<_SZ, _A>::bits) return;
    size_t full_segs = index / seg_t_bits;
    size_t residue = index - full_segs * seg_t_bits;
    _segments[full_segs] |= ((seg_t)1 << residue);
}

template<size_t _SZ, typename _A>
inline void _signed<_SZ, _A>::unset_bit
    (size_t index)
{
    using seg_t = typename _A::seg_t;
    constexpr size_t seg_t_bits = _A::seg_t_bits;

    if(index >= _signed<_SZ, _A>::bits) return;
    size_t full_segs = index / seg_t_bits;
    size_t residue = index - full_segs * seg_t_bits;
    _segments[full_segs] &= ~((seg_t)1 << residue);
}

template<size_t _SZ, typename _A>
inline void _signed<_SZ, _A>::toggle_bit
    (size_t index)
{
    using seg_t = typename _A::seg_t;
    constexpr size_t seg_t_bits = _A::seg_t_bits;

    if(index >= _signed<_SZ, _A>::bits) return;
    size_t full_segs = index / seg_t_bits;
    size_t residue = index - full_segs * seg_t_bits;
    _segments[full_segs] ^= ((seg_t)1 << residue);
}

template<size_t _SZ, typename _A>
inline uint8_t _signed<_SZ, _A>::get_flags() const
{
    return flags;
}


template<size_t _SZ, typename _A>
inline int8_t _signed<_SZ, _A>::sign() const
{
    return (flags & NEGATIVE) ? -1 : 1;
}


template<size_t _SZ, typename _A>
inline int8_t _signed<_SZ, _A>::sign_bool() const
{
    return (flags & NEGATIVE) ?  1 : 0;
}


template<size_t _SZ, typename _A>
inline void _signed<_SZ, _A>::set_sign_bool(bool b)
{
    flags &= ~NEGATIVE; flags |= NEGATIVE * b;
}


template<size_t _SZ, typename _A>
inline void _signed<_SZ, _A>::set_sign(int s)
{
    flags &= ~NEGATIVE; flags |= NEGATIVE * (s < 0);
}


template<size_t _SZ, typename _A>
inline bool _signed<_SZ, _A>::trucated() const
{
    return (flags & TRUNCATED);
}


template<size_t _SZ, typename _A>
inline void _signed<_SZ, _A>::negate()
{
    flags ^= NEGATIVE;
}


template<size_t _SZ, typename _A>
inline bool _signed<_SZ, _A>::is_zero() const
{
    for(auto e : _segments) if(e != 0) return false;
    return true; 
}


// count trailing zeros
template<size_t _SZ, typename A>
size_t _signed<_SZ, A>::ctz
    () const
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
size_t _signed<_SZ, A>::clz
    () const
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
_signed<_SZ, _A>::_signed
    ()
{
}


template<size_t _SZ, typename _A>
template<typename I, typename /* SFINAE */>
_signed<_SZ, _A>::_signed
    (I val)
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
_signed<_SZ, _A>::_signed
    (F val)
{
    using info = float_info<F>;
    float_decomp<F> d = float_decompose<F>(val);
    if(d.exp >= 0) {
        import(d.mnt, info::mnt_full_words + 1);
        //NOTE(kacper):
        //      floating points are assumed to be normalized
        //      thus setting implicit first bit
        set_bit(info::mnt_bits);

        // *this >>= 52;
        //NOTE(kacper):
        //      either I do not understand floats or
        //      my exponent is always one too big - hence "+/-1"
        if( (size_t)d.exp > info::mnt_bits )
            *this <<= (d.exp - info::mnt_bits - 1);
        else
            *this >>= (info::mnt_bits - d.exp + 1);
    }
}


template<size_t _SZ, typename _A>
template<size_t SZ>
_signed<_SZ, _A>::_signed
    (const _signed<SZ, _A>& other)
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

// FIXME(kacper): this is dumb byte-byte copy, works only on little-endian
template<size_t _SZ, typename _A>
template<typename T>
bool _signed<_SZ, _A>::import
    (T * data, size_t count)
{
    using byte = unsigned char;
    using seg_t = typename _A::seg_t;
    byte * data_byte_ptr = reinterpret_cast<byte*>(data);
    byte *  seg_byte_ptr = reinterpret_cast<byte*>(&_segments[0]);
    size_t data_bytes = sizeof(T)     * count;
    size_t  seg_bytes = sizeof(seg_t) * _signed<_SZ, _A>::segments_count;

    for(size_t i = 0;
        i < data_bytes && i < seg_bytes; i++) {
        *(seg_byte_ptr + i) = *(data_byte_ptr + i);
    }
    return true;
}

// ======================================= ASSIGNMENT
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

// ============================= ADDITION

// internal add unsigned helper
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

// _signed<SZ1, A>  +  _signed<SZ2, A>
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

// _signed<SZ, A>  +  ArithmeticType
template<size_t SZ, typename A, typename T,
         typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
inline _signed<std::max(SZ,sizeof(T)*8)+1, A> operator+
    (const _signed<SZ, A>& lhs, T rhs)
{
    return operator+(lhs, _signed<sizeof(T) * 8, A>(rhs));
}

// ArithmeticType  +  _signed<SZ, A>
template<size_t SZ, typename A, typename T,
         typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
inline _signed<std::max(SZ,sizeof(T)*8)+1, A> operator+
    (T lhs, const _signed<SZ, A>& rhs)
{
    return operator+(_signed<sizeof(T) * 8, A>(lhs), rhs);
}

// _signed<SZ1, A>  +=  _signed<SZ2, A>
template<size_t SZ1, size_t SZ2, typename A>
inline _signed<SZ1, A>& operator+=
    (_signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    lhs = lhs + rhs;
    return lhs;
}

// _signed<SZ1, A>  +=  ArithmeticType
template<size_t SZ, typename T,
         typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
inline _signed<std::max(SZ,sizeof(T)*8)+1>& operator+=
    (_signed<SZ>& lhs, T rhs)
{
    return lhs += _signed<sizeof(T) * 8>(rhs);
}


// ============================= SUBTRACTION

// internal subtract unsigned helper
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

// _signed<SZ1, A>  -  _signed<SZ2, A>
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

// ArithmeticType  -  _signed<SZ, A>
template<size_t SZ, typename A, typename T>
inline _signed<std::max(SZ, sizeof(T) * 8) + 1, A> operator-
    (T lhs, const _signed<SZ, A>& rhs)
{
    return _signed<sizeof(T)*8, A>(lhs) -  rhs;
}

// _signed<SZ, A>  -  ArithmeticType
template<size_t SZ, typename A, typename T>
inline _signed<std::max(SZ, sizeof(T) * 8) + 1, A> operator-
    (const _signed<SZ, A>& lhs, T rhs)
{
    return lhs - _signed<sizeof(T)*8, A>(rhs);
}

// _signed<SZ1, A>  -=  _signed<SZ2, A>
template<size_t SZ1, size_t SZ2, typename A>
inline _signed<SZ1, A>& operator-=
    (_signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    lhs = lhs - rhs;
    return lhs;
}

// _signed<SZ1, A>  -=  ArithmeticType
template<size_t SZ, typename T,
         typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
inline _signed<std::max(SZ,sizeof(T)*8)+1>& operator-=
    (_signed<SZ>& lhs, T rhs)
{
    return lhs -= _signed<sizeof(T) * 8>(rhs);
}

// ============================= MULTIPLICATION

// FFT

template<typename T>
static constexpr T bit_reverse
    (T in, uint8_t count = sizeof(T) * 8)
{
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
static_assert(bit_reverse(0u) == 0); // is constexpr check

template<typename A, typename I>
static void fast_fourier_transform
    (I first, I last, bool inverse = false)
{
    using fp_t = typename A::fp_t;

    size_t size = last - first;

    for(size_t i = 0; i < size; i++) {
        size_t temp = bit_reverse(i, std::log2(MSB(size - 1)) + 1);

        if(temp < i)
            std::swap( *(first + temp), *(first + i));
    }

    for(int64_t i = std::log2(size) - 1;  i >= 0; i--) {
        for(uint64_t j = 0; j < (1 << i); j++) {
            auto part_sz = size / (1 << i);
            for(uint64_t k = 0; k < size / (1 << (i + 1)); k++) {
                auto w = std::exp(std::complex<fp_t>
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

// internal multiplicatoin helper
template<size_t SZ1, size_t SZ2, typename A>
static _signed<SZ1+SZ2, A> mul_u
    (const _signed<SZ1, A>& lhs, _signed<SZ2, A> rhs)
{
    using fp_t  = typename A::fp_t;
    using seg_t = typename A::seg_t;
    constexpr size_t seg_t_bits = A::seg_t_bits;
    constexpr size_t RET_SZ = SZ1 + SZ2;
    constexpr size_t tmp_arr_sz = _signed<RET_SZ, A>::segments_count;
    constexpr size_t arr_sz = is_pow2(tmp_arr_sz)
        ? tmp_arr_sz : MSB(tmp_arr_sz) << 1;

    static_assert(is_pow2(arr_sz));

    _signed<RET_SZ, A> ret;
    std::array<std::complex<fp_t>, arr_sz> X, Y, Z;

    for(size_t i = 0; i < arr_sz; i++) {
        X[i] = lhs.get_segment(i);
        Y[i] = rhs.get_segment(i);
    }

    fast_fourier_transform<A>(X.begin(), X.end(), false);
    fast_fourier_transform<A>(Y.begin(), Y.end(), false);

    for(size_t i=0; i<arr_sz; i++) Z[i] = X[i] * Y[i];

    fast_fourier_transform<A>(Z.begin(), Z.end(), true);

    // during 'reassemly' of the number a little bit more space is needed
    constexpr size_t TMP_SZ = sizeof(seg_t) * RET_SZ;

    _signed<TMP_SZ, A> temp = 0;

    for(size_t i = 0; i < arr_sz; i++) {
        //TODO(kacper): tidy it up after implementing floating point assignment
        temp = _signed<TMP_SZ, A>( std::round(Z[i].real()) );
        //TODO(kacper): make it one bitshift (12% in perf report)
        temp <<= i * seg_t_bits;
        //TODO(kacper): constexpr this log
        temp >>= std::log2(arr_sz);
        ret += temp;
    }
    return ret;
}

// _signed<SZ1, A>  *  _signed<SZ2, A>
template<size_t SZ1, size_t SZ2, typename A>
inline _signed<SZ1+SZ2, A> operator*
    (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    auto ret = mul_u(lhs,rhs);
    ret.set_sign_bool( lhs.sign() xor rhs.sign() );
    return ret;
}

// _signed<SZ1, A>  *  ArithmeticType
template<size_t SZ, typename A, typename T>
inline _signed<SZ+sizeof(T)*8, A> operator*
    (const _signed<SZ, A>& lhs, T rhs)
{
    return operator*(lhs, _signed<sizeof(T)*8, A>(rhs));
}

// _signed<SZ1, A>  *  ArithmeticType
template<size_t SZ, typename A, typename T>
inline _signed<SZ+sizeof(T)*8, A> operator*
    (T lhs, const _signed<SZ, A>& rhs)
{
    return operator*(_signed<sizeof(T)*8, A>(lhs), rhs);
}


// ============================= DIVISION

// ======================================= UNIARY OPERATORS

// operator~
template<size_t _SZ, typename A>
inline _signed<_SZ, A> operator~
    (_signed<_SZ, A> lhs)
{
    //TODO: truncation?
    for(auto& s : lhs._segments) s = ~s;
    return lhs;
}

// operator prefix++;
template<size_t _SZ, typename A>
inline _signed<_SZ, A>& operator++
    (_signed<_SZ, A>& lhs)
{
    lhs = (lhs+1);
    return lhs;
}

// operator postfix++;
template<size_t _SZ, typename A>
inline _signed<_SZ, A> operator++
    (_signed<_SZ, A>& lhs, int)
{
    _signed<_SZ, A> ret(lhs);
    ++lhs;
    return ret;
}

// operator prefix--
template<size_t _SZ, typename A>
inline _signed<_SZ, A>& operator--
    (_signed<_SZ, A> lhs)
{
    lhs = (lhs-1);
    return lhs;
}

// operator postfix--
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


// ======================================= BINARY OPERATORS
// operator<<
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
        // NOTE(kacper):
        //       this conditional is here due to bitshift operator restriction
        //       to range [0, sizeof(T)*8) on x86_64 (Well its
        //       the shift modulo register size but offically undefined behaviour)
        // TODO(kacper):
        //      take it out of the loop ( though compiler did it for me probably )
        seg_t seglo = seglo_shift < seg_t_bits
            ? lhs.get_segment(i - seg_dist - 1) >> seglo_shift
            : 0;

        seg_t seghi = lhs.get_segment(i - seg_dist) << seghi_shift;

        ret._segments[i] = seglo | seghi;
    }
    return ret;
}


// operator<<=
template<size_t SZ, typename A>
inline _signed<SZ, A>& operator<<=
    (_signed<SZ, A>& lhs, size_t shift)
{
    lhs = lhs << shift;
    return lhs;
}


// operator>>
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

        // NOTE(kacper): lookup comment above, the same thing applies
        seg_t seghi = seghi_shift < seg_t_bits
            ? lhs.get_segment(i + seg_dist + 1) << seghi_shift
            : 0;
        ret._segments[i] = seglo | seghi;
    }
    return ret;
}


// operator>>=
template<size_t SZ, typename A>
inline _signed<SZ, A>& operator>>=
    (_signed<SZ, A>& lhs, size_t shift)
{
    lhs = lhs >> shift;
    return lhs;
}


// operator&
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





// _signed<SZ1, A>  |  _signed<SZ2, A>
template<size_t SZ1, size_t SZ2, typename A>
inline _signed<std::max(SZ1, SZ2)> operator|
    (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    _signed<std::max(SZ1, SZ2), A> ret;
    for(size_t i=0; i < std::min(lhs.segments_count, rhs.segments_count); i++)
        ret._segments[i] = lhs.get_segment(i) | rhs.get_segment(i);
    return ret;
}

// _signed<SZ1, A>  |  ArithmeticType
template<size_t SZ, typename T, typename A>
inline _signed<std::max(SZ, sizeof(T)*8), A> operator|
    (const _signed<SZ, A>& lhs, const T rhs)
{
    return lhs & _signed<(sizeof(T) * 8)>(rhs);
}

// ArithmeticType  |  _signed<SZ2, A>
template<size_t SZ, typename T, typename A>
inline _signed<std::max(SZ, sizeof(T) * 8), A> operator|
    (const T rhs, const _signed<SZ, A>& lhs)
{
    return lhs & _signed<(sizeof(T) * 8)>(rhs);
}

// _signed<SZ1, A>  ^  _signed<SZ2, A>
template<size_t SZ, typename T, typename A>
inline _signed<std::max(SZ, sizeof(T)*8)> operator&
    (const _signed<SZ, A>& lhs, const T rhs)
{
    return lhs & _signed<(sizeof(T) * 8), A>(rhs);
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


// ======================================= RELATIONAL OPERATORS

// internal comparison helper
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


// _signed<SZ1, A>  >  _signed<SZ2, A>
template<size_t SZ1, size_t SZ2, typename A>
inline bool operator>
    (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    if(lhs.sign_bool() xor rhs.sign_bool())
        return lhs.sign() > 0;
    else
        return (lhs.sign() * comp_u(lhs, rhs)) > 0;
}


// _signed<SZ1, A>  <  _signed<SZ2, A>
template<size_t SZ1, size_t SZ2, typename A>
inline bool operator<
    (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    if(lhs.sign_bool() xor rhs.sign_bool())
        return lhs.sign() < 0;
    else
        return (lhs.sign() * comp_u(lhs, rhs)) < 0;
}


// _signed<SZ1, A>  ==  _signed<SZ2, A>
template<size_t SZ1, size_t SZ2, typename A>
inline bool operator==
    (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    if(lhs.sign() == rhs.sign())
        return comp_u(lhs, rhs) == 0;
    else
        return lhs.is_zero() && rhs.is_zero();
}


// _signed<SZ1, A>  !=  _signed<SZ2, A>
template<size_t SZ1, size_t SZ2, typename A>
inline bool operator!=
    (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    return !(lhs == rhs);
}


// _signed<SZ1, A>  <=  _signed<SZ2, A>
template<size_t SZ1, size_t SZ2, typename A>
inline bool operator<=
    (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    return (lhs < rhs || lhs == rhs);
}


// _signed<SZ1, A>  >=  _signed<SZ2, A>
template<size_t SZ1, size_t SZ2, typename A>
inline bool operator>=
    (const _signed<SZ1, A>& lhs, const _signed<SZ2, A>& rhs)
{
    return (lhs > rhs || lhs == rhs);
}

} //namespace apa

