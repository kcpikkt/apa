

#include "bigint.h"
#include <complex>
#include <vector>

namespace bigint{
// Constructors
//
template<size_t _SZ> 
Signed<_SZ>::Signed()
{}

template<size_t _SZ>
template<typename T> 
Signed<_SZ>::Signed(T val){
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

// utility
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

template<size_t SZ1, size_t SZ2> 
/*static */Signed<std::max(SZ1, SZ2)+1> add_u(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    constexpr size_t ret_sz = std::max(SZ1, SZ2)+1;
    Signed<ret_sz> ret;
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
}

template<size_t SZ1, size_t SZ2> 
Signed<std::max(SZ1,SZ2)+1> sub_u(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    constexpr size_t ret_sz = std::max(SZ1, SZ2)+1;
    Signed<ret_sz> ret;
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
        ret = add_u(ret,Signed<8>(1) );
        ret.toggle_sign();
    }
    return ret;
}
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
    if(lhs.sign() == rhs.sign())    { return add_u<SZ1,SZ2>(lhs, rhs);}
    else                            {
        if(rhs.is_negative())         return sub_u<SZ1,SZ2>(lhs, rhs);
        else                          return sub_u<SZ1,SZ2>(rhs, lhs);
    }
}
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
    if(lhs.sign() == rhs.sign())    { return sub_u<SZ1,SZ2>(lhs, rhs);}
    else                            { return add_u<SZ1,SZ2>(lhs, rhs);}
}

   
template<typename Iter>
void fft(Iter first, Iter last, bool inverse = false){
    size_t size = last - first;
    if(size >= 2){
        std::vector<std::complex<double>> temp(size/2);
        for(size_t i=0; i<size/2; i++){
            temp [i] = first[i * 2 + 1];
            first[i] = first[i * 2];
        }
        for(size_t i=0; i<size/2; i++)
            first[i + size/2] = temp[i];

        auto split = first + size/2;
        fft(first,split);
        fft(split,last);
        for(size_t k=0; k<size/2; k++){
            auto w = std::exp(std::complex<double>(0, (inverse ? 1 : 1) * 2.0 * M_PI * k / size));
            auto bottom = first[k];
            auto top = first[k + size/2];

            first[k]          = bottom + w * top;
            first[k + size/2] = bottom - w * top;
        }
    }
}

template<typename Iter>
void ifft(Iter first, Iter last){
    fft(first, last, true);
    size_t size = last - first;
    for(auto& el=first; el!=last; el++)
        (*el).real((*el).real()/size);
}
template<size_t SZ1, size_t SZ2> 
Signed<SZ1+SZ2> mul_u(const Signed<SZ1>& lhs, Signed<SZ2> rhs){
    //TODO: check for overflow
    //TODO: pow2
    constexpr size_t pow2_sz = Signed<SZ1+SZ2>::segments_count << 1;
    lhs.print_segs();

    Signed<pow2_sz * impl_t_bit_sz> ret; 
    std::array<std::complex<double>, pow2_sz> X, Y, Z;

    for(size_t i=0; i<pow2_sz; i++) X[i] = lhs.get_segment(i);
    for(size_t i=0; i<pow2_sz; i++) Y[i] = rhs.get_segment(i);

    auto pp = [](std::array<std::complex<double>, pow2_sz> A){ for(size_t i=0; i<pow2_sz; i++){ std::cout << A[i] << " "; } std::puts("\n");};

    pp(X);
    fft(X.begin(), X.end(), false);
    pp(X);

    std::cout << std::endl;
    pp(Y);
    fft(Y.begin(), Y.end(), false);
    pp(Y);

    for(size_t i=0; i<pow2_sz; i++) Z[i] = X[i] * Y[i];

    std::cout << std::endl;
    pp(Z);
    ifft(Z.begin(), Z.end());
    pp(Z);

    for(size_t i=0; i<pow2_sz; i++) ret._segments[i] = ((impl_t)(Z[i].real() / 8));

    std::cout << ret.binary_string() << std::endl;
    // std::cout << lhs.binary_string() << std::endl;
    // std::cout << rhs.binary_string() << std::endl;


    return ret;
}
//template<size_t SZ, typename T>
//inline Signed<SZ+sizeof(T)*8> operator*(T lhs, const Signed<SZ>& rhs){
//    return operator*(Signed<sizeof(T)*8>(lhs), rhs);
//}
//
//template<size_t SZ, typename T>
//inline Signed<SZ+sizeof(T)*8> operator*(const Signed<SZ>& lhs, T rhs){
//    return operator*(lhs, Signed<sizeof(T)*8>(rhs));
//}
template<size_t SZ1, size_t SZ2>
inline Signed<SZ1+SZ2> operator*(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    //FIXME
    //ret.set_sign(lhs.sign() xor rhs.sign());

    return mul_u(lhs,rhs);
}

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
////operator uniary+
//template<size_t SZ>
//inline Signed<SZ> operator+(Signed<SZ> lhs){
//    lhs.set_sign(false);
//    return lhs;
//}
////operator uniary-
//template<size_t SZ>
//inline Signed<SZ> operator-(Signed<SZ> lhs){
//    lhs.set_sign(true);
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
inline bool comp_u(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){ //is lhs greater
    for(size_t i = std::max(lhs.segments_count, rhs.segments_count); i > 0; i--){
        auto lhs_seg = lhs.get_segment(i-1);
        auto rhs_seg = rhs.get_segment(i-1);
        if((bool)lhs_seg xor (bool)rhs_seg) return (lhs_seg == 0);
        else                                if (lhs_seg > rhs_seg) return true;
    }
    return false;
}
//operator>
template<size_t SZ1, size_t SZ2>
inline bool operator>(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    if(lhs.sign() == rhs.sign())
        return (lhs.is_positive() ? comp_u(lhs, rhs) : comp_u(rhs, lhs));
    else{
        if(lhs.is_zero() && rhs.is_zero()) return false;
        else return lhs.is_positive();
    }
}
//operator<
template<size_t SZ1, size_t SZ2>
inline bool operator<(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    if(lhs.sign() == rhs.sign())
        return (lhs.is_positive() ? comp_u(rhs, lhs) : comp_u(lhs, rhs));
    else{
        if(lhs.is_zero() && rhs.is_zero()) return false;
        else return rhs.is_positive();
    }
}
//operator==
template<size_t SZ1, size_t SZ2>
inline bool operator==(const Signed<SZ1>& lhs, const Signed<SZ2>& rhs){
    if( lhs.sign() != rhs.sign()) 
        return lhs.is_zero() && rhs.is_zero();

    for(size_t i = 0; i < std::max(lhs.segments_count, rhs.segments_count); i++)
        if(lhs.get_segment(i-1) != rhs.get_segment(i-1)) return false;
    return true;
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
}
