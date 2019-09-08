# APA -Arbitrary Precision Arithmetic
## Signed Integer
APA is yet another library that lets you perform calculations on numbers beyond cpu word size.
It is self contained except for stl dependency which is planned to be optional or dropped entirely as its not crucial at all. It is tested against well known gmp with catch2.
Only available type at the moment is signed integer implemented as template with 
compile-time known size which I thought may be a good idea, since it gives more information
to the compiler and because usually, for example in crypto and hashing,
final result size is known beforehand. Because of that operators return values 
big enough contain maximal result of given operation, for example adding integers of 
size m and n will return integer of size max(m,n) + 1 and multiplying them will result in 
integer of size m+n thus be careful with auto.

### Floating Point Numbers
Low level floating point types conversion is provided without intermediate cast to integer
so you can assign floating point value greater than maximum integer value.
This basically means that floating point type is punned, 
mantissa is read directly from memory and shifted by exponent value. 
This was tested only on little-endian cpu though should work also on big-endian.


Addition and subtraction are performed in O(n) with schoolbook algorithms.
Multiplication uses non-recursive Fast Fourier Transform to run in O(n*log(n)) 
though for small integer sizes worse time complexity algorithms may perform better (TODO!).
Division


#TODO
* Implement worse time complexity algorithms that perform better for smaller N.
* Allow for non power of two integer sizes.
* Acquire compile-time randomness for testing so integer sizes can be randomized.
* Implement function returning last operation maximum error. (Numerical error are mostly negleble but present)
* Big-endian support
* operator overloads between types with different attributes

 
