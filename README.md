# APA - C++17 Arbitrary Precision Arithmetic Library
APA is yet another library that lets you perform calculations on numbers beyond cpu word size.
It is self contained except for stl dependency which is planned to be optional or dropped entirely as its not crucial at all. It is tested against well known gmp with catch2.

What is perhaps different about it is that Signed Integer type is implemented as a class template 
with bit width being template parameter and therefore compile-time known.
Motivated by the fact that in many use cases, like cryptography, 
interger sizes are fixed and known a priori, but its yet to be tested how good/bad of an idea that was
(if optimized, probably saving few conditional jumps and heap allocation in comparison to standard way 
of doing things, though it wasn't created with that premise but rather just as an experiment).

## Casts from Floating Points
Low level floating point types conversion is provided without intermediate cast to integer
so you can assign floating point value greater than maximum integer value.
This basically means that floating point type is punned, 
mantissa is read directly from memory and shifted by exponent value (IEC 559/IEEE 754 assumed). 
This was tested only on little-endian cpu though should work also on big-endian.

Addition and subtraction are performed in O(n) with schoolbook algorithms.
Multiplication uses non-recursive Fast Fourier Transform to run in O(n*log(n)) 
though for small integer sizes worse time complexity algorithms may perform better (TODO).
Division, Modulo and thus also casts from/to string literals are not yet implemented. 
Challange I gave to myself here was to implement the best complexity algorithms (at least at first)
but papers about fast division algorithms proved to be quite hard to get throgh for me (TODO!).

# Building
```
    $ cd build/
    $ ccmake .. # set compiler to clang++
    $ cmake .
    $ make 
```

# TODO
* Implement worse time complexity algorithms that perform better for smaller N.
* Fake non-power of 2 integer widths with truncation and stuff.
* Acquire compile-time randomness for testing so integer sizes can be randomized. (kind of macro hacky right now)
* Add function returning last operation maximum error. (Numerical errors are mostly neglible but present, e.g. in FFT)
* Big-endian compatibility.
* operator overloads between types with different attributes (which are basically underlying unsigned integer type 
and floating point used for calculations) (yeah, good luck with that template mess).
* Don't use CMake

# Notes
GMP and Catch2 are both fully included in repo under test/include because I have not yet found simple way to tell CMake
to find their headers in platform independent way.
