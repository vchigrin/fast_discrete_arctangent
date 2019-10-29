## Fast discrete version of std::atan2 function.

Outputs integral "sector number" for floating-point (x,y) coords of point.
Pre-computes table at initialization time for fast processing time.

Tests results on Intel Core i5-2450M CPU @ 2.50GHz and clang 6.0 compiler.
Tested on 20 millions points, randomly uniformly spreaded in
[(-5,-5), (5,5)] square.
Two different methods can give responses, that differ by 1 in very rare cases,
because of rounding errors.
In our tests we divide 360-degrees circle to 2000 sectors.

#### std::atan2 based solution (double)
- Time to process 20 millions points 1227 ms.

#### std::atan2 based solution (float)
- Time to process 20 millions points 814 ms.
- Number of rounding errors compared to simple std::atan2(double) based solution - 1657.

#### Our table based solution (float)
- Time to process 20 millions points 227 ms **72% faster**.
- Number of rounding errors compared to simple std::atan2(double) based solution - *97* -
    even *better* then std::atan2 based solution for single-precision float.

#### Our table based solution (double)
- Time to process 20 millions points 240 ms **80% faster**.
- Number of rounding errors compared to simple std::atan2(double) based solution - *0*.
