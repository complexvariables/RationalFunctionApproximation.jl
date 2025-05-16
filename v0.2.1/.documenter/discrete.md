
# Discrete data {#Discrete-data}

The original AAA algorithm ([Nakatsukasa, Sète, Trefethen 2018](https://epubs.siam.org/doi/abs/10.1137/16M1106122)) works with a fixed set of points on the domain of approximation. There is a legacy `aaa` function that can work with this type of data:

```julia
using RationalFunctionApproximation, ComplexRegions
x = -1:0.01:1
f = x -> tanh(5 * (x - 0.2))
r = aaa(x, f.(x))
```


```
Barycentric{Float64, Float64} rational interpolant of type (35, 35):
    -0.999988=>-1.0,  0.999329=>1.0,  -0.99186=>-0.35,  …  -0.999727=>-0.69
```


However, it&#39;s preferable to use the `approximate` function for this purpose, as the result type is more useful within the package. Simply pass the function and, in the form of a vector, the domain.

```julia
r = approximate(f, x)
```


```
Barycentric{Float64, Float64} rational interpolant of type (11, 11) on the domain: -1.0:0.01:1.0
```


As long as there are no singularities as close to the domain as the sample points are to one another, this fully discrete approach should be fine:

```julia
I = unit_interval
println("nearest pole is $(minimum(dist(z, I) for z in poles(r))) away")
_, err = check(r);
println("max error on the given domain: ", maximum(abs, err))
```


```
nearest pole is 0.31415926535543987 away
[ Info: Max error is 2.65e-14
max error on the given domain: 2.653433028854124e-14
```


But if the distance to a singularity is comparable to the sample spacing, the quality of the approximation may suffer. Even worse, the method may not be aware that it has failed.

```julia
f = x -> tanh(400 * (x - 0.2))
r = approximate(f, x)
println("nearest pole is $(minimum(dist(z, I) for z in poles(r))) away")
_, err = check(r);
println("max error on the given domain: ", maximum(abs, err))
err = maximum(abs(f(x)- r(x)) for x in range(-1, 1, 3000))
println("max error on finer test points: ", err)
```


```
n = 2  idx_max = 1  i_node = CartesianIndex(1, 1)
n = 3  idx_max = 125  i_node = CartesianIndex(127, 1)
n = 4  idx_max = 15  i_node = CartesianIndex(16, 1)
n = 5  idx_max = 150  i_node = CartesianIndex(154, 1)
n = 6  idx_max = 196  i_node = CartesianIndex(201, 1)
n = 7  idx_max = 120  i_node = CartesianIndex(122, 1)
n = 8  idx_max = 177  i_node = CartesianIndex(183, 1)
n = 9  idx_max = 92  i_node = CartesianIndex(94, 1)
n = 10  idx_max = 117  i_node = CartesianIndex(120, 1)
n = 11  idx_max = 113  i_node = CartesianIndex(116, 1)
n = 12  idx_max = 117  i_node = CartesianIndex(123, 1)
n = 13  idx_max = 108  i_node = CartesianIndex(111, 1)
n = 14  idx_max = 131  i_node = CartesianIndex(141, 1)
n = 15  idx_max = 113  i_node = CartesianIndex(118, 1)
n = 16  idx_max = 77  i_node = CartesianIndex(79, 1)
n = 17  idx_max = 120  i_node = CartesianIndex(132, 1)
n = 18  idx_max = 147  i_node = CartesianIndex(162, 1)
n = 19  idx_max = 100  i_node = CartesianIndex(104, 1)
n = 20  idx_max = 59  i_node = CartesianIndex(61, 1)
n = 21  idx_max = 108  i_node = CartesianIndex(115, 1)
n = 22  idx_max = 171  i_node = CartesianIndex(191, 1)
n = 23  idx_max = 45  i_node = CartesianIndex(47, 1)
n = 24  idx_max = 114  i_node = CartesianIndex(130, 1)
n = 25  idx_max = 105  i_node = CartesianIndex(113, 1)
n = 26  idx_max = 118  i_node = CartesianIndex(137, 1)
n = 27  idx_max = 146  i_node = CartesianIndex(169, 1)
n = 28  idx_max = 32  i_node = CartesianIndex(34, 1)
n = 29  idx_max = 99  i_node = CartesianIndex(107, 1)
n = 30  idx_max = 106  i_node = CartesianIndex(121, 1)
n = 31  idx_max = 171  i_node = CartesianIndex(200, 1)
n = 32  idx_max = 82  i_node = CartesianIndex(88, 1)
n = 33  idx_max = 120  i_node = CartesianIndex(145, 1)
n = 34  idx_max = 104  i_node = CartesianIndex(119, 1)
nearest pole is 0.004481098151799039 away
[ Info: Max error is 1.67e-14
max error on the given domain: 1.6653345369377348e-14
max error on finer test points: 0.0507185469433139
```


In the continuous mode, the adaptive sampling of the domain attempts to ensure that the approximation is accurate everywhere.

```julia
r = approximate(f, I; tol=1e-12)
err = maximum(abs(f(x)- r(x)) for x in range(-1, 1, 3000))
println("max error on finer test points: ", err)
```


```
max error on finer test points: 7.924771949774367e-13
```

