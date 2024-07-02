# Calling from Python

You can call the functions in this package from Python using the [`PythonCall`/`JuliaCall`](https://juliapy.github.io/PythonCall.jl/stable/) package. 

## Installation

It's recommended to create a new virtual environment to try this out. In Python, you need to install `juliacall` via

```bash
pip install juliacall
```

Then you start Python and run:

```python
from juliacall import Main as jl
```

This will download and initialize a copy of Julia. Finally, you need to install this package in that Julia environment:

```python
jl.seval('using Pkg; Pkg.add("RationalFunctionApproximation")')
```

That should be all you need in the Python environment.

## Usage

In each new Python session, you need to load the packages:

```python
from juliacall import Main as jl
jl.seval('using RationalFunctionApproximation, PythonCall')
```

All the functions and constants exposed to Julia by this package are available using the `jl` object. For example, to use the discrete AAA algorithm:

```python
import numpy as np    # if installed in Python
x = np.linspace(-1, 1, 1000)
y = np.tanh(5 * (x - 0.2))
r = jl.aaa(x, y)
print(r)
```

```
Barycentric rational function of type (11,11)
```

This will return a wrapped Julia object that you can use in Python as if it were a Python object. For example, you can evaluate the approximation at a point:

```python
r(0.5)
```

```
0.9051482536448658
```

If you want to apply the function at multiple points, you can use comprehensions or vectorize it in numpy:

```python
rv = np.vectorize(r)
rv(np.array([0.5, 0.6, 0.7]))
```

```
array([0.90514825, 0.96402758, 0.9866143 ])
```

You can get information about the approximation using any documented function in the package, e.g.:

```python
print(jl.poles(r))    # returns wrapped Julia type
```

```
ComplexF64[0.20000000000544785 - 0.31415926535542893im, 0.20000000000544788 + 0.31415926535542893im, 0.20000207991810143 - 0.942477292594254im, 0.20000207991810143 + 0.9424772925942541im, 0.20308324780986833 - 1.5724812056318853im, 0.20308324780986833 + 1.5724812056318853im, 0.29268586746842673 - 2.3408220889660796im, 0.29268586746842673 + 2.34082208896608im, 0.9695028397625358 + 4.390786420000105im, 0.969502839762536 - 4.390786420000105im, 21.59156666159181 + 0.0im]
```

``` python
print(np.array(jl.residues(r)))    # converts to numpy array
```

```
[  0.2       +2.72029915e-11j   0.2       -2.72031942e-11j
   0.19999893+6.55140711e-06j   0.19999893-6.55140637e-06j
   0.20352821+4.86975387e-03j   0.20352821-4.86975387e-03j
   0.33454619+6.91112099e-02j   0.33454619-6.91112099e-02j
   1.25164001-5.59634589e-01j   1.25164001+5.59634589e-01j
 -32.51419889+0.00000000e+00j]
```

## Passing Python functions

To use continuous approximation, you can pass a Python function to the `approximate` function.

```python
def f(x):
    return np.tanh(5 * (x - 0.2))

r = jl.approximate(f, jl.unit_interval)
r(.5)
```

```
0.9051482536448647
```
