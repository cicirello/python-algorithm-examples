# python-algorithm-examples: Examples of algorithms in Python
# 
# Copyright (c) 2023 Vincent A Cicirello
# https://www.cicirello.org/
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import timeit

def fibonacci_exponential(n):
    """Computes the n-th number in the Fibonacci sequence, using the naive
    recursive algorithm implied by the definition of the Fibonacci sequence.

    Note: For simplicity, this ignores negative n, for which the basic form of
    the sequence is undefined.

    Runtime (assuming additions are constant time):
    A simplistic analysis of the runtime leads to O(2**n). The number of
    additions is O(2**n). If we make the naive assumption that the integers
    we are adding are small (e.g., 32-bit integers), then we can assume that
    addition is a constant time operation. Thus, the runtime of O(2**n) under
    this assumption. Note that this assumption holds for n <= 46.

    Runtime (no assumptions on magnitude of n):
    For large enough n (specifically n > 46), the additions we'll encounter
    will not be a constant time operation. Python ints are objects that can
    represent integers of arbitrary magnitude. Addition of two m-bit integers
    requires O(m) time. The Fibonacci sequence grows exponentially. We can
    bound the n-th Fibonacci number Fib_n above by O(2**n). The number of bits
    m needed to represent Fib_n is thus m = O(lg Fib_n) = O(lg (2**n)) = O(n).
    Thus, the runtime for the additions is linear in n, and not a constant time
    operation.
    
    Thus the runtime is O(n * 2**n), if we allow for computing Fibonacci
    numbers of arbitrary magnitude.

    Keyword arguments:
    n - computes the n-th Fibonacci number.
    """
    if n <= 1:
        return n
    return fibonacci_exponential(n-1) + fibonacci_exponential(n-2)

def fibonacci_linear(n):
    """Computes the n-th number in the Fibonacci sequence, using a simple
    linear iteration over the sequence.

    Note: For simplicity, this ignores negative n, for which the basic form of
    the sequence is undefined.

    Runtime (assuming additions are constant time):
    A simplistic analysis of the runtim is O(n). It performs O(n) additions.
    If we make the naive assumption that the integers we are adding are small
    (e.g., 32-bit integers), then we can assume that addition is a constant
    time operation. Thus, the runtime of O(n) under this assumption. Note
    that this assumption holds for n <= 46.
    
    Runtime (no assumptions on magnitude of n):
    For large enough n (specifically n > 46), the additions we'll encounter
    will not be a constant time operation. Python ints are objects that can
    represent integers of arbitrary magnitude. Addition of two m-bit integers
    requires O(m) time. The Fibonacci sequence grows exponentially. We can
    bound the n-th Fibonacci number Fib_n above by O(2**n). The number of bits
    m needed to represent Fib_n is thus m = O(lg Fib_n) = O(lg (2**n)) = O(n).
    Thus, the runtime for the additions is linear in n, and not a constant time
    operation.

    Thus, the runtime of this "linear" algorithm is not so linear, and is
    actually quadratic: O(n**2), since we are doing O(n) additions, and
    additions cost O(n). Only the low-n case is linear.

    Keyword arguments:
    n - computes the n-th Fibonacci number.
    """
    if n <= 1:
        return n
    fib_i_minus_1, fib_i = 0, 1
    for i in range(1, n):
        fib_i_minus_1, fib_i = fib_i, fib_i_minus_1 + fib_i
    return fib_i

def fibonacci_logarithmic(n):
    """Computes the n-th number in the Fibonacci sequence in logarithmic time,
    using an algorithm based on the following matrix equation:
    [F_{n+1}  F_n    ] = [1  1] ** n
    [F_n      F_{n-1}]   [1  0]
    from Donald Knuth's The Art of Computer Programming, Volume 1, Fundamental
    Algorithms, 3rd Edition.

    Note: For simplicity, this ignores negative n, for which the basic form of
    the sequence is undefined.

    Runtime (assuming addition and multiplication are constant time):
    A simplistic analysis of the runtim is O(lg n). It performs O(lg n)
    additions and mulitplications. If we make the naive assumption that the
    integers we are adding and multiplying are small (e.g., 32-bit integers),
    then we can assume that addition and multiplication are constant time
    operations. Thus, the runtime is O(lg n) under this assumption. Note
    that this assumption holds for n <= 46.

    Runtime (no assumptions on magnitude of n):
    For large enough n (specifically n > 46), the additions and multiplications
    will not be constant time operations. Python ints are objects that can
    represent integers of arbitrary magnitude. Addition of two m-bit integers
    requires O(m) time. Mulitplication is more costly than addition. The runtime
    to multiply two m-bit integers depends upon the algorithm used. Currently,
    Python uses Karatsuba's multiplication algorithm when multiplying very
    large integers. The runtime of Karatsuba's algorithm is approximately
    O(m**1.58) to multiply m-bit integers. The Fibonacci sequence grows
    exponentially. We can bound the n-th Fibonacci number Fib_n above by
    O(2**n). The number of bits m needed to represent Fib_n is thus m =
    O(lg Fib_n) = O(lg (2**n)) = O(n). Thus, the runtime for the
    multiplications is O(n**1.58).

    Thus, the runtime of this "logarithmic" time algorithm is not so
    logarithmic. Its runtime is O(n**1.58 lg(n)). It performs O(lg(n))
    additions that cost O(n), for total cost of additions: O(n log n). It
    performs O(lg n) multiplications that cost O(n**1.58) for total cost of
    multiplications: O(n**1.58 lg(n)). Thus, the runtime is O(n**1.58 lg(n)).
    Only the low-n case is logarithmic.

    Note that the O(n**1.58 lg(n)) of this algorithm is lower order than the
    quadratic runtime of the "linear" algorithm. Thus, this "logarithmic"
    algorithm should be the fastest of the three considered here for large
    values of n. However, the constant factors that O() hides are more
    substantial for this "logarithmic" algorithm than for the "linear" one, so
    you will likely find the "linear" algorithm faster than this one for the
    low n case.

    Keyword arguments:
    n - computes the n-th Fibonacci number.
    """
    if n <= 1:
        return n
    # Initialize the following implied matrix to the identity:
    # [g  f] = [1  0]
    # [f  h] = [0  1]
    # Note that we don't actually need an array-like structure.
    # When this function returns, f will be the n-th Fibonacci
    # number.
    g = h = 1
    f = 0
    # Initialize the following implied matrix with:
    # [x  y] = [1  1]
    # [y  z] = [1  0]
    # Note that this corresponds to Knuth's matrix equation
    # for n = 1. We're going to obtain the O(lg n) time by
    # computing the n-th power using a square and multiple
    # approach. Also note that we only need 3 variables because
    # all powers of this specific matrix are such that the minor
    # diagonal elements are identical.
    x = y = 1
    z = 0
    while True:
        if n & 1 == 1:
            # if right-most bit a 1, multiply the matrices
            fy = f * y
            f, g, h = f * x + h * y, fy + g * x, fy + h * z
        # iterating over bits of n
        n >>= 1
        if n == 0:
            break
        # square the Knuth matrix
        y2 = y * y
        x, y, z = x * x + y2, x * y + y * z, z * z + y2
    return f

def time_fibonacci_low_n(repetitions = 1000000, time_cutoff_exponential = 2.0):
    """Runs a timing experiment to deminstrate runtime of the three algorithms
    on low values of n, i.e., values of n such that the corresponding Fibonacci
    number is representable by a 32-bit signed int. Specifically, this considers
    values of n from 1 to 46.

    Keyword arguments:
    repetitions - this value is passed as number to timeit. The "linear" and
        "logarithmic" algorithms are quite fast, so the default is set high
        at 1000000 to get better measurements.
    time_cutoff_exponential - once the time for the exponential algorithm
        exceeds this value, it will just compare the other two.
    """
    time_fibonacci(
        n_min = 1,
        n_max = 46,
        n_update = lambda n_old : n_old + 1,
        repetitions = repetitions,
        time_cutoff_exponential = time_cutoff_exponential
    )
    
def time_fibonacci(
    n_min = 1,
    n_max = 1048576,
    n_update = lambda n_old : 2 * n_old,
    repetitions = 1,
    time_cutoff_exponential = 0.4):
    """Runs a timing experiment to demonstrate runtime of the three
    algorithms.

    Keyword arguments:
    n_min - first value of n to consider.
    n_max - consider cases of n no greater than this value.
    n_update - a function to determine next value of n based on current
        value of n. This function is assumed to increase values of n.
    repetitions - this value is passed as number to timeit. Be careful what
        you use here as the times grow quickly. The default of 1 should be
        fine if you are considering large values of n.
    time_cutoff_exponential - once the time for the exponential algorithm
        exceeds this value, it will just compare the other two.
    """
    n = n_min
    include_exponential = True
    template_all = "{0:11d} {1:11.8f} {2:11.8f} {3:11.8f}"
    template_no_exp = "{0:11d} {1:>11s} {2:11.8f} {3:11.8f}"
    print("{0:>11s} {1:>11s} {2:>11s} {3:>11s}".format(
        "n",
        "exponential",
        "linear",
        "logarithmic"
        )
    )
    while n <= n_max:
        template = template_no_exp 
        t_exponential = "N/A"
        if include_exponential:
            t_exponential = timeit.timeit(
                lambda : fibonacci_exponential(n),
                number = repetitions)
            template = template_all
            include_exponential = t_exponential < time_cutoff_exponential
        t_linear = timeit.timeit(
            lambda : fibonacci_linear(n),
            number=repetitions)
        t_logarithmic = timeit.timeit(
            lambda : fibonacci_logarithmic(n),
            number=repetitions)
        print(template.format(n, t_exponential, t_linear, t_logarithmic))
        n = n_update(n)
