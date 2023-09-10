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

from timeit import timeit

def fibonacci_exponential(n):
    """Computes the n-th number in the Fibonacci sequence, using the naive
    recursive algorithm implied by the definition of the Fibonacci sequence.

    Note: For simplicity, this ignores negative n, for which the basic form of
    the sequence is undefined.

    Runtime:
    A simplistic analysis of the runtime leads to O(2**n). The number of
    additions is exponential in n, and making the naive assumption that the
    values we are dealing with are small (e.g., 32-bit integers), additions are
    constant, so we have a runtime of O(2**n).

    For large enough n, the additions will be of integers not representable in the
    word-size of your system. The Fibonacci sequence grows rapidly, so addition
    shouldn't be assumed a constant time operation. The runtime to add integers
    is O(log m), where m is the value of the sum. Since the Fibonacci numbers grow
    exponentially, log(Fib_n) = O(n), thus the additions that take place have a
    runtime that grows linearly in n.

    Thus the runtime is actually O(n * 2**n).

    Keyword arguments:
    n - computes the n-th Fibonacci number.
    """
    if n <= 1:
        return n
    return fibonacci_exponential(n-1) + fibonacci_exponential(n-2)

def fibonacci_linear(n):
    """Computes the n-th number in the Fibonacci sequence, using a simple linear
    iteration over the sequence.

    Note: For simplicity, this ignores negative n, for which the basic form of
    the sequence is undefined.

    Runtime:
    Under the naive assumption that the values we are dealing with are small (e.g.,
    32-bit integers), this algorithm's runtime is O(n). It performs O(n) additions,
    and this naive assumption means that addition is constant time. Thus, runtime
    under this naive assumption is O(n).

    For large enough n, the additions will be of integers not representable in the
    word-size of your system. The Fibonacci sequence grows rapidly, so addition
    shouldn't be assumed a constant time operation. The runtime to add integers
    is O(log m), where m is the value of the sum. Since the Fibonacci numbers grow
    exponentially, log(Fib_n) = O(n), thus the additions that take place have a
    runtime that grows linearly in n.

    Thus, the runtime of this "linear" algorithm is not so linear, and is instead
    quadratic: O(n**2), since we are doing O(n) additions, and additions cost O(n).

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

    Runtime:
    Under the naive assumption that the values we are dealing with are small (e.g.,
    32-bit integers), this algorithm's runtime is O(lg n).

    It performs O(lg n) additions and mulitplications, and this naive assumption
    means that both operations are constant time. Thus, runtime under this naive
    assumption is O(lg n).

    For large enough n, the additions will be of integers not representable in the
    word-size of your system. The Fibonacci sequence grows rapidly, so addition
    shouldn't be assumed a constant time operation. The runtime to add integers
    is O(log m), where m is the value of the sum. Since the Fibonacci numbers grow
    exponentially, log(Fib_n) = O(n), thus the additions that take place have a
    runtime that grows linearly in n. 

    The runtime for the multiplications is a bit more complex. To multiply integers
    of magnitude m, the runtime is approximately O((log m)**1.58), assuming Karatsuba's
    multiplication algorithm, which is what Python uses for very large numbers. Due
    to the square and multiply approach, we never multiply integers larger than
    the sqrt(Fib_n). Thus the worst multiplication runs in O((log sqrt(Fib_n))**1.58)
    time. An upperbound for this is O(n**1.58) for the worst of the multiplications.

    Thus, the runtime of this "logarithmic" time algorithm is not so logarithmic. Its
    runtime is O(n**1.58 lg(n)). It performs O(lg(n)) additions that cost O(n), for
    total cost of additions: O(n log n). It performs O(lg n) multiplications that cost
    O(n**1.58) for total cost of multiplications: O(n**1.58 lg(n)).

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

def time_fibonacci_low_n(repetitions = 1, time_cutoff_exponential = 1.0):
    """Runs a timing experiment to deminstrate runtime of the three algorithms
    on low values of n, i.e., values of n such that the corresponding Fibonacci
    number is representable by a 32-bit signed int. Specifically, this considers
    values of n from 1 to 46.

    Keyword arguments:
    repetitions - this value is passed as number to timeit. Be careful what you use
        here as the times grow quickly. The default of 1 should be fine.
    time_cutoff_exponential - once the time for the exponential algorithm exceeds
        this value, it will just compare the other two.
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
    repetitions - this value is passed as number to timeit. Be careful what you use
        here as the times grow quickly. The default of 1 should be fine.
    time_cutoff_exponential - once the time for the exponential algorithm exceeds
        this value, it will just compare the other two.
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
            t_exponential = timeit(
                lambda : fibonacci_exponential(n),
                number = repetitions)
            template = template_all
            include_exponential = t_exponential < time_cutoff_exponential
        t_linear = timeit(
            lambda : fibonacci_linear(n),
            number=repetitions)
        t_logarithmic = timeit(
            lambda : fibonacci_logarithmic(n),
            number=repetitions)
        print(template.format(n, t_exponential, t_linear, t_logarithmic))
        n = n_update(n)
