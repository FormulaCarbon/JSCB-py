# Benchmarking `sig()` Against SciPy

import math
import random
import time
from statistics import mean

import numpy as np
import scipy.special

from src.jscb import sig


# =========================================================
# SCIPY VERSION
# =========================================================

def sig_scipy(rnn: float, djsmx: float) -> float:
    snn = rnn * math.log(2.0) * djsmx
    return scipy.special.gammainc(30.0, snn)


# =========================================================
# VECTORIZED SCIPY VERSION
# =========================================================

def sig_scipy_vectorized(rnn_array, djsmx_array):
    snn = rnn_array * math.log(2.0) * djsmx_array
    return scipy.special.gammainc(30.0, snn)


# =========================================================
# GENERATE TEST INPUTS
# =========================================================

N = 1_000_000

random.seed(0)

inputs = [
    (
        random.uniform(1.0, 1000.0),
        random.uniform(0.0, 5.0)
    )
    for _ in range(N)
]

print(f"Generated {N:,} test cases")
print()


# =========================================================
# BENCHMARK HELPER
# =========================================================

def benchmark_scalar(func, name):

    start = time.perf_counter()

    results = [func(rnn, djsmx) for rnn, djsmx in inputs]

    end = time.perf_counter()

    elapsed = end - start

    print(name)
    print(f"  Total time : {elapsed:.4f} sec")
    print(f"  Per call   : {elapsed / N * 1e6:.3f} µs")
    print()

    return results, elapsed


# =========================================================
# RUN SCALAR BENCHMARKS
# =========================================================

orig_results, orig_time = benchmark_scalar(
    sig,
    "Original sig()"
)

scipy_results, scipy_time = benchmark_scalar(
    sig_scipy,
    "SciPy sig()"
)


# =========================================================
# ACCURACY COMPARISON
# =========================================================

errors = [
    abs(a - b)
    for a, b in zip(orig_results, scipy_results)
]

print("Accuracy vs SciPy")
print(f"  Max error  : {max(errors)}")
print(f"  Mean error : {mean(errors)}")
print()


# =========================================================
# SPEEDUP
# =========================================================

print("Speedup")
print(f"  SciPy is {orig_time / scipy_time:.2f}x faster")
print()


# =========================================================
# VECTORIZED SCIPY BENCHMARK
# =========================================================

rnn_np = np.array([x[0] for x in inputs])
djs_np = np.array([x[1] for x in inputs])

start = time.perf_counter()

vec_results = sig_scipy_vectorized(rnn_np, djs_np)

end = time.perf_counter()

vec_time = end - start

print("Vectorized SciPy")
print(f"  Total time : {vec_time:.4f} sec")
print(f"  Per call   : {vec_time / N * 1e6:.6f} µs")
print()

print(f"Vectorized SciPy speedup vs original: {orig_time / vec_time:.2f}x")
