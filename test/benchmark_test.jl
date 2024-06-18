# SET OF BENCHMARK TESTS TO EVALUATE PERFORMANCES #
# =============================================== #

# Initialize database  - new way
using MAGEMin_C
using BenchmarkTools

data        =   Initialize_MAGEMin("ig", verbose=false; solver=2);
test        =   0         #KLB1
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   12.0
T           =   1100.0
@benchmark out         =   point_wise_minimization(P,T, data) seconds=10

# ================================= #
#           BENCHMARK RESULTS       #
# ================================= #
# Platform Info:
#   OS: Linux (x86_64-linux-gnu)
#   CPU: 12 × 11th Gen Intel(R) Core(TM) i5-11400H @ 2.70GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-15.0.7 (ORCJIT, tigerlake)
#   Threads: 1 on 12 virtual cores
# ================================= #
# BenchmarkTools.Trial: 73 samples with 1 evaluation.
#  Range (min … max):  65.803 ms … 71.038 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     69.352 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   69.282 ms ±  1.278 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================= #
# BenchmarkTools.Trial: 72 samples with 1 evaluation.
#  Range (min … max):  67.947 ms …  71.192 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     69.718 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   69.659 ms ± 602.103 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
# =========== NO NLopt_op ========= #
# julia> @benchmark out         =   point_wise_minimization(P,T, data)
# BenchmarkTools.Trial: 71 samples with 1 evaluation.
#  Range (min … max):  69.490 ms … 81.145 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     70.756 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   71.211 ms ±  1.592 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================ #



data        =   Initialize_MAGEMin("ig", verbose=false; solver=2);
test        =   6         #Hydrated MORB
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   4.0
T           =   450.0
@benchmark out         =   point_wise_minimization(P,T, data) seconds=10
# ================================= #
#           BENCHMARK RESULTS       #
# ================================= #
# Platform Info:
#   OS: Linux (x86_64-linux-gnu)
#   CPU: 12 × 11th Gen Intel(R) Core(TM) i5-11400H @ 2.70GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-15.0.7 (ORCJIT, tigerlake)
#   Threads: 1 on 12 virtual cores
# ================================= #
# BenchmarkTools.Trial: 144 samples with 1 evaluation.
#  Range (min … max):  66.900 ms … 83.233 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     69.018 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   69.575 ms ±  1.945 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================= #
# BenchmarkTools.Trial: 72 samples with 1 evaluation.
#  Range (min … max):  67.050 ms … 85.772 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     69.194 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   69.535 ms ±  2.256 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%
# =========== NO NLopt_op ========= #
# BenchmarkTools.Trial: 71 samples with 1 evaluation.
#  Range (min … max):  69.830 ms …  72.420 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     70.731 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   70.808 ms ± 457.227 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================ #



data        =   Initialize_MAGEMin("mb", verbose=false; solver=2);
test        =   0           #Amphibolite
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   4.0
T           =   450.0
@benchmark out         =   point_wise_minimization(P,T, data) seconds=10
# ================================= #
#           BENCHMARK RESULTS       #
# ================================= #
# Platform Info:
#   OS: Linux (x86_64-linux-gnu)
#   CPU: 12 × 11th Gen Intel(R) Core(TM) i5-11400H @ 2.70GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-15.0.7 (ORCJIT, tigerlake)
#   Threads: 1 on 12 virtual cores
# ================================= #
# BenchmarkTools.Trial: 198 samples with 1 evaluation.
# Range (min … max):  48.582 ms …  57.317 ms  ┊ GC (min … max): 0.00% … 0.00%
# Time  (median):     50.559 ms               ┊ GC (median):    0.00%
# Time  (mean ± σ):   50.723 ms ± 848.232 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================= #
# BenchmarkTools.Trial: 98 samples with 1 evaluation.
#  Range (min … max):  49.404 ms … 55.558 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     50.979 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   51.281 ms ±  1.153 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%
# =========== NO NLopt_op ========= #
# BenchmarkTools.Trial: 96 samples with 1 evaluation.
#  Range (min … max):  51.014 ms …  54.847 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     52.276 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   52.393 ms ± 797.302 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================ #




data        =   Initialize_MAGEMin("mb", verbose=false; solver=2);
test        =   0           #Amphibolite
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   4.0
T           =   850.0
@benchmark out         =   point_wise_minimization(P,T, data) seconds=10
# ================================= #
#           BENCHMARK RESULTS       #
# ================================= #
# Platform Info:
#   OS: Linux (x86_64-linux-gnu)
#   CPU: 12 × 11th Gen Intel(R) Core(TM) i5-11400H @ 2.70GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-15.0.7 (ORCJIT, tigerlake)
#   Threads: 1 on 12 virtual cores
# ================================= #
# BenchmarkTools.Trial: 152 samples with 1 evaluation.
#  Range (min … max):  63.868 ms … 68.964 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     65.947 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   66.074 ms ±  1.183 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================= #
# BenchmarkTools.Trial: 76 samples with 1 evaluation.
#  Range (min … max):  64.724 ms …  67.171 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     66.075 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   66.057 ms ± 521.906 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
# =========== NO NLopt_op ========= #
# BenchmarkTools.Trial: 74 samples with 1 evaluation.
#  Range (min … max):  66.121 ms … 75.342 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     67.113 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   67.638 ms ±  1.651 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================ #




data        =   Initialize_MAGEMin("mp", verbose=false; solver=2);
test        =   4         #Garnet migmatite
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   4.0
T           =   450.0
@benchmark out         =   point_wise_minimization(P,T, data) seconds=10
# ================================= #
#           BENCHMARK RESULTS       #
# ================================= #
# Platform Info:
#   OS: Linux (x86_64-linux-gnu)
#   CPU: 12 × 11th Gen Intel(R) Core(TM) i5-11400H @ 2.70GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-15.0.7 (ORCJIT, tigerlake)
#   Threads: 1 on 12 virtual cores
# ================================= #
# BenchmarkTools.Trial: 441 samples with 1 evaluation.
#  Range (min … max):  21.632 ms …  25.343 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     22.606 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   22.700 ms ± 494.144 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================= #
# BenchmarkTools.Trial: 217 samples with 1 evaluation.
#  Range (min … max):  21.898 ms …  25.564 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     22.970 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   23.044 ms ± 613.545 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
# =========== NO NLopt_op ========= #
# BenchmarkTools.Trial: 215 samples with 1 evaluation.
#  Range (min … max):  22.323 ms …  25.432 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     23.192 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   23.269 ms ± 471.992 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================ #




data        =   Initialize_MAGEMin("um", verbose=false; solver=2);
test        =   0         #serpentinite
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   24.0
T           =   450.0
@benchmark out         =   point_wise_minimization(P,T, data) seconds=10
# ================================= #
#           BENCHMARK RESULTS       #
# ================================= #
# Platform Info:
#   OS: Linux (x86_64-linux-gnu)
#   CPU: 12 × 11th Gen Intel(R) Core(TM) i5-11400H @ 2.70GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-15.0.7 (ORCJIT, tigerlake)
#   Threads: 1 on 12 virtual cores
# ================================= #
# BenchmarkTools.Trial: 1172 samples with 1 evaluation.
#  Range (min … max):  7.951 ms …  23.076 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     8.439 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   8.525 ms ± 522.466 μs  ┊ GC (mean ± σ):  0.02% ± 0.46%
# ================================= #
# BenchmarkTools.Trial: 594 samples with 1 evaluation.
#  Range (min … max):  7.974 ms …  10.230 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     8.384 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   8.408 ms ± 217.712 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
# =========== NO NLopt_op ========= #
# BenchmarkTools.Trial: 561 samples with 1 evaluation.
#  Range (min … max):  8.492 ms …  10.215 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     8.863 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   8.911 ms ± 215.344 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================ #
