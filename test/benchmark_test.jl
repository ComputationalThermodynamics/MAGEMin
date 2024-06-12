# SET OF BENCHMARK TESTS TO EVALUATE PERFORMANCES #
# =============================================== #

# Initialize database  - new way
using MAGEMin_C
using BenchmarkTools

data        =   Initialize_MAGEMin("ig", verbose=false);
test        =   0         #KLB1
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   12.0
T           =   1100.0
@benchmark out         =   point_wise_minimization(P,T, data)

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
# BenchmarkTools.Trial: 115 samples with 1 evaluation.
#  Range (min … max):  42.994 ms …  45.864 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     43.618 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   43.795 ms ± 529.447 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================ #




data        =   Initialize_MAGEMin("ig", verbose=false);
test        =   6         #Hydrated MORB
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   4.0
T           =   450.0
@benchmark out         =   point_wise_minimization(P,T, data)
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
# BenchmarkTools.Trial: 41 samples with 1 evaluation.
#  Range (min … max):  121.237 ms … 127.681 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     124.073 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   124.533 ms ±   1.510 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================ #




data        =   Initialize_MAGEMin("mb", verbose=false);
test        =   0           #Amphibolite
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   4.0
T           =   450.0
@benchmark out         =   point_wise_minimization(P,T, data)
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
# BenchmarkTools.Trial: 53 samples with 1 evaluation.
#  Range (min … max):  93.409 ms … 104.757 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     95.159 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   95.608 ms ±   1.932 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================ #




data        =   Initialize_MAGEMin("mb", verbose=false);
test        =   0           #Amphibolite
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   4.0
T           =   850.0
@benchmark out         =   point_wise_minimization(P,T, data)
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
# BenchmarkTools.Trial: 71 samples with 1 evaluation.
#  Range (min … max):  67.973 ms … 78.388 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     70.239 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   70.633 ms ±  1.793 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================ #




data        =   Initialize_MAGEMin("mp", verbose=false);
test        =   4         #Garnet migmatite
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   4.0
T           =   450.0
@benchmark out         =   point_wise_minimization(P,T, data)
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
# BenchmarkTools.Trial: 300 samples with 1 evaluation.
#  Range (min … max):  15.707 ms …  17.794 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     16.648 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   16.664 ms ± 391.919 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================ #




data        =   Initialize_MAGEMin("um", verbose=false);
test        =   0         #serpentinite
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   24.0
T           =   500.0
@benchmark out         =   point_wise_minimization(P,T, data)
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
# BenchmarkTools.Trial: 701 samples with 1 evaluation.
#  Range (min … max):  6.541 ms …  10.026 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     7.044 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   7.119 ms ± 333.453 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
# ================================ #