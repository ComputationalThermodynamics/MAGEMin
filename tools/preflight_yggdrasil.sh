#!/usr/bin/env bash
set -uo pipefail

usage() {
cat <<'EOF'
Usage: tools/preflight_yggdrasil.sh [options] [MAGEMin_C package dir]

Release gate to run before opening a Yggdrasil PR for MAGEMin_jll.

  1. Builds libMAGEMin inside Docker with the Yggdrasil flags (gcc -O3 -std=c99).
  2. Runs the MAGEMin_C test suite against that library, several times, with
     --check-bounds=yes and coverage (as julia-actions/julia-runtest does).
     Fails on any test failure, crash, OOM kill, or peak memory above --max-mem-mb.
       - linux/amd64: the platform where the v2.0.6 regression showed up.
         On an x86_64 host: full test/runtests.jl for every --julia version.
         On an Apple Silicon host it runs under emulation, which cannot run
         multi-threaded Julia or Julia >= 1.11 reliably, so only the
         single-threaded db_infos + serial suites run, with --x86-julia.
       - linux/arm64 (Apple Silicon host only): full test/runtests.jl, including
         the threaded suite, for every --julia version, natively.
  3. Builds the MAGEMin CLI with AddressSanitizer + UndefinedBehaviorSanitizer and
     runs every database over several test bulks and P-T points. Fails on any
     sanitizer report. New allocations are filled with a non-zero pattern, so reads
     of uninitialized heap memory (e.g. the v2.0.6 DEW n_w bug) fail every time.

Requires Docker. The MAGEMin_jll compat entry is relaxed in a temporary copy of the
package, so the gate can run before the new jll exists.

  MAGEMin_C package dir  Julia package to test (default: this repository)
  --c-src DIR            MAGEMin C source tree to build (default: this repository)
  --julia "V ..."        Julia versions for full native runs (default: "1.10 1")
  --x86-julia "V ..."    Julia versions for emulated x86_64 runs (default: "1.10")
  --runs N               test-suite repetitions per job (default: 2)
  --max-mem-mb N         peak memory limit for the test stage (default: 4000)
  --docker-mem SIZE      container memory limit (default: 7g)
  --skip-tests           skip stages 1-2
  --skip-sanitizers      skip stage 3
  --keep                 keep the work directory
  -h, --help             show this help
EOF
}

REPO=$(cd "$(dirname "$0")/.." && pwd)
C_SRC="$REPO"
PKG="$REPO"
JULIA_VERSIONS="1.10 1"
X86_JULIA="1.10"
RUNS=2
MAX_MEM_MB=4000
DOCKER_MEM=7g
DO_TESTS=1
DO_SAN=1
KEEP=0
case "$(uname -m)" in arm64|aarch64) NATIVE=linux/arm64 ;; *) NATIVE=linux/amd64 ;; esac

while [ $# -gt 0 ]; do
    case "$1" in
        --c-src)           C_SRC=$(cd "$2" && pwd); shift 2 ;;
        --julia)           JULIA_VERSIONS="$2"; shift 2 ;;
        --x86-julia)       X86_JULIA="$2"; shift 2 ;;
        --runs)            RUNS="$2"; shift 2 ;;
        --max-mem-mb)      MAX_MEM_MB="$2"; shift 2 ;;
        --docker-mem)      DOCKER_MEM="$2"; shift 2 ;;
        --skip-tests)      DO_TESTS=0; shift ;;
        --skip-sanitizers) DO_SAN=0; shift ;;
        --keep)            KEEP=1; shift ;;
        -h|--help)         usage; exit 0 ;;
        -*)                echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
        *)                 PKG=$(cd "$1" && pwd); shift ;;
    esac
done

command -v docker >/dev/null || { echo "docker not found" >&2; exit 2; }
docker info >/dev/null 2>&1 || { echo "docker daemon not reachable" >&2; exit 2; }
[ -f "$C_SRC/Makefile" ] && [ -d "$C_SRC/src" ] || { echo "not a MAGEMin C tree: $C_SRC" >&2; exit 2; }
[ -f "$PKG/Project.toml" ] && [ -f "$PKG/test/runtests.jl" ] || { echo "not a MAGEMin_C package: $PKG" >&2; exit 2; }

WORK=$(mktemp -d "${TMPDIR:-/tmp}/magemin_preflight.XXXXXX")
cleanup() { [ "$KEEP" = 1 ] || rm -rf "$WORK"; }
trap cleanup EXIT

mkdir -p "$WORK/csrc" "$WORK/pkg" "$WORK/out"
rsync -a --exclude .git --exclude '*.o' --exclude 'libMAGEMin.*' --exclude '/MAGEMin' --exclude '/MAGEMin_asan*' \
      --exclude '/output' --exclude '/src/saves' "$C_SRC/" "$WORK/csrc/"
rsync -a --exclude .git --exclude 'Manifest*.toml' --exclude 'libMAGEMin.*' --exclude '*.o' \
      --exclude '/output' "$PKG/" "$WORK/pkg/"
awk '/^\[/{sec=$0} !(sec=="[compat]" && $1=="MAGEMin_jll"){print}' "$WORK/pkg/Project.toml" > "$WORK/pkg/Project.toml.new" \
    && mv "$WORK/pkg/Project.toml.new" "$WORK/pkg/Project.toml"

FAIL=0
SUMMARY=()
note() { SUMMARY+=("$1"); echo "$1"; }

APT='apt-get -qq update >/dev/null 2>&1 && apt-get -qq install -y gcc make libnlopt-dev libopenblas-dev liblapacke-dev >/dev/null 2>&1'

JOBS=()
if [ "$NATIVE" = linux/amd64 ]; then
    for V in $JULIA_VERSIONS; do JOBS+=("linux/amd64|$V|full"); done
else
    for V in $X86_JULIA;      do JOBS+=("linux/amd64|$V|serial"); done
    for V in $JULIA_VERSIONS; do JOBS+=("linux/arm64|$V|full"); done
fi

if [ "$DO_TESTS" = 1 ]; then
    cat > "$WORK/build_lib.sh" <<EOF
$APT || exit 3
cp -R /csrc /b && cd /b
make -j"\$(nproc)" USE_MPI=0 CC=gcc CCFLAGS="-O3 -g -fPIC -std=c99" LIBS="-lm -lopenblas -llapacke -lnlopt" INC="" lib > /out/build_lib_\$ARCH.log 2>&1 || exit 4
cp libMAGEMin.dylib /lib_out/
EOF
    for PLAT in $(printf '%s\n' "${JOBS[@]}" | cut -d'|' -f1 | sort -u); do
        ARCH=${PLAT#linux/}
        echo "==> stage 1: building libMAGEMin for $PLAT (gcc -O3)"
        mkdir -p "$WORK/lib_$ARCH"
        if docker run --rm --platform "$PLAT" -e ARCH="$ARCH" -v "$WORK/csrc:/csrc:ro" -v "$WORK/out:/out" -v "$WORK/lib_$ARCH:/lib_out" \
                -v "$WORK/build_lib.sh:/build_lib.sh:ro" julia:1.10 bash /build_lib.sh; then
            note "PASS  build libMAGEMin ($PLAT, gcc -O3)"
        else
            note "FAIL  build libMAGEMin $PLAT (see $WORK/out/build_lib_$ARCH.log)"; FAIL=1
        fi
    done
fi

if [ "$DO_TESTS" = 1 ]; then
    { echo "$APT || exit 3"; cat <<'EOF'
mkdir -p /run && cp /lib_in/libMAGEMin.dylib /run/ && cd /run
export JULIA_PROJECT=/pkg
export JULIA_PKG_CONCURRENT_DOWNLOADS=1
ok=0
for a in 1 2 3; do
    julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); using MAGEMin_C' > /out/instantiate_$TAG.log 2>&1 \
        && julia --check-bounds=yes --code-coverage=user -e 'using Test, MAGEMin_C' >> /out/instantiate_$TAG.log 2>&1 \
        && { ok=1; break; }
done
[ $ok = 1 ] || { echo "instantiate failed"; exit 3; }
status=0
for i in $(seq 1 "$RUNS"); do
    peakf=/out/peak_${TAG}_$i.kb; echo 0 > $peakf
    ( p=0; while :; do s=$(awk '/^VmRSS/{s+=$2} END{print s+0}' /proc/[0-9]*/status 2>/dev/null); [ "$s" -gt "$p" ] && p=$s && echo $p > $peakf; sleep 0.5; done ) &
    sampler=$!
    if [ "$MODE" = full ]; then
        julia --check-bounds=yes --code-coverage=user /pkg/test/runtests.jl > /out/tests_${TAG}_$i.log 2>&1
    else
        julia --check-bounds=yes --code-coverage=user -e 'using Test
            @testset "db_infos" begin include("/pkg/test/test_db_infos.jl") end
            @testset "serial"   begin include("/pkg/test/tests.jl") end' > /out/tests_${TAG}_$i.log 2>&1
    fi
    rc=$?
    kill $sampler 2>/dev/null
    peak_mb=$(( $(cat $peakf) / 1024 ))
    used_local=$(grep -c "Using locally compiled version of libMAGEMin" /out/tests_${TAG}_$i.log)
    echo "RESULT run=$i rc=$rc peak_mb=$peak_mb local_lib=$used_local"
    [ "$rc" -ne 0 ] && status=1
    [ "$used_local" -eq 0 ] && status=1
    [ "$peak_mb" -gt "$MAX_MEM_MB" ] && status=1
done
exit $status
EOF
    } > "$WORK/run_tests.sh"
    for JOB in "${JOBS[@]}"; do
        IFS='|' read -r PLAT JV MODE <<<"$JOB"
        ARCH=${PLAT#linux/}
        TAG="${ARCH}_${JV}_${MODE}"
        if [ ! -f "$WORK/lib_$ARCH/libMAGEMin.dylib" ]; then note "FAIL  tests $TAG skipped: no library for $PLAT"; FAIL=1; continue; fi
        echo "==> stage 2: MAGEMin_C tests $TAG, $RUNS run(s)"
        docker run --rm --platform "$PLAT" -m "$DOCKER_MEM" -e TAG="$TAG" -e MODE="$MODE" -e RUNS="$RUNS" -e MAX_MEM_MB="$MAX_MEM_MB" \
            -v "$WORK/pkg:/pkg" -v "$WORK/lib_$ARCH:/lib_in:ro" -v "$WORK/out:/out" -v "$WORK/run_tests.sh:/run_tests.sh:ro" \
            -v "magemin_preflight_depot_${ARCH}_$JV:/root/.julia" "julia:$JV" bash /run_tests.sh > "$WORK/out/stage2_$TAG.txt" 2>&1
        rc=$?
        cat "$WORK/out/stage2_$TAG.txt"
        while read -r line; do
            r=$(sed -E 's/.*run=([0-9]+).*/\1/' <<<"$line"); trc=$(sed -E 's/.*rc=([0-9]+).*/\1/' <<<"$line")
            pk=$(sed -E 's/.*peak_mb=([0-9]+).*/\1/' <<<"$line"); ll=$(sed -E 's/.*local_lib=([0-9]+).*/\1/' <<<"$line")
            verdict=PASS; why=""
            [ "$trc" != 0 ] && verdict=FAIL && why="exit $trc$( [ "$trc" = 137 ] && echo ' (OOM kill)'); "
            [ "$ll" = 0 ] && verdict=FAIL && why="${why}local lib not loaded; "
            [ "$pk" -gt "$MAX_MEM_MB" ] && verdict=FAIL && why="${why}peak ${pk} MB > ${MAX_MEM_MB} MB; "
            note "$verdict  tests $TAG run $r  peak ${pk} MB  ${why}(log: out/tests_${TAG}_$r.log)"
        done < <(grep '^RESULT' "$WORK/out/stage2_$TAG.txt")
        if [ "$rc" -ne 0 ]; then FAIL=1; fi
        grep -q '^RESULT' "$WORK/out/stage2_$TAG.txt" || note "FAIL  tests $TAG did not run (see out/stage2_$TAG.txt)"
    done
fi

if [ "$DO_SAN" = 1 ]; then
    echo "==> stage 3: CLI under ASan + UBSan"
    cat > "$WORK/run_san.sh" <<EOF
$APT || exit 3
SAN="-fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=undefined"
cp -R /csrc /b && cd /b
make -j"\$(nproc)" USE_MPI=0 EXE_NAME=MAGEMin CC=gcc CCFLAGS="-O1 -g -fPIC -std=c99 \$SAN" LIBS="-lm -lopenblas -llapacke -lnlopt \$SAN" INC="" all > /out/build_san.log 2>&1 || exit 4
export ASAN_OPTIONS=detect_leaks=0:halt_on_error=1:malloc_fill_byte=190:max_malloc_fill_size=1048576
export UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1
mkdir -p /out/san
n=0; bad=0; skipped=0
run(){ tag=\$1; shift; n=\$((n+1)); timeout 600 ./MAGEMin "\$@" > /out/san/\$tag.log 2>&1; rc=\$?
       if grep -q "Unknown test" /out/san/\$tag.log; then skipped=\$((skipped+1)); return; fi
       if [ \$rc -ne 0 ] || grep -qE "ERROR: AddressSanitizer|runtime error:" /out/san/\$tag.log; then bad=\$((bad+1)); echo "SANFAIL \$tag rc=\$rc"; fi; }
for db in mp mb mbe ig igad igd um ume mtl mpe all po sb11 sb21 sb24 xMELTS pMELTS rMELTS; do
  for t in 0 1; do
    for pt in "2 500" "10 800" "25 1200"; do
      set -- \$pt
      run "\${db}_t\${t}_P\$1_T\$2" --Verb=0 --db=\$db --test=\$t --Pres=\$1 --Temp=\$2
    done
  done
done
run all_ds62   --Verb=0 --db=all --ds=62  --test=0 --Pres=8 --Temp=800
run mp_ds636   --Verb=0 --db=mp  --ds=636 --test=0 --Pres=8 --Temp=600
run mb_cpx1    --Verb=0 --db=mb  --mbCpx=1 --test=0 --Pres=8 --Temp=700
run mbe_cpx1   --Verb=0 --db=mbe --mbCpx=1 --test=0 --Pres=8 --Temp=700
echo "SANSUMMARY runs=\$n failed=\$bad skipped_unknown_test=\$skipped"
[ \$bad -eq 0 ]
EOF
    docker run --rm --platform "$NATIVE" -m "$DOCKER_MEM" -v "$WORK/csrc:/csrc:ro" -v "$WORK/out:/out" -v "$WORK/run_san.sh:/run_san.sh:ro" \
        julia:1.10 bash /run_san.sh > "$WORK/out/stage3.txt" 2>&1
    rc=$?
    grep -E '^SANFAIL|^SANSUMMARY' "$WORK/out/stage3.txt"
    if [ "$rc" -eq 0 ]; then
        note "PASS  sanitizers ($(grep '^SANSUMMARY' "$WORK/out/stage3.txt" | sed 's/SANSUMMARY //'))"
    else
        FAIL=1
        if grep -q '^SANSUMMARY' "$WORK/out/stage3.txt"; then
            note "FAIL  sanitizers ($(grep '^SANSUMMARY' "$WORK/out/stage3.txt" | sed 's/SANSUMMARY //'); logs: out/san/)"
            for f in $(grep '^SANFAIL' "$WORK/out/stage3.txt" | awk '{print $2}' | head -5); do
                echo "---- $f"; grep -m1 -A12 -E "ERROR: AddressSanitizer|runtime error:" "$WORK/out/san/$f.log" | head -14
            done
        else
            note "FAIL  sanitizer build (see $WORK/out/build_san.log)"
        fi
    fi
fi

echo
echo "================ preflight summary ================"
echo "C source : $C_SRC"
echo "package  : $PKG"
printf '%s\n' "${SUMMARY[@]}"
if [ "$FAIL" = 0 ]; then
    echo "RESULT: PASS - OK to open the Yggdrasil PR"
else
    KEEP=1
    echo "RESULT: FAIL - work directory kept at $WORK"
fi
exit $FAIL
