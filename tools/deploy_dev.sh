#!/bin/sh
#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMin
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : nriel[at]uni-mainz.de
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#
#
#   Rebuilds libMAGEMin, regenerates the Julia bindings and the database tables,
#   deploys both to the MAGEMin_C development package and to MAGEMinApp, and runs
#   the test suite.
#
#   Usage:   tools/deploy_dev.sh [MAGEMin_C dev dir] [MAGEMinApp dir]
#
set -e

REPO=$(cd "$(dirname "$0")/.." && pwd)
DEV=${1:-"$HOME/.julia/dev/MAGEMin_C"}
APP=${2:-"$HOME/seph/MAGEMinApp.jl_v1.7.2"}

cd "$REPO"

echo "==> 1/6 regenerating C bindings (gen/magemin_library.jl)"
julia --project=. gen/generator.jl

echo "==> 2/6 rebuilding libMAGEMin"
make clean
make USE_MPI=0 lib

echo "==> 3/6 regenerating julia/db_infos_generated.jl from the fresh library"
julia --project=. gen/generate_db_infos.jl

echo "==> 4/6 running the test suite"
julia --project=. test/runtests.jl

echo "==> 5/6 deploying to $DEV"
mkdir -p "$DEV"
for d in gen julia src test; do
    rsync -a --delete \
          --exclude '*.o' --exclude '*.o.tmp' --exclude '.DS_Store' --exclude '.vscode' \
          "$REPO/$d/" "$DEV/$d/"
done
cp "$REPO/Project.toml" "$DEV/Project.toml"
cp "$REPO/libMAGEMin.dylib" "$DEV/libMAGEMin.dylib"

echo "==> 6/6 deploying libMAGEMin.dylib to $APP"
if [ -d "$APP" ]; then
    cp "$REPO/libMAGEMin.dylib" "$APP/libMAGEMin.dylib"
else
    echo "    skipped: $APP does not exist"
fi

echo "done."
