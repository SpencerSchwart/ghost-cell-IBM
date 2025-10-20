#!/usr/bin/env bash
set -e

SRC="$1"

if [[ -z "$SRC" ]]; then
    echo "Usage: $0 <file.c> [-mpi] [-events] [...]" >&2
    exit 1
fi

CC="gcc"
QCCFLAGS="-source -autolink -disable-dimensions"
CFLAGS="-std=c99 -D_FORTIFY_SOURCE=2 -pipe -Wall -D_XOPEN_SOURCE=700 -O2"
LINKS="-L$BASILISK/gl -lglutils -lfb_tiny -lm"

for arg in "$@"; do
    case "$arg" in
        -mpi)
            CC="mpicc"
            QCCFLAGS+=" -D_MPI=1"
            ;;
        -openmp)
            QCCFLAGS+=" -fopenmp"
            CFLAGS+=" -fopenmp"
            ;;
        -display=1)
            QCCFLAGS+=" -DDISPLAY=1"
            ;;
        -display=-1)
            QCCFLAGS+=" -DDISPLAY=-1"
            ;;
        -events)
            QCCFLAGS+=" -events" 
            ;;
        -trace)
            QCCFLAGS+=" -DTRACE=2"
            ;;
        -debug)
            CFLAGS+=" -g -Wextra"
            ;;
    esac
done

OUT="${SRC##*/}"        # remove any leading path
OUT="${OUT%.c}"         # remove .c extension
OUTD="${OUT%.c}"        # out directory
OUT="$OUTD/$OUTD"       # move executable to directory with same name

GRN='\e[32m'
YEL='\e[33m'
BOLD='\e[1m'
RESET='\e[0m'

echo -e "\n${BOLD}${YEL}qcc${RESET}: flags = $QCCFLAGS"
echo -e "${BOLD}${YEL}qcc${RESET}: compiling ${BOLD}$SRC -> _$SRC ...\n"
qcc $QCCFLAGS "$SRC" 

mkdir -p "$OUTD"        # make new directory if one doesn't exists
echo -e "\n=================================================================================\n"

echo -e "${BOLD}${GRN}$CC${RESET}: flags = $CFLAGS"
echo -e "${BOLD}${GRN}$CC${RESET}: compiling ${BOLD}_$SRC -> $OUT ...\n"
$CC $CFLAGS "_$SRC" -o "$OUT" $LINKS

mv "_$SRC" "$OUTD"

echo -e "\nDone. Saved to $OUT"

