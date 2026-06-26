#!/usr/bin/env bash
set -e

SRC="$1"

if [[ -z "$SRC" ]]; then
    echo "Usage: $0 <file.c> [-mpi] [-events] [...]" >&2
    exit 1
fi

OUTD="$2"

check=${OUTD:0:1}
if [[ "$check" == "-" ]]; then
    OUT="${SRC##*/}"        # remove any leading path
    OUT="${OUT%.c}"         # remove .c extension
    OUTD="${OUT%.c}"        # out directory
fi
    
CC="gcc"
QCCFLAGS="-source -autolink -disable-dimensions"
CFLAGS="-std=c99 -D_GNU_SOURCE -D_FORTIFY_SOURCE=2 -pipe -Wall -D_XOPEN_SOURCE=700"
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
            LINKS+=" -L$BASILISK/wsServer -lws"
            ;;
        -display=-1)
            QCCFLAGS+=" -DDISPLAY=-1"
            LINKS+=" -L$BASILISK/wsServer -lws"
            ;;
        -events)
            QCCFLAGS+=" -events" 
            ;;
        -trace)
            QCCFLAGS+=" -DTRACE=2"
            ;;
        -debug)
            CFLAGS+=" -g -O0"
            ;;
        -O2)
            CFLAGS+=" -O2"
            ;;
        -g)
            CFLAGS+=" -g"
            ;;
        -strict) # Necessary to use with intel compiler
            CFLAGS+=" -fp-model=strict"
            ;;
    esac
done

OUT="${SRC##*/}"        # remove any leading path
OUT="${OUT%.c}"         # remove .c extension
OUT="$OUTD/$OUT"       # move executable to directory with same name

GRN='\e[32m'
YEL='\e[33m'
BOLD='\e[1m'
RESET='\e[0m'

echo -e "\n${BOLD}${YEL}qcc${RESET}: flags = $QCCFLAGS"
echo -e "${BOLD}${YEL}qcc${RESET}: compiling ${BOLD}$SRC -> _$SRC ...\n"
qcc $QCCFLAGS "$SRC" 

mkdir -p $OUTD
echo -e "\n=================================================================================\n"

echo -e "${BOLD}${GRN}$CC${RESET}: flags = $CFLAGS"
echo -e "${BOLD}${GRN}$CC${RESET}: compiling ${BOLD}_$SRC -> $OUT ...\n"
$CC $CFLAGS "_$SRC" -o "$OUT" $LINKS

mv "_$SRC" "$OUTD"

mkdir -p "$OUTD/plots"

if [[ "$SRC" == "impact-sphere.c"  || "$SRC" == "sphere.c" || "$SRC" == "impact-sphere3D.c" ]]; then
    mkdir -p "$OUTD/imgs"
    mkdir -p "$OUTD/dumps"
    mkdir -p "$OUTD/data"
    mkdir -p "$OUTD/outs"
elif [[ "$SRC" == "sessile-ibm.c" || "$SRC" == "sessile-ibm.c" || "$SRC" == "sessile-ibm.c" || "$SRC" == "sessile-ibm.c" ]]; then
    mkdir -p "$OUTD/outs"
    mkdir -p "$OUTD/results"
    mkdir -p "$OUTD/shapes"
elif [[ "$SRC" == "couette.c" ]]; then
    mkdir -p "$OUTD/data"
    mkdir -p "$OUTD/imgs"
    mkdir -p "$OUTD/profiles"
    mkdir -p "$OUTD/outs"
elif [[ "$SRC" == "sessile-ibm3D.c" ||  "$SRC" == "sessile-inclined-ibm3D.c" || "$SRC" == "sessile-sphere-ibm3D.c" || "$SRC" == "suspended-ibm3D.c"  || "$SRC" == "suspended-inclined-ibm3D.c" || "$SRC" == "sessile-oscillating-ibm3D.c" || "$SRC" == "sessile-oscillating-inclined-ibm3D.c" || "$SRC" == "suspended3D.c" || "$SRC" == "sessile-oscillating3D.c" ]]; then
    mkdir -p "$OUTD/results"
    mkdir -p "$OUTD/dumps"
    mkdir -p "$OUTD/imgs"
    mkdir -p "$OUTD/outs"
    mkdir -p "$OUTD/data"
fi

echo -e "\nDone. Saved to $OUT"
