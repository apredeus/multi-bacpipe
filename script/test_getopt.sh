#!/bin/bash 

set -euo pipefail

! getopt --test > /dev/null
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    >&2 echo "ERROR: \"getopt --test\" failed in this environment"
    >&2 echo "Please make sure you have the most up-to-date version of getopt!" 
    exit 1
fi

! PARSED=$(getopt --options=p:r: --longoptions=simple,multi --name "$0" -- "$@")

# read getopt’s output this way to handle the quoting right:
eval set -- "$PARSED"

REF=""
CPUS=""
SIMPLE=false 
MULTI=false

while true; do
    case "$1" in
        -r)
            REF="$2" 
            shift 2
            ;;
        -p)
            CPUS="$2"
            shift 2
            ;;
        --simple)
            SIMPLE=true
            shift
            ;;
        --multi)
            MULTI=true
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done

echo "DEBUG2: $@" 

echo "Cores: $CPUS Ref: $REF Simple: $SIMPLE Multi: $MULTI"


