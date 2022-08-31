#!/bin/bash
# Script to run LePhare's zphota program
PARAPATH=$LEPHARE/sel-4hi-q/lephare_parameters/
PARANAME=dr10_test_inputpara.para


while getopts "l:i:o:g:f:" opt; do
  case $opt in
    l) LIB_STEM=${OPTARG}
      ;;
    i) INPUT_STEM=${OPTARG}
      ;;
    o) OUTPUT_STEM=${OPTARG}
      ;;
    g) GLB_CONTEXT=${OPTARG}
      ;;
    f) FORB_CONTEXT=${OPTARG}
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done
# TODO: Program the input invocations more elegantly
LIB_STEM=${LIB_STEM:-baseline_templates}
LIB_STEM=${LIB_STEM}_pointlike
INPUT_STEM=${INPUT_STEM:-baseline_input}
OUTPUT_STEM=${OUTPUT_STEM:-test}
GLB_CONTEXT=${GLB_CONTEXT:--1}
FORB_CONTEXT=${FORB_CONTEXT:--1}

LIB_PATH=$LEPHARE/data/lephare_files/templates/
INPUT_PATH=$LEPHARE/data/lephare_input/
OUTPUT_PATH=$LEPHARE/data/lephare_output/
OUT_FPATH=${OUTPUT_PATH}${OUTPUT_STEM}_pointlike

printf $OUT_FPATH

# Check whether there already is an output file:
if [ -f "${OUT_FPATH}.fits" ]
    then
        printf "$OUT_FPATH already exists. Continue? (y/n)\n"
        read -n 1 -r;
        echo  # move to new line
    if [[ $REPLY =~ ^[Yy]$ ]]
        then
            echo  # move to new line
        else
            exit 0
    fi
fi
printf "\nLib-stem: $LIB_STEM\nInput-stem: $INPUT_STEM\nOutput-stem: $OUTPUT_STEM\nGlobal context: $GLB_CONTEXT, Forbidden context: $FORB_CONTEXT\n\n"

## Code to run photometric redshifts (pointlike sample) -> Forb_context = 24576 to exclude the HSC bands
$LEPHAREDIR/source/zphota -c ${PARAPATH}${PARANAME} \
-CAT_IN ${INPUT_PATH}${INPUT_STEM}_pointlike.in \
-CAT_OUT $OUT_FPATH.out \
-PARA_OUT ${PARAPATH}outputpara.para \
-ZPHOTLIB ${LIB_STEM}_mag_lib,star_mag_lib \
-MAG_REF 7 \
-MAG_ABS -30,-20 \
-GLB_CONTEXT $GLB_CONTEXT -FORB_CONTEXT $FORB_CONTEXT \
-PDZ_OUT ${OUT_FPATH}

# Convert the output to a .fits file, considering which output cols were requested:
jystilts="java -jar /home/hslxrsrv3/p1wx150/jystilts.jar"
$jystilts $LEPHARE/sel-4hi-q/jystilts_scripts/rewrite_fits_header.py ${PARAPATH}outputpara.para $OUT_FPATH.out


# Analyze the output:
python3 $LEPHARE/sel-4hi-q/assess_photoz_run.py $OUT_FPATH.fits pointlike -v