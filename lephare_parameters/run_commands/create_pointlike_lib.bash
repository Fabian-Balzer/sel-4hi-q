#!/bin/bash
# Script to run LePhare's MAG_GAL command for the pointlike sources.
PARAPATH=$LEPHARE/lephare_scripts/lephare_parameters/

# Create the galaxy pointlike magnitude library (results in 'pointlike${LIB_STEM}_mag_lib'):
LIB_STEM=${1:-baseline_templates}
LIB_STEM=${LIB_STEM}_pointlike
LISTFILE=${LIB_STEM}.list
echo "Using ${LISTFILE} as the stem for the pointlike magnitude library."


LISTPATH=$LEPHARE/lephare_scripts/lephare_parameters/template_lists/

# Check whether the .list file even exists:
[ ! -f "${LISTPATH}${LISTFILE}" ] && {
    echo "${LISTFILE} was not found in ${LISTPATH}. Exiting."; exit 0
    }

LIBPATH=$LEPHARE/data/lephare_files/templates/
MAGLIBNAME=${LIB_STEM}_mag_lib
MAGLIBFPATH=${LIBPATH}${MAGLIBNAME}.dat

# Check whether there already is a .dat file:
if [ -f "${LIBPATH}${MAGLIBNAME}.fits" ]
    then
        printf "$MAGLIBFPATH already exists. Continue? (y/n)\n"
        read -n 1 -r;
        echo  # move to new line
    if [[ $REPLY =~ ^[Yy]$ ]]
        then
            echo  # move to new line
        else
            exit 0
    fi
fi

$LEPHAREDIR/source/sedtolib -c ${PARAPATH}inputpara.para -t G \
-GAL_SED ${LISTPATH}${LISTFILE} \
-GAL_LIB ${LIB_STEM}_compiled_file

$LEPHAREDIR/source/mag_gal -c ${PARAPATH}inputpara.para -t G \
-LIB_ASCII YES \
-GAL_LIB_IN ${LIB_STEM}_compiled_file \
-GAL_LIB_OUT $MAGLIBNAME \
-EXTINC_LAW  SMC_prevot.dat -MOD_EXTINC  11,23 \
-EB_V 0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4 -EM_LINES NO

# Move the ASCII file to the destined path:
mv ${MAGLIBNAME}.dat $LIBPATH

# Convert it to a .fits file:
jystilts="java -jar /home/hslxrsrv3/p1wx150/jystilts.jar"
$jystilts $LEPHARE/lephare_scripts/jystilts/rewrite_fits_header.py 0 $MAGLIBFPATH
rm $MAGLIBFPATH
