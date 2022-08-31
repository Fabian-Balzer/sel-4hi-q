#!/bin/bash
# Script to run LePhare's MAG_GAL command for the star sources.
PARAPATH=$LEPHARE/sel-4hi-q/lephare_parameters/
PARANAME=dr10_test_inputpara.para

# Create the galaxy star magnitude library (results in 'star_${LIB_STEM}_mag_lib'):
LIB_STEM=${1:-baseline_templates}
LIB_STEM=${LIB_STEM}_star
LISTFILE=${LIB_STEM}.list
echo "Using ${LISTFILE} as the stem for the star magnitude library."


LISTPATH=$LEPHARE/sel-4hi-q/lephare_parameters/template_lists/

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
            echo "Proceeding..." # move to new line
        else
            exit 0
    fi
fi

# Create the star magnitude library (results in 'star_mag_lib'):
$LEPHAREDIR/source/sedtolib -c ${PARAPATH}${PARANAME} -t S \
-STAR_LIB ${LIB_STEM}_compiled_file \
-STAR_SED ${LISTPATH}${LISTFILE}

$LEPHAREDIR/source/mag_gal -c ${PARAPATH}${PARANAME} -t S \
    -LIB_ASCII YES \
    -STAR_LIB_IN ${LIB_STEM}_compiled_file \
    -STAR_LIB_OUT $MAGLIBNAME


# Move the ASCII file to the destined path:
mv ${MAGLIBNAME}.dat $LIBPATH

# Convert it to a .fits file:
jystilts="java -jar /home/hslxrsrv3/p1wx150/jystilts.jar"
$jystilts $LEPHARE/sel-4hi-q/jystilts_scripts/rewrite_fits_header.py 0 $MAGLIBFPATH
rm $MAGLIBFPATH
