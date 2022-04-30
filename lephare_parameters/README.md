# Create input data:
export STEM=newtest
assemble_cat -h


# Create the extinction file:
$LEPHAREDIR_200509/source/filter_extinc -f compiled_filter_file  >exctinction.filt

# Create the star magnitude library (results in 'star_mag_lib'):
$LEPHAREDIR/source/sedtolib -c $LEPHARE/data/lephare_files/inputpara.para -t S \
-STAR_LIB star_compiled_file \
-STAR_SED $LEPHARE/data/lephare_files/template_lists/STAR_MOD_ALL.list

$LEPHAREDIR/source/mag_gal -c $LEPHARE/data/lephare_files/inputpara.para -t S \
    -LIB_ASCII YES \
    -STAR_LIB_IN star_compiled_file \
    -STAR_LIB_OUT star_mag_lib

# Create the galaxy pointlike magnitude library (results in 'pointlike_mag_lib'):
export POINTLIKE_LIST=pointlike_suraj.list # _baseline_templates.list
export POINTLIKE_LIB=pointlike_suraj # _baseline_templates.list
$LEPHAREDIR/source/sedtolib -c $LEPHARE/data/lephare_files/inputpara.para -t G \
-GAL_SED $LEPHARE/data/lephare_files/templates/lists/${POINTLIKE_LIST} \
-GAL_LIB ${POINTLIKE_LIB}_compiled_file

$LEPHAREDIR/source/mag_gal -c $LEPHARE/data/lephare_files/inputpara.para -t G \
-LIB_ASCII YES \
-GAL_LIB_IN ${POINTLIKE_LIB}_compiled_file \
-GAL_LIB_OUT ${POINTLIKE_LIB}_mag_lib \
-EXTINC_LAW  SMC_prevot.dat -MOD_EXTINC  11,23 \
-EB_V 0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4 -EM_LINES NO

# Create the galaxy extended magnitude library (results in 'extended_mag_lib'):
export EXTENDED_LIST=extended_suraj.list # baseline_templates.list
export EXTENDED_LIB=extended_suraj # _baseline_templates.list
$LEPHAREDIR/source/sedtolib -c $LEPHARE/data/lephare_files/inputpara.para -t G \
-GAL_SED $LEPHARE/data/lephare_files/templates/lists/${EXTENDED_LIST} \
-GAL_LIB ${EXTENDED_LIB}_compiled_file

$LEPHAREDIR/source/mag_gal -c $LEPHARE/data/lephare_files/inputpara.para -t G \
-LIB_ASCII YES \
-GAL_LIB_IN ${EXTENDED_LIB}_compiled_file \
-GAL_LIB_OUT ${EXTENDED_LIB}_mag_lib \
-EM_LINES NO

# Unfortunately, the ASCII template files are directly stored in the directory the command is run from.
mv ${EXTENDED_LIB}_mag_lib $LEPHARE/data/lephare_files/templates/
mv ${POINTLIKE_LIB}_mag_lib $LEPHARE/data/lephare_files/templates/
mv $star_mag_lib $LEPHARE/data/lephare_files/templates/
# The template files may need to be shifted to $LEPHARE/data/lephare_files/templates/.
jystilts $LEPHARE/lephare_scripts/jystilts/rewrite_fits_header.py 0 $LEPHARE/data/lephare_files/templates/star_mag_lib.dat
jystilts $LEPHARE/lephare_scripts/jystilts/rewrite_fits_header.py 0 $LEPHARE/data/lephare_files/templates/${POINTLIKE_LIB}_mag_lib.dat
jystilts $LEPHARE/lephare_scripts/jystilts/rewrite_fits_header.py 0 $LEPHARE/data/lephare_files/templates/${EXTENDED_LIB}_mag_lib.dat


export STEM=new_test
export STEM_OUT=without_i_suraj
export POINTLIKECONTEXT=8191 # 32767# 40959 # 8188
### Code to run photometric redshifts (pointlike sample) -> Forb_context = 24576 to exclude the HSC bands
$LEPHAREDIR/source/zphota -c $LEPHARE/data/lephare_files/inputpara.para \
-CAT_IN ${LEPHARE}/data/lephare_input/${STEM}_pointlike.in \
-CAT_OUT ${LEPHARE}/data/lephare_output/${STEM_OUT}_pointlike.out \
-PARA_OUT $LEPHARE/data/lephare_files/outputpara.para \
-ZPHOTLIB ${POINTLIKE_LIB}_mag_lib,star_mag_lib \
-MAG_REF 7 \
-MAG_ABS -30,-20 \
-ERR_SCALE 0.2,0.2,0.05,0.05,0.05,0.1,0.1,0.1,0.1,0.2,0.2,0.3,0.4,0.05,0.05,0.05 \
-GLB_CONTEXT $POINTLIKECONTEXT

jystilts $LEPHARE/lephare_scripts/jystilts/rewrite_fits_header.py $LEPHARE/data/lephare_files/outputpara.para ${LEPHARE}/data/lephare_output/${STEM_OUT}_pointlike.out


export EXTENDEDCONTEXT=8191# 32767 # 40959 # 8188
### Code to run photometric redshifts (extended sample) -> Forb_context = 30720 to exclude W3, W4 and the HSC bands
$LEPHAREDIR/source/zphota -c $LEPHARE/data/lephare_files/inputpara.para \
-CAT_IN ${LEPHARE}/data/lephare_input/${STEM}_extended.in \
-CAT_OUT ${LEPHARE}/data/lephare_output/${STEM_OUT}_extended.out \
-PARA_OUT $LEPHARE/data/lephare_files/outputpara.para \
-ZPHOTLIB ${EXTENDED_LIB}_mag_lib,star_mag_lib \
-ERR_SCALE 0.2,0.2,0.05,0.05,0.05,0.1,0.1,0.1,0.1,0.2,0.2,0.3,0.4,0.05,0.05,0.05 \
-MAG_REF 7 \
-MAG_ABS -24,-8 \
-GLB_CONTEXT $EXTENDEDCONTEXT

jystilts $LEPHARE/lephare_scripts/jystilts/rewrite_fits_header.py $LEPHARE/data/lephare_files/outputpara.para ${LEPHARE}/data/lephare_output/${STEM_OUT}_extended.out

