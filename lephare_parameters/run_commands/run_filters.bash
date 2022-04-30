#!/bin/bash
# Script to run LePhare's filter command and provide extinction information if desired.
PARAPATH=$LEPHARE/lephare_scripts/lephare_parameters/

# Create the filter file:
FILTFILENAME=${1:-compiled_filter_file}
$LEPHAREDIR/source/filter -c ${PARAPATH}inputpara.para >${PARAPATH}${FILTFILENAME}.filt

echo "A filter file has been written to '${PARAPATH}${FILTFILENAME}.'"

while getopts ":e" opt; do
  case $opt in
    e) # Create the extinction file:
        EXTINCTIONFILENAME=extinction_info_file.filt
        $LEPHAREDIR_200509/source/filter_extinc -f ${FILTFILENAME}  >${PARAPATH}${EXTINCTIONFILENAME}

        echo "An extinction info file has been written to '${PARAPATH}${EXTINCTIONFILENAME}.'"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done

