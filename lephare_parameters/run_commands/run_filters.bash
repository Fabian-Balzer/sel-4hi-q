#!/bin/bash
# Script to run LePhare's filter command and provide extinction information if desired.
PARAPATH=$LEPHARE/lephare_scripts/lephare_parameters/
PARANAME=dr10_test_inputpara.para

# Create the filter file:
FILTFILENAME=${1:-compiled_filter_file}
# Check whether this filter file already exists:
if [ -f "${PARAPATH}${FILTFILENAME}.filt" ]
    then
        printf "$FILTFILENAME.filt already exists. Continue? (y/n)\n"
        read -n 1 -r;
        echo  # move to new line
    if [[ $REPLY =~ ^[Yy]$ ]]
        then
            echo "Proceeding..." # move to new line
        else
            exit 0
    fi
fi
$LEPHAREDIR/source/filter -c ${PARAPATH}${PARANAME} \
-FILTER_REP $LEPHARE/data/lephare_files/filters \
>${PARAPATH}${FILTFILENAME}.filt

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

