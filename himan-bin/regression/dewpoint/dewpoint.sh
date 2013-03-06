#!/bin/sh

set -x

if [ -z "$HIMAN" ]; then
	export HIMAN="../../himan"
fi

rm -f himan*.grib

$HIMAN -d 3 -f dewpoint_arome.json -t grib source.grib

grib_compare result.grib dewpoint_arome.json.grib

if [ $? -eq 0 ];then
  echo dewpoint success on CPU!
else
  echo dewpoint failed on CPU
#  exit 1
fi

if [ "$CUDA_TOOLKIT_PATH" != "" ]; then

  mv dewpoint_arome.json.grib dewpoint_arome.json-CPU.grib

  $HIMAN -d 5 -f dewpoint_arome.json -t grib source.grib

  grib_compare -A 0.001 dewpoint_arome.json.grib dewpoint_arome.json-CPU.grib

  if [ $? -eq 0 ];then
    echo dewpoint success on GPU!
  else
    echo dewpoint failed on GPU
 #   exit 1
  fi

fi

