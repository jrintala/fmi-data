#!/bin/sh

set -x

if [ -z "$HIMAN" ]; then
	export HIMAN="../../himan-bin/build/release/himan"
fi

rm -f RHO-KGM3*.grib

$HIMAN -d 5 -f density_ECMWF.json -t grib --no-cuda ecmwf_pressure_lvl_source.grib

grib_compare ./RHO-KGM3_pressure_850_rll_233_231_0_006.grib result_ecmwf_pressure_lvl.grib 

if [ $? -eq 0 ];then
  echo density/hirlam success on CPU!
else
  echo density/hirlam failed on CPU
  exit 1
fi

rm -f RHO-KGM3*.grib
