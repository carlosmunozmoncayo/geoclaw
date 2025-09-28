#!/bin/bash

#Sets some environment variables, run with source
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# echo $script_dir
export TEMPCLAW="$(cd ${script_dir}/.. && pwd)"
export TEMPGEOCLAW="$(cd ${script_dir} && pwd)"
export TEMPGEOLIB="$(cd ${script_dir}/src/1d_classic/shallow && pwd)"
export TEMPBOUSSLIB="$(cd ${script_dir}/src/1d_classic/hyp_bouss && pwd)"

export TEMPBOUSSLIB2D="$(cd ${script_dir}/src/2d/hyp_bouss && pwd)"
export TEMPGEOLIB2D="$(cd ${script_dir}/src/2d/shallow && pwd)"
