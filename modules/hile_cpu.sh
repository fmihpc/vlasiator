#!/bin/bash

export HTTP_PROXY=http://www-cache.cs.helsinki.fi:3128
export HTTPS_PROXY=http://www-cache.cs.helsinki.fi:3128
export https_proxy=http://www-cache.cs.helsinki.fi:3128
export http_proxy=http://www-cache.cs.helsinki.fi:3128

module load papi
module load cray-pmi
module load libfabric

