#!/bin/bash

rm -rf double float

(
    cp -r prototypes double
    cd double
    sed -i 's/TEST_TYPE_PROTOTYPE/double/g' *
)

(
    cp -r prototypes float
    cd float
    sed -i 's/TEST_TYPE_PROTOTYPE/float/g' *
)