#!/bin/bash

pwd
ls -l

cd current
mkdir build
cd build
cmake ..
make -j

ls -l
