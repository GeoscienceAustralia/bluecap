#!/usr/bin/env bash

# Run the Bluecap test suite.
#===========================

# NB this example script assumes that:
# 1. you have installed pytest;
# 2. the "main.py" script is located in the directory above the present directory; 
# 3. the bluecap_test_baseline folder is located two directories above the present directory. 

# If this is not the case, make a copy of this file and change the environmental variables
# "BLUECAP_PATH" to point to the directory containing the main.py script
# and 
# "BLUECAP_BASELINE_PATH" to point to the test baseline directory.


# Environment variables used by the bluecap test scripts

# path to main.py script in the bluecap folder
BLUECAP_PATH=$(pwd)/../

# path to test baseline
BLUECAP_BASELINE_PATH=$BLUECAP_PATH/../../bluecap_test_baseline/


# Running the tests

export BLUECAP_PATH
export BLUECAP_BASELINE_PATH

pytest ./