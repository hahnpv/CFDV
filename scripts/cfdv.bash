#!/bin/bash
PARENT_PATH=$( cd "$(dirname $(dirname "${BASH_SOURCE[0]}"))" ; pwd -P )
PATH=$PATH:$PARENT_PATH/src/build:$PARENT_PATH/scripts
alias residual="residual.py"
alias sensor="sensor.py"
