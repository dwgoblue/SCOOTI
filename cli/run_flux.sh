#!/bin/bash
# Usage: ./example/run_flux.sh example/demo_config.json

CONFIG_JSON="$1"

if [[ -z "$CONFIG_JSON" ]]; then
  echo "Usage: $0 path/to/config.json"
  exit 1
fi

if [[ ! -f "$CONFIG_JSON" ]]; then
  echo "Error: Config file not found: $CONFIG_JSON"
  exit 1
fi

# Resolve absolute path to config
CONFIG_ABS=$(realpath "$CONFIG_JSON")


matlab -nodisplay -nosplash -r "addpath(genpath(fullfile(pwd, 'metabolicModel'))); \
  config = jsondecode(fileread('$CONFIG_ABS')); \
  validate_config(config); \
  multiObj_CBM(config); \
  exit;"




#matlab -nodisplay -nosplash -r "try, \
#  addpath(genpath(fullfile(pwd, 'metabolicModel'))); \
#  config = jsondecode(fileread('$CONFIG_ABS')); \
#  validate_config(config); \
#  multiObj_CBM(config); \
#  catch ME, disp(getReport(ME)); exit(1); end; exit;"

