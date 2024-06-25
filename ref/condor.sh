#!/bin/bash

# Source the ROOT environment
source /opt/root61404/bin/thisroot.sh

# Ensure an index is provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 <index>"
  exit 1
fi

Index=$1

# Config file should be in the current directory and named as <index>.conf
Configuration="./${Index}.conf"

# Path to the executable
Execute="../bin/Coalescence"

# Check if the configuration file exists
if [ ! -f "$Configuration" ]; then
  echo "Configuration file ${Configuration} not found!"
  exit 1
fi

# Check if the executable exists and is executable
if [ ! -x "$Execute" ]; then
  echo "Executable ${Execute} not found or not executable!"
  exit 1
fi

# Execute the binary with the configuration file
$Execute "$Configuration"

