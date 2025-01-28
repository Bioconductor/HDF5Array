#!/bin/bash
#

set -e  # exit immediately if a simple command exits with a non-zero status

pid="$1"
logfile="$2"
interval="$3"
./ps_infinite_loop.sh "$pid" "$interval" >"$logfile" &
echo $!

