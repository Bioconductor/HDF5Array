#!/bin/bash
#

set -e  # exit immediately if a simple command exits with a non-zero status

while true; do
	ps u -p $1
	sleep $2
done

