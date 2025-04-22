#!/bin/bash

# Default options
parallel=false
show_number=false

# Parse options
while [[ "$1" == -* ]]; do
  case "$1" in
    -p) parallel=true ;;
    -n) show_number=true ;;
    --) shift; break ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
  shift
done

# Expect N and command
count=$1
shift

if ! [[ "$count" =~ ^[0-9]+$ ]]; then
  echo "Usage: runntimes.sh [-p] [-n] <N> <command> [args...]"
  exit 1
fi

# Main loop
for ((i=1; i<=count; i++)); do
  if $parallel; then
    (
      $show_number && echo "Run #$i:"
      "$@"
    ) &
  else
    $show_number && echo "Run #$i:"
    "$@" || exit 1
  fi
done

if $parallel; then
  wait
fi
