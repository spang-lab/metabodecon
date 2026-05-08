#!/usr/bin/env bash
# kill-nfs-lock.sh
#
# Find and kill all processes holding open a .nfs* silly-rename lock file,
# then retry the rm -rf of the LOCK directory.
#
# Usage:
#   ./scripts/kill-nfs-lock.sh [LOCK_DIR]
#
# Default LOCK_DIR: ~/.Rpackages/00LOCK-metabodecon

set -euo pipefail

LOCK_DIR="${1:-$HOME/.Rpackages/00LOCK-metabodecon}"

if [[ ! -d "$LOCK_DIR" ]]; then
    echo "Lock directory not found: $LOCK_DIR"
    exit 0
fi

echo "Searching for processes locking files in: $LOCK_DIR"

# lsof and fuser only see open file descriptors, not memory-mapped files.
# NFS silly-rename (.nfs*) files are typically mmap'd by shared libraries
# loaded into running processes.  Scan /proc/*/maps to catch those too.
PIDS=""

# 1. /proc/*/maps — catches mmap'd .so files (the common NFS case)
for maps in /proc/[0-9]*/maps; do
    pid="${maps%/maps}"; pid="${pid#/proc/}"
    if grep -q "$LOCK_DIR" "$maps" 2>/dev/null; then
        PIDS="$PIDS $pid"
    fi
done

# 2. /proc/*/fd — catches open file descriptors (belt-and-suspenders)
for fd_dir in /proc/[0-9]*/fd; do
    pid="${fd_dir%/fd}"; pid="${pid#/proc/}"
    if ls -la "$fd_dir" 2>/dev/null | grep -q "$LOCK_DIR"; then
        PIDS="$PIDS $pid"
    fi
done

# Deduplicate
PIDS=$(echo "$PIDS" | tr ' ' '\n' | sort -u | tr '\n' ' ' | xargs)

if [[ -z "$PIDS" ]]; then
    echo "No processes found holding files open in $LOCK_DIR"
else
    echo "Found PIDs: $PIDS"
    for PID in $PIDS; do
        echo "  Killing PID $PID ($(ps -p "$PID" -o comm= 2>/dev/null || echo 'unknown'))"
        kill "$PID" 2>/dev/null || true
    done

    # Give processes a moment to release file handles.
    sleep 1

    # If any survived, force-kill.
    REMAINING=""
    for maps in /proc/[0-9]*/maps; do
        pid="${maps%/maps}"; pid="${pid#/proc/}"
        if grep -q "$LOCK_DIR" "$maps" 2>/dev/null; then
            REMAINING="$REMAINING $pid"
        fi
    done
    REMAINING=$(echo "$REMAINING" | tr ' ' '\n' | sort -u | tr '\n' ' ' | xargs)
    if [[ -n "$REMAINING" ]]; then
        echo "Force-killing remaining PIDs: $REMAINING"
        for PID in $REMAINING; do
            kill -9 "$PID" 2>/dev/null || true
        done
        sleep 1
    fi
fi

echo "Removing: $LOCK_DIR"
rm -rf "$LOCK_DIR"
echo "Done."
