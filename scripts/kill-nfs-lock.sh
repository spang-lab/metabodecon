#!/usr/bin/env bash
# kill-nfs-lock.sh
#
# Help recover from a stale R-package install lock left behind on NFS.
#
# By default (without --force) this script does NOT kill anything. It
# prints diagnostic info about each lock directory and instructions for
# the safe install workaround (`R CMD INSTALL --no-lock
# --no-staged-install .`), which writes directly to the final library
# path and avoids the lock-dir check entirely. This avoids killing R
# processes (such as long-running radian sessions) that have the
# package's shared library mmap'd.
#
# Pass --force to actually find and kill the processes holding the
# silly-rename .nfs* files open, then remove the LOCK directory.
#
# Usage:
#   ./scripts/kill-nfs-lock.sh [--force] [LOCK_DIR ...]
#
# Default: every directory matching ~/.Rpackages/00LOCK-*

set -euo pipefail

FORCE=0
LOCK_DIRS=()
for arg in "$@"; do
    case "$arg" in
        --force|-f) FORCE=1 ;;
        -h|--help)
            sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        --*)
            echo "Unknown option: $arg" >&2
            exit 2
            ;;
        *)
            LOCK_DIRS+=("$arg")
            ;;
    esac
done

if [[ ${#LOCK_DIRS[@]} -eq 0 ]]; then
    shopt -s nullglob
    LOCK_DIRS=("$HOME"/.Rpackages/00LOCK-*)
    shopt -u nullglob
fi

if [[ ${#LOCK_DIRS[@]} -eq 0 ]]; then
    echo "No lock directories found matching $HOME/.Rpackages/00LOCK-*"
    exit 0
fi

# Collect every PID with a /proc/<pid>/maps entry pointing into LOCK_DIR.
# These are the processes holding the .so files mmap'd; killing them
# releases the silly-rename .nfs* files so the directory can be removed.
find_pids() {
    local LOCK_DIR="$1"
    local out=""
    for maps in /proc/[0-9]*/maps; do
        local pid="${maps%/maps}"; pid="${pid#/proc/}"
        if grep -q "$LOCK_DIR" "$maps" 2>/dev/null; then
            out="$out $pid"
        fi
    done
    for fd_dir in /proc/[0-9]*/fd; do
        local pid="${fd_dir%/fd}"; pid="${pid#/proc/}"
        if ls -la "$fd_dir" 2>/dev/null | grep -q "$LOCK_DIR"; then
            out="$out $pid"
        fi
    done
    echo "$out" | tr ' ' '\n' | sort -u | tr '\n' ' ' | xargs
}

report_lock_dir() {
    local LOCK_DIR="$1"
    if [[ ! -d "$LOCK_DIR" ]]; then
        echo "Lock directory not found: $LOCK_DIR"
        return 0
    fi
    local PIDS
    PIDS=$(find_pids "$LOCK_DIR")
    echo "Lock directory: $LOCK_DIR"
    if [[ -z "$PIDS" ]]; then
        echo "  No live processes are holding files in this directory."
        echo "  You can simply: rm -rf '$LOCK_DIR'"
    else
        echo "  Processes holding files open in this directory:"
        for PID in $PIDS; do
            local cmd
            cmd=$(ps -p "$PID" -o comm= 2>/dev/null || echo 'unknown')
            echo "    PID $PID  ($cmd)"
        done
        echo
        echo "  To install without killing these processes (recommended),"
        echo "  run from the package root:"
        echo "    R CMD INSTALL --no-lock --no-staged-install ."
        echo
        echo "  '--no-lock' skips the lock-directory check; "
        echo "  '--no-staged-install' writes directly to the final library "
        echo "  path. Together they let you install while the above "
        echo "  processes keep the old shared library mmap'd."
        echo
        echo "  If you really need to free the lock directory, rerun this"
        echo "  script with --force to kill the listed PIDs."
    fi
}

kill_lock_dir() {
    local LOCK_DIR="$1"
    if [[ ! -d "$LOCK_DIR" ]]; then
        echo "Lock directory not found: $LOCK_DIR"
        return 0
    fi
    echo "Searching for processes locking files in: $LOCK_DIR"
    local PIDS
    PIDS=$(find_pids "$LOCK_DIR")
    if [[ -z "$PIDS" ]]; then
        echo "No processes found holding files open in $LOCK_DIR"
    else
        echo "Found PIDs: $PIDS"
        for PID in $PIDS; do
            echo "  Killing PID $PID ($(ps -p "$PID" -o comm= 2>/dev/null || echo 'unknown'))"
            kill "$PID" 2>/dev/null || true
        done
        sleep 1
        local REMAINING
        REMAINING=$(find_pids "$LOCK_DIR")
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
}

if [[ $FORCE -eq 1 ]]; then
    for LOCK_DIR in "${LOCK_DIRS[@]}"; do
        kill_lock_dir "$LOCK_DIR"
    done
    echo "Done."
else
    for LOCK_DIR in "${LOCK_DIRS[@]}"; do
        report_lock_dir "$LOCK_DIR"
    done
fi
