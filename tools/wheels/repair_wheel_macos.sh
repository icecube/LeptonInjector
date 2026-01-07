#!/bin/bash
# Repair macOS wheel using delocate
# Usage: repair_wheel_macos.sh <wheel> <dest_dir> <delocate_archs>

set -ex

WHEEL="$1"
DEST_DIR="$2"
DELOCATE_ARCHS="$3"

# Use delocate to bundle dependencies
# The --require-archs flag ensures we only bundle libraries for the target architecture
delocate-wheel --require-archs "$DELOCATE_ARCHS" -w "$DEST_DIR" -v "$WHEEL"
