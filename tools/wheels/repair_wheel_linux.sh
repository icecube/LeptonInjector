#!/bin/bash
# Repair Linux wheel: fix RPATH and then run auditwheel
# Usage: repair_wheel_linux.sh <wheel> <dest_dir>

set -ex

WHEEL="$1"
DEST_DIR="$2"

# Create temp directory
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Unpack the wheel
unzip -q "$WHEEL" -d "$TMPDIR"

# Find all .so files in the LeptonInjector directory
PACKAGE_DIR="$TMPDIR/LeptonInjector"

# Fix RPATH for all shared libraries
fix_rpath() {
    local BINARY="$1"
    if [ ! -f "$BINARY" ]; then
        return
    fi

    echo "Setting RPATH in: $BINARY"
    patchelf --set-rpath '$ORIGIN' "$BINARY" 2>/dev/null || true
}

# Fix all .so files in the package
for f in "$PACKAGE_DIR"/*.so "$PACKAGE_DIR"/*.so.*; do
    if [ -f "$f" ]; then
        fix_rpath "$f"
    fi
done

# Also check the lib directory if it exists
if [ -d "$TMPDIR/lib" ]; then
    for f in "$TMPDIR/lib"/*.so "$TMPDIR/lib"/*.so.*; do
        if [ -f "$f" ]; then
            fix_rpath "$f"
        fi
    done
fi

# Repack the wheel with original name (must match dist-info directory)
WHEEL_NAME=$(basename "$WHEEL")
cd "$TMPDIR"
zip -q -r "$TMPDIR/$WHEEL_NAME" ./*

# Run auditwheel on the fixed wheel
# Set LD_LIBRARY_PATH so auditwheel can find external libraries (photospline, etc.)
# and bundled libraries from the wheel
export LD_LIBRARY_PATH="/usr/local/lib:/usr/local/lib64:$PACKAGE_DIR:$TMPDIR/lib:${LD_LIBRARY_PATH:-}"
echo "LD_LIBRARY_PATH set to: $LD_LIBRARY_PATH"
echo "Contents of $PACKAGE_DIR:"
ls -la "$PACKAGE_DIR" || true
echo "Contents of $TMPDIR/lib (if exists):"
ls -la "$TMPDIR/lib" 2>/dev/null || echo "lib directory does not exist"
echo "Contents of /usr/local/lib (system libraries installed by before-all):"
ls -la /usr/local/lib/*.so* 2>/dev/null | head -10 || echo "No .so files in /usr/local/lib"

auditwheel repair -w "$DEST_DIR" "$TMPDIR/$WHEEL_NAME"
