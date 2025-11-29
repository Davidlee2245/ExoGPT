#!/bin/bash
# Download DGL CUDA 11.6 wheel file for Python 3.10

set -e

WHEEL_FILE="dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl"
WHEEL_URL="https://data.dgl.ai/wheels/repo.html"

echo "Downloading DGL wheel file: $WHEEL_FILE"
echo ""

# Try to download from DGL wheels repository
echo "Attempting to download from DGL wheels repository..."
if curl -s "$WHEEL_URL" | grep -q "$WHEEL_FILE"; then
    echo "Found wheel file in repository index."
    # Try direct download
    DIRECT_URL="https://data.dgl.ai/wheels/repo.html/$WHEEL_FILE"
    if curl -f -o "$WHEEL_FILE" "$DIRECT_URL" 2>/dev/null; then
        echo "✓ Successfully downloaded: $WHEEL_FILE"
        exit 0
    fi
fi

echo ""
echo "Direct download failed. Please download manually:"
echo ""
echo "Option 1: Visit DGL releases page:"
echo "  https://github.com/dmlc/dgl/releases"
echo "  Search for: dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl"
echo ""
echo "Option 2: Check DGL wheels repository:"
echo "  https://data.dgl.ai/wheels/repo.html"
echo "  Look for the wheel file matching: $WHEEL_FILE"
echo ""
echo "Option 3: Build from source (if wheel not available):"
echo "  See: https://github.com/dmlc/dgl"
echo ""
echo "After downloading, place the wheel file in the current directory and run:"
echo "  pip install $WHEEL_FILE"




