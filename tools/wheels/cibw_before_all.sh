#!/bin/bash
# Install dependencies before building wheels
# This runs once per platform, not per Python version

set -ex

PROJECT_DIR="$1"

if [ "$RUNNER_OS" = "Linux" ]; then
    # On manylinux, install build dependencies
    # manylinux_2_28 is based on AlmaLinux 8

    # Install basic build tools
    yum install -y gcc gcc-c++ make cmake wget tar bzip2 zlib-devel zip unzip

    # Install HDF5 - package name varies by arch
    # On aarch64, we may need to enable PowerTools/CRB repo or build from source
    yum install -y epel-release || true
    yum install -y hdf5-devel || {
        # Try enabling PowerTools repo for aarch64
        yum install -y dnf-plugins-core || true
        dnf config-manager --set-enabled powertools || dnf config-manager --set-enabled crb || true
        yum install -y hdf5-devel || {
            # Build HDF5 from source as last resort
            echo "Building HDF5 from source..."
            cd /tmp
            wget https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.5/hdf5-1.14.5.tar.gz
            tar xzf hdf5-1.14.5.tar.gz
            cd hdf5-1.14.5
            ./configure --prefix=/usr/local --enable-cxx
            make -j$(nproc)
            make install
            ldconfig
            cd /tmp
        }
    }

    # Install Boost headers (not Boost.Python, just headers)
    yum install -y boost-devel || true

    # Install CFITSIO (needed by photospline)
    yum install -y cfitsio-devel || yum install -y cfitsio || true

    # Install SuiteSparse (needed by photospline)
    yum install -y suitesparse-devel || true

    # Build and install photospline from source
    PHOTOSPLINE_VERSION="2.3.1"
    cd /tmp
    wget https://github.com/icecube/photospline/archive/refs/tags/v${PHOTOSPLINE_VERSION}.tar.gz
    tar xzf v${PHOTOSPLINE_VERSION}.tar.gz
    cd photospline-${PHOTOSPLINE_VERSION}

    # Patch CMakeLists.txt to fix cmake_minimum_required for newer CMake versions
    # photospline uses VERSION 3.1.0 which is too old for modern CMake (3.27+)
    sed -i 's/cmake_minimum_required *(VERSION 3\.1/cmake_minimum_required(VERSION 3.5/' CMakeLists.txt || true
    sed -i 's/cmake_policy(VERSION 3\.1/cmake_policy(VERSION 3.5/' CMakeLists.txt || true

    mkdir build && cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release
    make -j$(nproc)
    make install
    ldconfig

    # Verify photospline installation
    echo "Photospline installed:"
    photospline-config --version || echo "photospline-config not in PATH"
    ls -la /usr/local/bin/photospline-config || true

    # Add to PATH
    export PATH="/usr/local/bin:$PATH"

elif [ "$RUNNER_OS" = "macOS" ]; then
    # On macOS, use Homebrew for basic dependencies
    # Note: Don't install suite-sparse from brew - it has header issues with photospline
    brew install hdf5 cfitsio boost

    # Build photospline from source (avoid Homebrew version due to header compatibility issues)
    # The Homebrew photospline + suite-sparse combo has extern "C" issues with C++ complex headers
    PHOTOSPLINE_VERSION="2.3.1"
    cd /tmp
    curl -L -o photospline.tar.gz https://github.com/icecube/photospline/archive/refs/tags/v${PHOTOSPLINE_VERSION}.tar.gz
    tar xzf photospline.tar.gz
    cd photospline-${PHOTOSPLINE_VERSION}

    # Patch CMakeLists.txt to fix cmake_minimum_required for newer CMake versions
    # photospline uses VERSION 3.1.0 which is too old for modern CMake (3.27+)
    sed -i.bak 's/cmake_minimum_required *(VERSION 3\.1/cmake_minimum_required(VERSION 3.5/' CMakeLists.txt || true
    sed -i.bak 's/cmake_policy(VERSION 3\.1/cmake_policy(VERSION 3.5/' CMakeLists.txt || true

    mkdir build && cd build
    # Build without SuiteSparse to avoid the extern "C" header conflict
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17
    make -j$(sysctl -n hw.ncpu)
    sudo make install

    echo "Photospline installed from source"
    /usr/local/bin/photospline-config --version
fi

echo "Dependencies installed successfully"
