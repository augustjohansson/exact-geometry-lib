FROM ubuntu:24.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ccache \
    cmake \
    libeigen3-dev \
    libcgal-dev \
    python3 \
    python3-pip \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

ENV PATH="/usr/lib/ccache:${PATH}"

# Install Python dependencies
RUN pip3 install --no-cache-dir --break-system-packages \
    pybind11[global] \
    pytest

WORKDIR /simpex

# Copy the source tree
COPY . .

# Configure and build
RUN cmake -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -Dpybind11_DIR=$(python3 -c "import pybind11; print(pybind11.get_cmake_dir())") \
    && cmake --build build --parallel $(nproc)

# Run C++ tests
RUN ctest --test-dir build --output-on-failure

# Run Python tests (simpex module is in build/python/)
RUN SIMPEX_BUILD_DIR=/simpex/build/python pytest python/tests -v
