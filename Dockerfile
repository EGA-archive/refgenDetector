FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /usr/src/app

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        pkg-config \
        make \
        libssl-dev \
        libffi-dev \
        libpq-dev && \
    rm -rf /var/lib/apt/lists/*

# Copy repository contents
COPY . .

# Install Python dependencies
RUN python -m pip install --upgrade pip setuptools wheel && \
    python -m pip install -r src/requirements.txt

# Install refgenDetector from the local repository
RUN python -m pip install -e .

# Directory for input files
WORKDIR /data

# Keep the container running for repeated executions
ENTRYPOINT ["tail", "-f", "/dev/null"]
