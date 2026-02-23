# GenomeViz Docker Image
# Genome Assembly Visualization Tool

FROM python:3.11-slim

LABEL maintainer="Aaron Thiel"
LABEL description="GenomeViz - Circular Genome Assembly Visualization Tool"
LABEL version="1.3.0"

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV QT_QPA_PLATFORM=offscreen

# Install system dependencies required for mappy (minimap2), and tools required by Nextflow
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    bash \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy source code
COPY src/ ./src/
COPY genomeViz.py .

# Make genomeViz.py executable and add /app to PATH so Nextflow processes can call it directly
RUN chmod +x genomeViz.py
ENV PATH="/app:${PATH}"
