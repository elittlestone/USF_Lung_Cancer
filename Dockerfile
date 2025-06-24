# Use the official Miniconda3 base image
FROM continuumio/miniconda3:latest

# Set working directory inside container
WORKDIR /workspace

# Install system dependencies (git needed to clone repos or for some packages)
RUN apt-get update && apt-get install -y git && rm -rf /var/lib/apt/lists/*

# Install Pixi in the base environment
RUN conda install -c conda-forge pixi -y

# Copy the Pixi environment definition file into the container
COPY pixi.toml /workspace/pixi.toml

# Use Pixi to create the Conda environment from the .toml file
RUN pixi env create -f pixi.toml

# Switch the default shell to activate the created environment automatically
SHELL ["conda", "run", "-n", "scenv", "/bin/bash", "-c"]

# Copy the rest of your project files into the container
COPY . /workspace

# Default command runs bash inside the activated conda environment
CMD ["bash"]
