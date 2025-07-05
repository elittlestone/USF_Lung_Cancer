# Use the official Miniconda base image
FROM continuumio/miniconda3:latest

# Set the working directory inside the container
WORKDIR /workspace

# Install system dependencies (you can expand this list as needed)
RUN apt-get update && apt-get install -y \
    git \
    zsh \
    && rm -rf /var/lib/apt/lists/*

# Install Pixi via conda
RUN conda install -c conda-forge pixi -y


# Copy your pixi.toml file into the container and Install environment based on the pixi.toml file
COPY . /workspace
RUN pixi install

# Copy the user's local .zshrc into the container's home directory for zsh to pick it up
# Ensure it's sourced. We'll append to .zshrc to ensure it's sourced after the pixi activation.
RUN if [ -f /workspace/.zshrc ]; then \
    cat /workspace/.zshrc >> /root/.zshrc; \
    fi
    
# Default command drops you into a shell inside the pixi env
CMD ["pixi", "run", "zsh"]
