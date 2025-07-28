# syntax=docker/dockerfile:1

# Use a base image that supports Conda and SIF
FROM ctglabcnr/excavator2

# Set the working directory
WORKDIR /app

# Pre-copy only the env files to leverage Docker caching
COPY conda_envs/ conda_envs/


# Create Conda environments from YAML files
RUN micromamba env create -f conda_envs/linux/HMM_conda_env.yml && \
    micromamba env create -f conda_envs/linux/miXer_ml_conda_env.yml && \
    micromamba env create -f conda_envs/linux/mixerPre.yml
# Copy files into the container
# Since these files change often, docker will run again the later instructions
# Without being able to use the layer caching mechanism (thus, longer image creation times)
# Moving them after the env creation step
COPY . /app

# Copy the entry point scripts into the container
RUN chmod +x /app/entrypoints/*.sh && \
  mkdir -p /app/tmp/

ENTRYPOINT ["/app/entrypoints/mixer.sh"]
