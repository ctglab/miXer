# syntax=docker/dockerfile:1

# Use a base image that supports Conda and SIF
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy files into the container
COPY . /app
#@TODO: remove the featureImportance.yml 
# Create Conda environments from YAML files
RUN conda env create -f conda_envs/linux/featureImportance.yml && \
    conda env create -f conda_envs/linux/HMM_conda_env.yml && \
    conda env create -f conda_envs/linux/mergedBam.yml && \
    conda env create -f conda_envs/linux/miXer_ml_conda_env.yml && \
    conda env create -f conda_envs/linux/mixerSingularity.yml && \
    conda env create -f conda_envs/linux/mixerPre.yml

# Make sure that jq is installed, used to load json parameters in .sh entrypoint files
# Copy the entry point scripts into the container
# Copy the entry point scripts into the container
# Set the main entry point script

# Since these files change often, docker will run again the later instructions
# Without being able to use the layer caching mechanism (thus, longer image creation times)
# Moving them after the env creation step
RUN apt-get update && \
    apt-get install -y jq

# Copy the entry point scripts into the container
RUN chmod +x /app/entrypoints/*.sh && \
  chmod +x /app/main_entrypoint.sh && \
  mkdir -p /app/tmp/