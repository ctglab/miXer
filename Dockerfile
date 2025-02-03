# syntax=docker/dockerfile:1

# Use a base image that supports Conda and SIF
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy files into the container
COPY conda_envs/ /app/conda_envs/

# Create Conda environments from YAML files
RUN conda env create -f conda_envs/linux/featureImportance.yml && \
    conda env create -f conda_envs/linux/HMM_conda_env.yml && \
    conda env create -f conda_envs/linux/mergedBam.yml && \
    conda env create -f conda_envs/linux/miXer_ml_conda_env.yml && \
    conda env create -f conda_envs/linux/mixerSingularity.yml && \
    conda env create -f conda_envs/linux/mixerPre.yml

# Make sure that jq is installed, used to load json parameters in .sh entrypoint files
RUN apt-get update && \
    apt-get install -y jq

# Copy the entry point scripts into the container
COPY entrypoints/ /app/entrypoints/
RUN chmod +x /app/entrypoints/*.sh

# Set the main entry point script
COPY main_entrypoint.sh /app/
RUN chmod +x /app/main_entrypoint.sh

# Since these files change often, docker will run again the later instructions
# Without being able to use the layer caching mechanism (thus, longer image creation times)
# Moving them after the env creation step

COPY preprocessingMixer/ /app/preprocessingMixer/
COPY processing/ /app/processing/
COPY tmp/ /app/tmp/
