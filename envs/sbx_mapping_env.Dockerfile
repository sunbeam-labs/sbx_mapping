FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_mapping_env

COPY envs/sbx_mapping_env.yml ./

# Install environment
RUN conda env create --file sbx_mapping_env.yml --name sbx_mapping

ENV PATH="/opt/conda/envs/sbx_mapping/bin/:${PATH}"

# Make sure the environment is working
SHELL ["conda", "run", "-n", "sbx_mapping", "/bin/bash", "-c"]

# Run
CMD "bash"