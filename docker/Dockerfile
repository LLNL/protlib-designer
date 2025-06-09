# Use the official Python 3.10 image, since protlib-designer requires Python >= 3.10
FROM python:3.10-slim-buster AS protlib_dev

# Install system packages for development:
# - git: for cloning the repository
# - build-essential: compilers (gcc/g++/make) in case any dependencies need compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    git build-essential \
 && apt-get remove -y build-essential \
 && apt-get autoremove -y \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

# (Removed package lists to keep the image slim)

# Set a working directory for the application code
WORKDIR /opt/protlib-designer

# Clone the protlib-designer repository from GitHub into the working directory
# Using git allows us to fetch the latest code for development
RUN git clone https://github.com/LLNL/protlib-designer.git .
# Now the directory /opt/protlib-designer contains the project source code

# Upgrade pip to the latest version and install protlib-designer with all its dependencies.
# We use ".[all,dev]" to include all optional features (PLM and inverse folding extras)
# and development tools (for testing and linting).
RUN pip install --upgrade pip && \
    pip install -e ".[all]"
# Explanation:
# - The "-e" flag installs in editable mode, so changes to the source reflect immediately (useful for development).
# - "[all,dev]" extra ensures all backend functionality and development/test dependencies are installed:
#    * "all" installs all additional scoring dependencies (e.g., language model and inverse folding models).
#    * "dev" installs tools like black, flake8, pytest for development and debugging.

# (Optional) Expose any ports if the application served a web service.
# protlib-designer is a CLI tool, so no ports need to be exposed by default.

# Set the default command to an interactive shell.
# This makes it easy to use the container for development (you can run protlib-designer or other tools inside).
# For example, once inside the container, you can run the sample command from the README:
#    protlib-designer ./example_data/trastuzumab_spm.csv 10
CMD ["bash"]
