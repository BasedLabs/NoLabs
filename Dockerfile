# Use an official PyTorch runtime as a parent image
FROM pytorch/pytorch:1.11.0-cuda11.3-cudnn8-devel

# Set the working directory in the container to /app
WORKDIR /app

# Add the current directory contents into the container at /app
ADD . /app

RUN apt-get -y update
RUN apt-get -y install git
# Install Flask and any other needed packages specified in requirements.txt
# RUN pip install --no-cache-dir -r requirements.txt
RUN install_deps.sh

# Install NVM
RUN wget -qO- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.5/install.sh | bash
RUN export NVM_DIR="$([ -z "${XDG_CONFIG_HOME-}" ] && printf %s "${HOME}/.nvm" || printf %s "${XDG_CONFIG_HOME}/nvm")" \
[ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"

RUN nvm install --lts
WORKDIR /app/src/server/frontend
RUN npm install
WORKDIR /app

# Make port 5000 available to the world outside this container
EXPOSE 5000
WORKDIR /app
ENV PYTHONPATH=/app
# Run app.py when the container launches
ENTRYPOINT ["dockerfile_entrypoint.sh"]