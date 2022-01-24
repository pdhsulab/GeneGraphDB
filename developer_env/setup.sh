#!/bin/bash

# Get the current dir.
if [ -n "$BASH_VERSION" ]; then
    DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
elif [ -n "$ZSH_VERSION" ]; then
    DIR=${0:a:h}  # https://unix.stackexchange.com/a/115431
else
	echo "Error: Unknown shell; cannot determine path to GeneGraphDB local repository"
fi
export GENEGRAPHDB_REPO_DIR="$(dirname $DIR)"

git --git-dir="$GENEGRAPHDB_REPO_DIR/.git" config --local core.autocrlf input
git --git-dir="$GENEGRAPHDB_REPO_DIR/.git" config filter.prepare_notebook_for_repository.clean 'developer_env/prepare_notebook_for_repository.py'

alias cd_genegraphdb="cd $GENEGRAPHDB_REPO_DIR"

DOCKER_BASH_HISTORY="$GENEGRAPHDB_REPO_DIR/data/docker.bash_history"
touch $DOCKER_BASH_HISTORY

DOCKER_IMAGE="genegraphdb"

# docker aliases
alias genegraphdb_docker_build="docker build -t $DOCKER_IMAGE $GENEGRAPHDB_REPO_DIR/developer_env"

alias genegraphdb_docker_run="docker run -it --rm \
    -v $GENEGRAPHDB_REPO_DIR:/GeneGraphDB \
    -v $HOME/.config/gcloud:/root/.config/gcloud \
    -v $GENEGRAPHDB_REPO_DIR/data/docker.bash_history:/root/.bash_history \
    -v /home/jupyter/GeneGraphDB/genegraph.db:/GeneGraphDB/data/genegraph.db \
    $DOCKER_IMAGE"

alias genegraphdb_docker_jupyter="docker run -it --rm \
    --hostname localhost \
    -v $GENEGRAPHDB_REPO_DIR:/GeneGraphDB \
    -v $HOME/.config/gcloud:/root/.config/gcloud \
    -v $GENEGRAPHDB_REPO_DIR/data/docker.bash_history:/root/.bash_history \
    -v /home/jupyter/GeneGraphDB/genegraph.db:/GeneGraphDB/data/genegraph.db \
    -p 0.0.0.0:8888:8888 \
    $DOCKER_IMAGE \
    jupyter notebook \
        --port=8888 \
        --ip=0.0.0.0 \
        --allow-root \
        --no-browser \
        --NotebookApp.custom_display_url=http://localhost:8888"