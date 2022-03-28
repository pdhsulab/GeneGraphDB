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

alias cd_ggdb="cd $GENEGRAPHDB_REPO_DIR"

DOCKER_BASH_HISTORY="$GENEGRAPHDB_REPO_DIR/data/docker.bash_history"
touch $DOCKER_BASH_HISTORY

DOCKER_IMAGE="genegraphdb"
NEO4J_NETWORK_NAME="neo4j_network"
NEO4J_DIR="$GENEGRAPHDB_REPO_DIR/data/neo4j"

docker network create $NEO4J_NETWORK_NAME

# docker aliases
alias ggdb_docker_build="docker build -t $DOCKER_IMAGE $GENEGRAPHDB_REPO_DIR/developer_env"

alias ggdb_docker_build_mac_m1="docker build --platform linux/amd64 -t $DOCKER_IMAGE $GENEGRAPHDB_REPO_DIR/developer_env"

alias ggdb_docker_run="docker run -it --rm \
    -v $GENEGRAPHDB_REPO_DIR:/GeneGraphDB \
    -v $HOME/.config/gcloud:/root/.config/gcloud \
    -v $GENEGRAPHDB_REPO_DIR/data/docker.bash_history:/root/.bash_history \
    --network $NEO4J_NETWORK_NAME \
    $DOCKER_IMAGE"
    # -v /home/jupyter/GeneGraphDB/genegraph.db:/GeneGraphDB/data/genegraph.db \

alias ggdb_docker_jupyter="docker run -it --rm \
    --hostname localhost \
    -v $GENEGRAPHDB_REPO_DIR:/GeneGraphDB \
    -v $HOME/.config/gcloud:/root/.config/gcloud \
    -v $GENEGRAPHDB_REPO_DIR/data/docker.bash_history:/root/.bash_history \
    -p 0.0.0.0:8888:8888 \
    --network $NEO4J_NETWORK_NAME \
    $DOCKER_IMAGE \
    jupyter notebook \
        --port=8888 \
        --ip=0.0.0.0 \
        --allow-root \
        --no-browser \
        --NotebookApp.custom_display_url=http://localhost:8888"

# note the NEO4J_AUTH command determines:
#   * user name (default "neo4j"; cannot be changed)
#   * password ("test"; can be changed)
alias ggdb_neo4j_run="docker run -d --rm \
    --name testneo4j \
    -p 0.0.0.0:7474:7474 \
    -p 0.0.0.0:7687:7687 \
    -v $NEO4J_DIR/data:/data \
    -v $NEO4J_DIR/logs:/logs \
    -v $NEO4J_DIR/import:/var/lib/neo4j/import \
    -v $NEO4J_DIR/plugins:/plugins \
    -v $GENEGRAPHDB_REPO_DIR/data/csv_exports:/var/lib/neo4j/import/csv_exports \
    --env NEO4J_AUTH=neo4j/test \
    --network $NEO4J_NETWORK_NAME \
    neo4j:latest"
