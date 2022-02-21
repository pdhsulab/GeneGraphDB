#!/bin/bash

#NEO4J_DIR="$GENEGRAPHDB_REPO_DIR/data/neo4j"
NEO4J_DIR="/home/john/GeneGraphDB/data/neo4j"

# https://neo4j.com/developer/docker-run-neo4j/

docker network create neo4j_network

docker run \
    --name testneo4j \
    -p7474:7474 -p7687:7687 \
    -d \
    --rm \
    -v $NEO4J_DIR/data:/data \
    -v $NEO4J_DIR/logs:/logs \
    -v $NEO4J_DIR/import:/var/lib/neo4j/import \
    -v $NEO4J_DIR/plugins:/plugins \
    --env NEO4J_AUTH=neo4j/test \
    --network "neo4j_network" \
    neo4j:latest

# network stuff
# https://stackoverflow.com/questions/42385977/accessing-a-docker-container-from-another-container
