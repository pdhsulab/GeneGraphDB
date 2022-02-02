#!/bin/bash

NEO4J_DIR="/Users/john/zoo_ventures/arc_institute/GeneGraphDB/data"

# https://neo4j.com/developer/docker-run-neo4j/

docker run \
    --name testneo4j \
    -p7474:7474 -p7687:7687 \
    -d \
    --rm \
    -v $NEO4J_DIR/neo4j/data:/data \
    -v $NEO4J_DIR/neo4j/logs:/logs \
    -v $NEO4J_DIR/neo4j/import:/var/lib/neo4j/import \
    -v $NEO4J_DIR/neo4j/plugins:/plugins \
    --env NEO4J_AUTH=neo4j/test \
    --network "neo4j_network" \
    neo4j:latest

# network stuff
# https://stackoverflow.com/questions/42385977/accessing-a-docker-container-from-another-container
