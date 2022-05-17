docker run \
    --rm -it \
    -p 8086:8086 --net neo4j_network \
    --name local_bigtable  \
    bigtruedata/gcloud-bigtable-emulator start --host-port=0.0.0.0:8086
