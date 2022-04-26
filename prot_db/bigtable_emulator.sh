docker run --rm -it -p 8086:8086 --net merantixnet --name local_bigtable  bigtruedata/gcloud-bigtable-emulator start --host-port=0.0.0.0:8086
