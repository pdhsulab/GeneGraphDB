THREADS=32

mkdir -p ../../clusters/OUTPUT/tmp

mmseqs easy-linclust ../../clusters/INPUT/stringent/mmseqs2_testdb_input.faa ../../clusters/OUTPUT/stringent/tmp/clu ../../clusters/OUTPUT/stringent/tmp/clu_tmp --threads ${THREADS} -e 0.001 --min-seq-id 0.9 -c 0.9 --cov-mode 0 --spaced-kmer-mode 0 --remove-tmp-files 1

