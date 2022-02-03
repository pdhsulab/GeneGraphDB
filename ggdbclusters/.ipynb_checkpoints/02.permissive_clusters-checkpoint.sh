THREADS=32

mkdir -p ../../clusters/OUTPUT/permissive/tmp

mmseqs easy-linclust ../../clusters/INPUT/permissive/clu_perm_mmseqs_input.faa ../../clusters/OUTPUT/permissive/tmp/clu ../../clusters/OUTPUT/permissive/tmp/clu_tmp --threads ${THREADS} -e 0.001 --min-seq-id 0.3 -c 0.3 --cov-mode 0 --spaced-kmer-mode 0 --remove-tmp-files 1


# mmseqs createdb OUTPUT/stringent/stringent_rep_seq.faa OUTPUT/permissive/tmp/seqdb --dbtype 1 --shuffle 1 --createdb-mode 1 --write-lookup 0 --id-offset 0 --compressed 0 -v 3
# mmseqs linclust OUTPUT/permissive/tmp/seqdb OUTPUT/permissive/tmp/clu OUTPUT/permissive/tmp/clu_tmp --threads ${THREADS} -e 0.001 --min-seq-id 0.3 -c 0.9 --cov-mode 1 --spaced-kmer-mode 0 --remove-tmp-files 1
# mmseqs createtsv OUTPUT/permissive/tmp/seqdb OUTPUT/permissive/tmp/seqdb OUTPUT/permissive/tmp/clu OUTPUT/permissive/permissive_clusters.tsv --threads ${THREADS} -v 3
# mmseqs result2repseq OUTPUT/permissive/tmp/seqdb OUTPUT/permissive/tmp/clu OUTPUT/permissive/tmp/clu_rep --db-load-mode 0 --compressed 0 --threads ${THREADS} -v 3
# mmseqs result2flat OUTPUT/permissive/tmp/seqdb OUTPUT/permissive/tmp/seqdb OUTPUT/permissive/tmp/clu_rep OUTPUT/permissive/permissive_rep_seq.faa --use-fasta-header -v 3

