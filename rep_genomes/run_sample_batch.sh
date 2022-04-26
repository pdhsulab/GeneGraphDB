THREADS=1
STEP_SIZE=1

let START_ROW=${1}*$STEP_SIZE+1
let END_ROW=${START_ROW}+${STEP_SIZE}

#echo $START_ROW
#echo $END_ROW

SAMPLES=($(cat samples.txt | awk -v start=${START_ROW} -v end=${END_ROW} '{if(NR>=start && NR < end) print $0}'))

for SAMPLE in "${SAMPLES[@]}"
do

echo $SAMPLE

## SINGLE THREADED ONLY
python functions.py download_genome genome=WORKDIR/${SAMPLE}.fna.gz sample=${SAMPLE} rerun=False || { echo 'download_genome failed' ; exit 1; }
#python functions.py pyfastx_genome genome=WORKDIR/${SAMPLE}.fna.gz genome_pyfastx=WORKDIR/${SAMPLE}.fna.gz.fxi sample=${SAMPLE} rerun=False || { echo 'pyfastx failed' ; exit 1; }

## SINGLE THREADED ONLY
#python functions.py genome_stats genome=WORKDIR/${SAMPLE}.fna.gz genome_stats=WORKDIR/${SAMPLE}.genome_stats.tsv contigs=WORKDIR/${SAMPLE}.contigs.tsv.gz sample=${SAMPLE} rerun=False || { echo 'genome_stats failed' ; exit 1; }

## MULTITHREAD
#python functions.py prodigal genome=WORKDIR/${SAMPLE}.fna.gz faa=WORKDIR/${SAMPLE}.prodigal.faa.gz gff=WORKDIR/${SAMPLE}.prodigal.gff.gz sample=${SAMPLE} rerun=False threads=$THREADS || { echo 'prodigal failed' ; exit 1; }
#rm WORKDIR/${SAMPLE}.prodigal.faa.gz
#python functions.py hmmsearch_crisprcastyper faa=WORKDIR/${SAMPLE}.prodigal.faa.gz gff=WORKDIR/${SAMPLE}.prodigal.gff.gz genome_stats=WORKDIR/${SAMPLE}.genome_stats.tsv hmmsearch_out=WORKDIR/${SAMPLE}.prodigal.hmmsearch.crisprcastyper.tsv.gz sample=${SAMPLE} rerun=False threads=$THREADS || { echo 'hmmsearch_crisprcastyper failed' ; exit 1; }
#python functions.py mash_sketch genome=WORKDIR/${SAMPLE}.fna.gz mash_out=WORKDIR/${SAMPLE}.fna.msh sample=${SAMPLE} rerun=False threads=$THREADS || { echo 'mash_sketch failed' ; exit 1; }
#python functions.py genomesearch_markers faa=WORKDIR/${SAMPLE}.prodigal.faa.gz genome_stats=WORKDIR/${SAMPLE}.genome_stats.tsv genomesearch_markers=WORKDIR/${SAMPLE}.genomesearch.markers.tsv.gz sample=${SAMPLE} rerun=False threads=$THREADS || { echo 'genomesearch_markers failed' ; exit 1; }
#python functions.py minced genome=WORKDIR/${SAMPLE}.fna.gz genome_stats=WORKDIR/${SAMPLE}.genome_stats.tsv minced=WORKDIR/${SAMPLE}.minced.gff.gz sample=${SAMPLE} threads=$THREADS rerun=False || { echo 'minced failed' ; exit 1; }
#python functions.py crisprcasfinder genome=WORKDIR/${SAMPLE}.fna.gz genome_stats=WORKDIR/${SAMPLE}.genome_stats.tsv crisprcasfinder=WORKDIR/${SAMPLE}.crisprcasfinder.gff.gz sample=${SAMPLE} threads=$THREADS rerun=False || { echo 'crisprcasfinder failed' ; exit 1; }
#python functions.py hmmsearch_isescan faa=WORKDIR/${SAMPLE}.prodigal.faa.gz gff=WORKDIR/${SAMPLE}.prodigal.gff.gz hmmsearch_isescan_out=WORKDIR/${SAMPLE}.prodigal.hmmsearch.domains.isescan.tsv.gz isescan_hmm_path=/home/mdurrant/data/isescan.hmm sample=${SAMPLE} rerun=False threads=$THREADS || { echo 'hmmsearch_isescan failed' ; exit 1; }
## SINGLE THREADED ONLY
#python functions.py find_cas_fusions hmmsearch_out=WORKDIR/${SAMPLE}.prodigal.hmmsearch.crisprcastyper.tsv.gz faa=WORKDIR/${SAMPLE}.prodigal.faa.gz gff=WORKDIR/${SAMPLE}.prodigal.gff.gz cas_fusions_bed=WORKDIR/${SAMPLE}.cas_fusions.bed.gz cas_fusions_summary=WORKDIR/${SAMPLE}.cas_fusions.tsv.gz crisprcastyper_domains=WORKDIR/${SAMPLE}.crisprcastyper.domains.tsv.gz crisprcastyper_hmm=/home/mdurrant/CRISPESTDB/02.ISOLATE_GENOME_WORKFLOW/data/crisprcastyper_profiles.hmm sample=${SAMPLE} rerun=False || { echo 'find_cas_fusions failed' ; exit 1; }
python functions.py makeblastdb genome=WORKDIR/${SAMPLE}.fna.gz blastdb_nhr=WORKDIR/${SAMPLE}.blastdb.nhr blastdb_nin=WORKDIR/${SAMPLE}.blastdb.nin blastdb_nog=WORKDIR/${SAMPLE}.blastdb.nog blastdb_nsd=WORKDIR/${SAMPLE}.blastdb.nsd blastdb_nsi=WORKDIR/${SAMPLE}.blastdb.nsi blastdb_nsq=WORKDIR/${SAMPLE}.blastdb.nsq sample=${SAMPLE} rerun=False || { echo 'makeblastdb failed' ; exit 1; }
#python functions.py self_targeting_spacers genome=WORKDIR/${SAMPLE}.fna.gz crisprcasfinder=WORKDIR/${SAMPLE}.crisprcasfinder.gff.gz minced=WORKDIR/${SAMPLE}.minced.gff.gz faa=WORKDIR/${SAMPLE}.prodigal.faa.gz gff=WORKDIR/${SAMPLE}.prodigal.gff.gz self_targeting_spacers=WORKDIR/${SAMPLE}.self_targeting_spacers.tsv.gz sample=${SAMPLE} rerun=False || { echo 'self_targeting_spacers failed' ; exit 1; }
#python functions.py make_sqldb sqldb=WORKDIR/${SAMPLE}.db sample=${SAMPLE} rerun=False || { echo 'make_sqldb failed' ; exit 1; }
#python functions.py crispr_neighbors genome=WORKDIR/${SAMPLE}.fna.gz hmmsearch_cct=WORKDIR/${SAMPLE}.prodigal.hmmsearch.crisprcastyper.tsv.gz minced=WORKDIR/${SAMPLE}.minced.gff.gz faa=WORKDIR/${SAMPLE}.prodigal.faa.gz gff=WORKDIR/${SAMPLE}.prodigal.gff.gz sqldb=WORKDIR/${SAMPLE}.db crispr_neighbors=WORKDIR/${SAMPLE}.crispr_neighbors.db sample=${SAMPLE} rerun=False || { echo 'crispr_neighbors failed' ; exit 1; }
#python functions.py blast_casdelta_flanks blastdb=WORKDIR/${SAMPLE}.blastdb genome=WORKDIR/${SAMPLE}.fna.gz blast_info=WORKDIR/${SAMPLE}.blast_casdelta_flanks.tsv regions_fna=WORKDIR/${SAMPLE}.blast_casdelta_flanks.fna sample=${SAMPLE} rerun=False || { echo 'blast_casdelta_flanks failed' ; exit 1; }
#python functions.py annotate_casdelta_flanks regions_fna=WORKDIR/${SAMPLE}.blast_casdelta_flanks.fna gff=WORKDIR/${SAMPLE}.blast_casdelta_flanks.gff sample=${SAMPLE} rerun=False || { echo 'annotate_casdelta_flanks failed' ; exit 1; }
#python functions.py cmsearch_omegarna genome=WORKDIR/${SAMPLE}.fna.gz cmsearch_out=WORKDIR/${SAMPLE}.cmsearch.omegarna.tsv sample=${SAMPLE} rerun=False threads=$THREADS || { echo 'cmsearch_omegarna failed' ; exit 1; }

rm -rf WORKDIR/${SAMPLE}.*

done
