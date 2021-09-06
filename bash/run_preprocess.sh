
## Init objects to vdjtools scripts

#############
me=$(whoami)
application_files=/Users/$me/Dropbox/aplastic_anemia_tcr
files=/Users/$me/Dropbox/lag3/

vdj=$application_files/applications/vdjtools-1.2.1/vdjtools-1.2.1.jar
vdjdb=$application_files/applications/vdjdb-1.1.5/vdjdb-1.1.5.jar
vdjdb_new=$application_files/data/selected_tcrb/databases/vdjdb_new


#############

cd $files/data/tcrseq/prior/
prior=$(ls -d "$PWD"/*);

cd $files/data/tcrseq/io_naive/
io_naive=$(ls -d "$PWD"/*);

cd $files/data/tcrseq/tx/
tx=$(ls -d "$PWD"/*);

cd $files/data/tcrseq/robert/
robert=$(ls -d "$PWD"/*);

cd /Users/$me/Dropbox/Emerson/data/vdjt/
emerson=$(ls -d "$PWD"/*);

cd $files/data/tcrseq/lag3_sorted/cd4/
lag3_cd4=$(ls -d "$PWD"/*);

cd $files/data/tcrseq/lag3_sorted/cd8/
lag3_cd8=$(ls -d "$PWD"/*);

###############
cd $files
###############
clear

java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $lag3_cd4 $files/data/tcrseq/lag3_sorted/cd4/
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $lag3_cd8 $files/data/tcrseq/lag3_sorted/cd8/

java -Xmx4G -jar $vdj FilterNonFunctional $lag3_cd4 $files/data/tcrseq/lag3_sorted/cd4/
java -Xmx4G -jar $vdj FilterNonFunctional $lag3_cd8 $files/data/tcrseq/lag3_sorted/cd8/


java -Xmx4G -jar $vdj FilterNonFunctional $tx $files/data/tcrseq/tx/new/


## Calc diversity stats for each sample, with downsampled, to 20k reads
lag3=$(echo $prior $io_naive)
unsort_meta=$(echo $io_naive $prior $emerson $robert)

java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa --downsample-to 20000 $robert $files/results/tcrb/diversity/20k_sampled_robert
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa --downsample-to 20000 $lag3 $files/results/tcrb/diversity/20k_sampled_lag3
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa --downsample-to 20000 $unsort_meta $files/results/tcrb/diversity/20k_sampled_total_meta
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa --downsample-to 20000 $unsort_meta $files/results/tcrb/diversity/20k_sampled_total_meta
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa --downsample-to 20000 $unsort_meta $files/results/tcrb/diversity/20k_sampled_total_meta

## Subsample to 20k
java -Xmx4G -jar $vdj DownSample --size 20000 $lag3 $files/data/tcrseq/downsampled/




## GLIPH
gliph=/Users/hru/Documents/Laaketieteen_tohtori/applications/gliph-master/gliph/bin/gliph-group-discovery.pl
$gliph --kmer_mindepth=3 --simdepth=1000 --tcr $files/tcrb_data/gliph/pr1_total_filtered_gliph.txt
$gliph --kmer_mindepth=3 --simdepth=1000 --tcr $files/tcrb_data/gliph/pr1_total_unfiltered_gliph.txt



## Pool all Emerson
java -Xmx120G -jar $vdj JoinSamples --intersect-type aa --times-detected 3 --compress $emerson /Users/$me/Dropbox/Emerson/results/public/emerson_public
