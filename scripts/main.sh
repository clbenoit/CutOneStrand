#!/usr/bin/env bash

## environment variables
#CONDAPATH=/home/ptngs/miniconda3/
CONDAPATH=/data/home/cbenoit3/miniconda3/
show_help() {
    echo -e ""
    echo -e "############################### HELP ####################################################\n"
    echo -e "                 CutOneStrand version = 1.0.0                        \n"
    echo "Usage: $0 [args...]" >&2
    echo -e "Available arguments : \n"
    echo "   -g, --gene,           gene to scan for positions to cut on one strand only, ex : RYR1"
    echo "   -c, --cas,            cas9 you want to use to cut your gene (Only spcas9 available on v1.0)" 
    echo "   -f, --frequence,      Minimal variant frequency in gnomAD v.3 population"
    echo "   -o, --output,        Output file name to store results in"
    echo "   -h, --help,          Show this help section"
    echo ""
    echo -e "This pipeline was developped at the CHU Grenoble Alpes"
    echo -e "Feel free to address any issue at : benoitclement.sand@gmail.com \n"
    echo -e "##########################################################################################\n"
    echo
    # echo some stuff here for the -a or --add-options 
    exit 1
}

date=$(date)
echo -e "\n"
echo -e "####### Starting CutOneStrand V1 (${date}) #######" 

PARAMS=""
while (( "$#" )); do
  case "$1" in
    -g|--gene)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        GENE=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -c|--cas)
      #sp_20 sa_21 cas12a
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
	#if [ $2 != "sp_20" ]  && [ $2 != "sa_21" ] && [ $2 != "cpf1_24" ];then
	#if [ $2 != "spcas9ngg" ]  && [ $2 != "cpf1" ];then
	if [ $2 != "spcas9ngg" ];then
		echo -e "FAILED : Unrecognized value parsed to -c|--cas argument. Available cas9 :"
		echo -e "spcas9ngg"
		#echo -e "sa_21"
		#echo -e "cpf1 \n"
		date=$(date)
		echo -e "Stopping CutOneStrand V1... (${date})\n " 
		exit 1
	else
        	CAS=$2
        shift 2
	fi
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1   
      fi
      ;; 
    -f|--frequence)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        FREQ=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -h|--help)
        show_help
        shift 2
      ;;
    -o|--output)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        OUTFILE="$( cd "$(dirname "$2")" && pwd )/$( basename $2)"
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    # *) # preserve positional arguments
    #  PARAMS="$PARAMS $1"
    #  shift
    #  ;;
  esac
done

if [ "$OUTFILE" == "" ]; then
   echo "ERROR: Please provide an output filename"
   exit 1
fi
if [ "$GENE" == "" ]; then
   echo "ERROR: Please provide a gene name"
   exit 1
fi
if [ "$FREQ" == "" ]; then
   echo "ERROR: Please provide a mutation frequency"
   exit 1
elif [[ ! $FREQ =~ ^-?[0-9]+$ ]] ;then
	echo -e "Error: Argument for $1 has to be an integer between 0 and 100 \nExiting analysis... \n" >&2
	exit 1
fi
if [ "$CAS" == "" ]; then
   echo "ERROR: Please provide a cas name"
   exit 1
fi

### Define parameters ####
WDIR=$(dirname $(dirname $(readlink -f "${BASH_SOURCE}" )))
echo -e "with parameters :"
echo -e "Working directory : ${WDIR}"
echo -e "Gene : ${GENE}"
echo -e "Cas9 : ${CAS}"
echo -e "Minimal frequence : ${FREQ} %"
echo -e "Output file : ${OUTFILE} \n"

echo -e "####### Checking software environment... #######"
mkdir -p ${WDIR}/tools
TMPDIR=/tmp/CutOneStrand
mkdir -p ${TMPDIR}

### Install crispor.tefor python3 version ####
if { conda env list | grep 'preprocessdata'; } >/dev/null 2>&1
then
	echo -e "Dependencies already available for data preprocessing steps"
else
  echo -e "Installing python dependencies for data preprocessing steps..."
  conda env create --name preprocessdata -f ${WDIR}/environments/preprocessdata.yml
  echo "DONE"
fi

if { conda env list | grep 'java'; } >/dev/null 2>&1
then
	echo -e "Java dependancies already available for data preprocessing steps"
else
  echo -e "Installing java dependencies for data preprocessing steps..."
  conda env create --name java -f ${WDIR}/environments/java.yml
  echo "DONE"
fi

if { conda env list | grep 'pyvcf'; } >/dev/null 2>&1
then
	echo "Python dependencies available for results formating steps"
else
	echo "Installing python dependencies for results formating steps ..."
	conda env create --name pyvcf -f ${WDIR}/environments/pyvcf.yml
	echo "DONE"
fi

cd ${WDIR}/tools

if [ -f ${WDIR}/tools/flashry/FlashFry-assembly-1.15.jar ]; then 
	echo -e "Flashry already installed"

else 
	echo -e "Installing Flashry..."
	mkdir -p ${WDIR}/tools/flashry; cd ${WDIR}/tools/flashry
	wget https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar
	echo "OK Flashry installed on version : "
fi
echo -e "SUCCESS : Software environments successfully set up. \n"


echo -e "####### Checking required annotations... #######"
CHR=$(grep ${GENE} ${WDIR}/annotations/hg38/chr_ens_names_matches.csv | cut -f1 -d ',')

mkdir -p ${WDIR}/sequences/hg38
cd ${WDIR}/sequences/hg38
source ${CONDAPATH}/etc/profile.d/conda.sh
conda activate preprocessdata

# download fasta genome
if [ ! -f ${WDIR}/sequences/hg38/${CHR}.fa ];then
	echo -e "Downloading chromosomic reference sequence..."
	wget --timestamping "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/${CHR}.fa.gz"
	gunzip ${WDIR}/sequences/hg38/${CHR}.fa.gz
	echo -e "Creating samtools index on chromosome fasta file..."
	samtools faidx ${WDIR}/sequences/hg38/${CHR}.fa
else
	echo -e "Chromosomic reference sequence already present on youyr computing system"
fi

# Download gencode genes annotations  # hg38
mkdir -p ${WDIR}/annotations/hg38/gtfs
cd ${WDIR}/annotations/hg38/gtfs
if [ ! -f ${WDIR}/annotations/hg38/gtfs/gencode.v43.annotation.gtf ];then
  echo -e "Downloading gencode annotation..."
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
  gunzip ${WDIR}/annotations/hg38/gtfs/gencode.v43.annotation.gtf.gz
else
	echo -e "Gencode annotation already present on your computing system"
fi

grep ${GENE} ${WDIR}/annotations/hg38/gtfs/gencode.v43.annotation.gtf > ${WDIR}/annotations/hg38/gtfs/${GENE}_gencode.v43.annotation.gtf

### Get all ${GENE} SNPS from GnomAD
# Download ${CHR} genomic vcf v3.1
mkdir -p ${WDIR}/SNPS/hg38
cd ${WDIR}/SNPS/hg38
if [ ! -f ${WDIR}/SNPS/hg38/gnomad.genomes.v3.1.2.sites.${CHR}.vcf.bgz ];then
	echo -e "Downloading gnomAD chromosome file..."
	cd ${WDIR}/SNPS/hg38/
	wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.${CHR}.vcf.bgz
else
	echo -e "GnomAD chromosome file already present on your computing system"
fi
echo -e "SUCCESS : Required annotations successfully set up.\n"

echo -e "####### Computing analysis... #######"
##################" Keep variants with AF between 0.3 and 0.8 and with the FILTER field annotated as PASS ##################

FREQ_BCF=$(jq -n ${FREQ}/100 | grep -o "^....")
if [ ! -f ${WDIR}/SNPS/hg38/gnomad.genomes.v3.1.2.sites.${CHR}_filtered.vcf.bgz ];then
	echo -e "Filtering out gnomAD variants below ${FREQ} % frequency..."
	bcftools view -f "PASS" --threads 8 -i "INFO/AF[0] > ${FREQ_BCF}" ${WDIR}/SNPS/hg38/gnomad.genomes.v3.1.2.sites.${CHR}.vcf.bgz -Oz -o ${WDIR}/SNPS/hg38/gnomad.genomes.v3.1.2.sites.${CHR}_filtered.vcf.bgz
fi

# get header
if [ ! -f ${WDIR}/SNPS/hg38/${GENE}_${FREQ}_SNPS_CANDIDATES_hg38.vcf ];then
	echo -e "Intersecting gnomAD positions with freq > ${FREQ}% and ${GENE} annotations..."
	bcftools view -h ${WDIR}/SNPS/hg38/gnomad.genomes.v3.1.2.sites.${CHR}_filtered.vcf.bgz > ${WDIR}/SNPS/hg38/${GENE}_${FREQ}_SNPS_CANDIDATES_hg38.vcf
	bedtools intersect -b ${WDIR}/annotations/hg38/gtfs/${GENE}_gencode.v43.annotation.gtf -a ${WDIR}/SNPS/hg38/gnomad.genomes.v3.1.2.sites.${CHR}_filtered.vcf.bgz -wa -u >> ${WDIR}/SNPS/hg38/${GENE}_${FREQ}_SNPS_CANDIDATES_hg38.vcf
	echo -e "SUCCESS\n"
fi

###### Get genomic contest of interesing SNPs #######
if [ ! -d ${WDIR}/tools/jvarkit/JVARKIT/ ];then
	echo -e "Installing JVARKIT..."
	mkdir -p ${WDIR}/tools/jvarkit
	cd ${WDIR}/tools/jvarkit 
	wget https://uncloud.univ-nantes.fr/index.php/s/4sL77oWR2BFzSBH/download/jvarkit.jar
	unzip -o jvarkit.jar
	echo -e "SUCCESS \n"
fi

echo -e "Creating picard dictionnary..."
conda activate java
rm ${WDIR}/sequences/hg38/${CHR}.dict
picard CreateSequenceDictionary -R ${WDIR}/sequences/hg38/${CHR}.fa -O ${WDIR}/sequences/hg38/${CHR}.dict --VERBOSITY DEBUG

echo -e "\n"
echo -e "Getting genomic context on candidates positions..."
java -jar ${WDIR}/tools/jvarkit/JVARKIT/jvarkit.jar biostar251649 -r ${WDIR}/sequences/hg38/${CHR}.fa ${WDIR}/SNPS/hg38/${GENE}_${FREQ}_SNPS_CANDIDATES_hg38.vcf -n 35 | \
java -jar ${WDIR}/tools/jvarkit/JVARKIT/jvarkit.jar bioalcidaejdk --nocode -F VCF -e 'stream().forEach(V->println(">"+V.getContig()+":"+V.getStart()+"\n"+V.getAttribute("SEQ5_35")+"["+V.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/"))+"]"+V.getAttribute("SEQ3_35")));' > ${WDIR}/SNPS/hg38/${GENE}_${FREQ}_SNPS_CANDIDATES_hg38_WTH_FLANKING_SEQUENCES.fasta

echo -e "\n"
echo -e "Detecting mutations which create or delete crispr coding sites..."
source ${CONDAPATH}/etc/profile.d/conda.sh
conda activate pyvcf
python3 ${WDIR}/scripts/selectSNPs.py --fasta ${WDIR}/SNPS/hg38/${GENE}_${FREQ}_SNPS_CANDIDATES_hg38_WTH_FLANKING_SEQUENCES.fasta --mutlist ${WDIR}/SNPS/hg38/Mutations_modifying_crispr_cutting_sites.tsv --candidates ${WDIR}/SNPS/hg38/${GENE}_${FREQ}_SNPS_CANDIDATES_hg38.vcf --cas ${CAS} --candidatesfasta ${WDIR}/SNPS/hg38/candidates_sites.fasta

##### FLASHRY #######
## Creates index for fasta file
echo -e "\n"
echo -e "Creating flashfry database..."
mkdir -p ${WDIR}/tools/flashry/db
java -jar ${WDIR}/tools/flashry/FlashFry-assembly-1.15.jar index -reference ${WDIR}/sequences/hg38/${CHR}.fa -tmpLocation /tmp --enzyme ${CAS} --database ${WDIR}/tools/flashry/db/${CHR}_${CAS}_database
echo -e "Discovering flashfry targets"
java -jar ${WDIR}/tools/flashry/FlashFry-assembly-1.15.jar discover --database ${WDIR}/tools/flashry/db/${CHR}_${CAS}_database --fasta ${WDIR}/SNPS/hg38/candidates_sites.fasta --output ${WDIR}/SNPS/hg38/${GENE}_${CAS}_SNPS_OFF_TARGETS_hg38.txt --flankingSequence 20

# Compute flashfry scores ##
echo -e "Computing flashfry scores..."
java -jar ${WDIR}/tools/flashry/FlashFry-assembly-1.15.jar score \
 --input ${WDIR}/SNPS/hg38/${GENE}_${CAS}_SNPS_OFF_TARGETS_hg38.txt \
 --output ${WDIR}/SNPS/hg38/${GENE}_${CAS}_SNPS_OFF_TARGETS_hg38_scores.tsv \
 --scoringMetrics doench2016cfd,hsu2013,dangerous,minot,moreno2015,doench2014ontarget \
 --database ${WDIR}/tools/flashry/db/${CHR}_${CAS}_database

python3 ${WDIR}/scripts/flashFryResultsSNPinfosMerging.py --outfile ${OUTFILE} --mutlist ${WDIR}/SNPS/hg38/Mutations_modifying_crispr_cutting_sites.tsv --flashfry ${WDIR}/SNPS/hg38/${GENE}_${CAS}_SNPS_OFF_TARGETS_hg38_scores.tsv

echo -e "\n"
echo -e "Final results stored in ${OUTFILE}"


