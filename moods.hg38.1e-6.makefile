SHELL=/bin/bash
THRESHOLD=1e-6
ASSEMBLY=hg38
REMOTE_TFBS_BED_GZ=/net/seq/data2/projects/ctrader/data/TF/increased_threshold/moods_1e-6/hg38.archetype_motifs.sorted.CT20240621.bed.gz
LOCAL_TFBS_PREFIX=${PWD}/moods.${ASSEMBLY}.${THRESHOLD}
LOCAL_TFBS_BED_GZ=${LOCAL_TFBS_PREFIX}.bed.gz
LOCAL_TFBS_BED=${LOCAL_TFBS_PREFIX}.bed
LOCAL_TFBS_WITH_SEQ_BED=${LOCAL_TFBS_PREFIX}.seq.bed
S3_DEST_URL=s3://areynolds-us-west-2/cd3plus/052524/tfbs-bam/

all:

copy:
	cp -f ${REMOTE_TFBS_BED_GZ} ${LOCAL_TFBS_BED_GZ}
	module add bedops && gunzip -c ${LOCAL_TFBS_BED_GZ} | sort-bed --max-mem 4G --tmpdir ${PWD} - > ${LOCAL_TFBS_BED}

addSeq:
	${PWD}/addSeq.py ${LOCAL_TFBS_BED} ${LOCAL_TFBS_WITH_SEQ_BED}

extract:
	${PWD}/extract_chroms.sh ${LOCAL_TFBS_WITH_SEQ_BED} ${THRESHOLD}

convert:
	${PWD}/convert_chroms.sh ${LOCAL_TFBS_WITH_SEQ_BED} ${THRESHOLD}

upload:
	aws s3 cp --dryrun . ${S3_DEST_URL} --recursive --exclude "*" --include "${LOCAL_TFBS_PREFIX}.seq.chr*.bam"
	aws s3 cp --dryrun . ${S3_DEST_URL} --recursive --exclude "*" --include "${LOCAL_TFBS_PREFIX}.seq.chr*.bam.bai"

upload-real:
	aws s3 cp . ${S3_DEST_URL} --recursive --exclude "*" --include "${LOCAL_TFBS_PREFIX}.seq.chr*.bam"
	aws s3 cp . ${S3_DEST_URL} --recursive --exclude "*" --include "${LOCAL_TFBS_PREFIX}.seq.chr*.bam.bai"
