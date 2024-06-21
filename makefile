SHELL=/bin/bash
THRESHOLD=1e-5
REMOTE_TFBS_BED_GZ=/net/seq/data2/projects/ctrader/data/TF/increased_threshold/moods_${THRESHOLD}/moods.combined.all.bed.gz
LOCAL_TFBS_PREFIX=${PWD}/moods.combined.all.${THRESHOLD}
LOCAL_TFBS_BED_GZ=${LOCAL_TFBS_PREFIX}.bed.gz
LOCAL_TFBS_BED=${LOCAL_TFBS_PREFIX}.bed
S3_DEST_URL=s3://areynolds-us-west-2/cd3plus/052524/tfbs-bam/

all:

copy:
	cp -f ${REMOTE_TFBS_BED_GZ} ${LOCAL_TFBS_BED_GZ}
	module add bedops && gunzip -c ${LOCAL_TFBS_BED_GZ} | sort-bed --max-mem 4G --tmpdir ${PWD} - > ${LOCAL_TFBS_BED}

extract:
	${PWD}/extract_chroms.sh ${LOCAL_TFBS_BED} ${THRESHOLD}

convert:
	${PWD}/convert_chroms.sh ${LOCAL_TFBS_BED} ${THRESHOLD}

upload:
	aws s3 cp --dryrun . ${S3_DEST_URL} --recursive --exclude "*" --include "${LOCAL_TFBS_PREFIX}.chr*.bam"
	aws s3 cp --dryrun . ${S3_DEST_URL} --recursive --exclude "*" --include "${LOCAL_TFBS_PREFIX}.chr*.bam.bai"

upload-real:
	aws s3 cp . ${S3_DEST_URL} --recursive --exclude "*" --include "${LOCAL_TFBS_PREFIX}.chr*.bam"
	aws s3 cp . ${S3_DEST_URL} --recursive --exclude "*" --include "${LOCAL_TFBS_PREFIX}.chr*.bam.bai"
