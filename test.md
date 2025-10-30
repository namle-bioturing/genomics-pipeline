```bash
bam=${1}
sample_id=${2}
region=autosomes
output_file_path="${sample_id}.out"
output_summary_file_path="${sample_id}.summary"
output_rdata="${sample_id}.rdata"

Rscript /data/GL/nextflow/nextflow/bin/ExomeDepth/export_bam_counts.R ${bam} ${sample_id} ${region} ${output_file_path} ${output_summary_file_path} ${output_rdata}
```

docker run -tiv .:/workspace -v /data/GL:/data/GL namxle/batchcnv-tidyverse:4.4.0 bash

nohup /bin/bash test.sh 3997.deduped.bam 3997 &
nohup /bin/bash test.sh 4549.deduped.bam 4549 > 4549.log &
nohup /bin/bash test.sh 4653.deduped.bam 4653 > 4653.log &

Rscript /data/GL/nextflow/nextflow/bin/ExomeDepth/call_exome_cnv_pca.R 4653.rdata autosomes 4653.csv 4653.call.rdata baseline