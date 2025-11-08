## Install VEP


### Install cache and plugins

```bash
docker pull ensemblorg/ensembl-vep
```

```bash
    docker run -u $(id -u):$(id -g) --rm -itd --name vep_installation \
    -v /mnt/nasdev2/namle/.vep:/opt/vep/.vep \
    ensemblorg/ensembl-vep \
    INSTALL.pl -a cpf --species homo_sapiens_merged --assembly GRCh38 -c /opt/vep/.vep --PLUGINS all -r /opt/vep/.vep/plugins
```

```bash
docker run -u $(id -u):$(id -g) -itd --name vep_installation \
    -v /mnt/nasdev2/namle/.vep:/opt/vep/.vep \
    ensemblorg/ensembl-vep \
    INSTALL.pl -a p --species homo_sapiens_merged --assembly GRCh38 -c /opt/vep/.vep --PLUGINS all -r /opt/vep/.vep/plugins
```

### Install plugins data

#### 1. AlphaMissense

```bash
# Create storage
mkdir -p /mnt/nasdev2/namle/.vep/plugins_data/AlphaMissense
cd /mnt/nasdev2/namle/.vep/plugins_data/AlphaMissense

# For hg38
wget https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz .
tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.tsv.gz
bgzip -dc -@ 10 AlphaMissense_hg38.tsv.gz \
| awk 'BEGIN{FS=OFS="\t"} { sub(/^chr/, "", $1) } 1' \
| bgzip -@ 10 -c > AlphaMissense_hg38.modified.tsv.gz
tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.modified.tsv.gz

# For hg19
wget https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz .
tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg19.tsv.gz
bgzip -dc -@ 10 AlphaMissense_hg19.tsv.gz \
| awk 'BEGIN{FS=OFS="\t"} { sub(/^chr/, "", $1) } 1' \
| bgzip -@ 10 -c > AlphaMissense_hg19.modified.tsv.gz
tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg19.modified.tsv.gz
```

#### 2. CADD

```bash
mkdir -p /mnt/nasdev2/namle/.vep/plugins_data/CADD
cd  /mnt/nasdev2/namle/.vep/plugins_data/CADD

# For hg38
# SNVs
nohup wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz -O whole_genome_SNVs_v1.7_hg38.tsv.gz &
wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi -O whole_genome_SNVs_v1.7_hg38.tsv.gz.tbi
# InDels
nohup wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz -O gnomad.genomes.r4.0.indel.v1.7.hg38.tsv.gz &
wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz.tbi -O gnomad.genomes.r4.0.indel.v1.7.hg38.tsv.gz.tbi

# For hg19
# SNVs
nohup wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz -O whole_genome_SNVs_v1.7_hg19.tsv.gz &
wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz.tbi -O whole_genome_SNVs_v1.7_hg19.tsv.gz.tbi
# InDels
wget https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/gnomad.genomes-exomes.r4.0.indel.tsv.gz -O gnomad.genomes-exomes.r4.0.indel.v1.7.hg19.tsv.gz
wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh37/gnomad.genomes-exomes.r4.0.indel.tsv.gz.tbi -O gnomad.genomes-exomes.r4.0.indel.v1.7.hg19.tsv.gz.tbi
```

#### 3. REVEL

```bash
# Create storage
mkdir -p /mnt/nasdev2/namle/.vep/plugins_data/REVEL
cd /mnt/nasdev2/namle/.vep/plugins_data/REVEL

# Download
wget https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip

# Extract data - follow VEP instructions
unzip revel-v1.3_all_chromosomes.zip
cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
bgzip new_tabbed_revel.tsv

tabix -f -s 1 -b 2 -e 2 new_tabbed_revel.tsv.gz

zcat new_tabbed_revel.tsv.gz | head -n1 > h
zgrep -h -v ^#chr new_tabbed_revel.tsv.gz | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat h - | bgzip -c > new_tabbed_revel_grch38.tsv.gz
tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.tsv.gz

# For hg19
mv new_tabbed_revel.tsv.gz revel_hg19.tsv.gz
mv new_tabbed_revel.tsv.gz.tbi revel_hg19.tsv.gz.tbi

# For hg38
mv new_tabbed_revel_grch38.tsv.gz revel_hg38.tsv.gz
mv new_tabbed_revel_grch38.tsv.gz.tbi revel_hg38.tsv.gz.tbi

# Clean up
rm h revel_with_transcript_ids tabbed_revel.tsv
```

#### 4. SpliceAI
```bash
# Create storage 
mkdir -p /mnt/nasdev2/namle/.vep/plugins_data/SpliceAI
cd /mnt/nasdev2/namle/.vep/plugins_data/SpliceAI

# For hg38
# For SNV
wget "https://basespace-data-east.s3.us-east-1.amazonaws.com/170dc484120a49f0b897aae301891840/spliceai_scores.raw.snv.hg38.vcf.gz?AWSAccessKeyId=AKIARPYQJSWQZKU3DBUF&Expires=1762920256&response-content-disposition=attachment%3Bfilename%3Dspliceai_scores.raw.snv.hg38.vcf.gz&response-content-type=application%2Fx-gzip&Signature=%2FTJn1tJDHUSnRKwxZJFzdXWdeek%3D" -O spliceai_scores.raw.snv.hg38.vcf.gz 
wget "https://basespace-data-east.s3.us-east-1.amazonaws.com/170dc484120a49f0b897aae301891840/spliceai_scores.raw.snv.hg38.vcf.gz.tbi?AWSAccessKeyId=AKIARPYQJSWQZKU3DBUF&Expires=1762920288&response-content-disposition=attachment%3Bfilename%3Dspliceai_scores.raw.snv.hg38.vcf.gz.tbi&response-content-type=application%2Foctet-stream&Signature=FcHKgWUSy2VS%2F3XH8T5eoZeN%2FIU%3D" -O spliceai_scores.raw.snv.hg38.vcf.gz.tbi
# For InDel
wget "https://basespace-data-east.s3.us-east-1.amazonaws.com/170dc484120a49f0b897aae301891840/spliceai_scores.raw.indel.hg38.vcf.gz?AWSAccessKeyId=AKIARPYQJSWQZKU3DBUF&Expires=1762920043&response-content-disposition=attachment%3Bfilename%3Dspliceai_scores.raw.indel.hg38.vcf.gz&response-content-type=application%2Fx-gzip&Signature=Cjw9P5bSJZCI2Rw13y2%2Fcx09RRg%3D" -O spliceai_scores.raw.indel.hg38.vcf.gz
wget "https://basespace-data-east.s3.us-east-1.amazonaws.com/170dc484120a49f0b897aae301891840/spliceai_scores.raw.indel.hg38.vcf.gz.tbi?AWSAccessKeyId=AKIARPYQJSWQZKU3DBUF&Expires=1762920110&response-content-disposition=attachment%3Bfilename%3Dspliceai_scores.raw.indel.hg38.vcf.gz.tbi&response-content-type=application%2Foctet-stream&Signature=xd3LGOVvlnbAma%2F4LoYgg18UWX8%3D" -O spliceai_scores.raw.indel.hg38.vcf.gz.tbi

# For hg19
# For SNV
# For InDel
```