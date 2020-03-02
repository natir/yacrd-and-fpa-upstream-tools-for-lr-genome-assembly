#!/bin/bash

mkdir -p references
mkdir -p data

# reference

## E. coli CFT073
curl "https://www.ebi.ac.uk/ena/data/view/AE014075&display=fasta" | seqtk seq -A - > references/SRR8494940_ont.fasta
ln -s references/SRR8494940_ont.fasta references/SRR8494911_pb.fasta
ln -s references/SRR8494940_ont.fasta references/ref_e_coli_cft073.fasta

## H. sapiens chr1
curl ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz | seqtk seq -A - > references/h_sapiens_chr1_ref.fasta

## D. melanogaster
curl "https://www.ebi.ac.uk/ena/data/view/AE014298.5,AE014134.6,AE013599.5,AE014296.5,AE014297.3,AE014135.4,CP007106.1,KJ947872.2&display=fasta" | seqtk seq -A - > references/d_melanogaster_reads_ont_ref.fasta

## C. elegans
curl "https://www.ebi.ac.uk/ena/data/view/BX284601.5,BX284602.5,BX284603.4,BX284604.4,BX284605.5,BX284606.5&display=fasta" | seqtk seq -A - > references/c_elegans_pb_ref.fasta

# reads

## D. melanogaster
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR670/003/SRR6702603/SRR6702603_1.fastq.gz | seqtk seq -A - > data/d_melanogaster_reads_ont.fasta

## H. sapiens
curl http://s3.amazonaws.com/nanopore-human-wgs/chr1.sorted.bam | samtools bam2fq | seqtk seq -A - > data/h_sapiens_chr1_ont.fasta

## C. elegans
rm data/c_elegans_pb.fasta
for i in $(curl http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/wget.html | grep "fasta" | cut -d\" -f2 | grep "40X")
do
    curl "http://datasets.pacb.com.s3.amazonaws.com${i}" >> data/c_elegans_pb.fasta
done
    
## PRJEB6403
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR253/003/ERR2531093/ERR2531093_subreads.fastq.gz | seqtk seq -A - > data/ERR2531093_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR253/004/ERR2531094/ERR2531094_subreads.fastq.gz | seqtk seq -A - > data/ERR2531094_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR253/005/ERR2531095/ERR2531095_subreads.fastq.gz | seqtk seq -A - > data/ERR2531095_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR253/006/ERR2531096/ERR2531096_subreads.fastq.gz | seqtk seq -A - > data/ERR2531096_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR256/004/ERR2564864/ERR2564864_subreads.fastq.gz | seqtk seq -A - > data/ERR2564864_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR256/005/ERR2564865/ERR2564865_subreads.fastq.gz | seqtk seq -A - > data/ERR2564865_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR260/003/ERR2603033/ERR2603033_subreads.fastq.gz | seqtk seq -A - > data/ERR2603033_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR261/003/ERR2615933/ERR2615933_subreads.fastq.gz | seqtk seq -A - > data/ERR2615933_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR261/004/ERR2615934/ERR2615934_subreads.fastq.gz | seqtk seq -A - > data/ERR2615934_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR263/009/ERR2632359/ERR2632359_subreads.fastq.gz | seqtk seq -A - > data/ERR2632359_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR265/006/ERR2651536/ERR2651536_subreads.fastq.gz | seqtk seq -A - > data/ERR2651536_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR267/004/ERR2672424/ERR2672424_subreads.fastq.gz | seqtk seq -A - > data/ERR2672424_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR267/005/ERR2672425/ERR2672425_subreads.fastq.gz | seqtk seq -A - > data/ERR2672425_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR268/000/ERR2681660/ERR2681660_subreads.fastq.gz | seqtk seq -A - > data/ERR2681660_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR269/007/ERR2695057/ERR2695057_subreads.fastq.gz | seqtk seq -A - > data/ERR2695057_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR269/008/ERR2695058/ERR2695058_subreads.fastq.gz | seqtk seq -A - > data/ERR2695058_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/006/ERR3253076/ERR3253076_subreads.fastq.gz | seqtk seq -A - > data/ERR3253076_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/008/ERR3253078/ERR3253078_subreads.fastq.gz | seqtk seq -A - > data/ERR3253078_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/009/ERR3253079/ERR3253079_subreads.fastq.gz | seqtk seq -A - > data/ERR3253079_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/000/ERR3253080/ERR3253080_subreads.fastq.gz | seqtk seq -A - > data/ERR3253080_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/001/ERR3253081/ERR3253081_subreads.fastq.gz | seqtk seq -A - > data/ERR3253081_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/002/ERR3253082/ERR3253082_subreads.fastq.gz | seqtk seq -A - > data/ERR3253082_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/003/ERR3253083/ERR3253083_subreads.fastq.gz | seqtk seq -A - > data/ERR3253083_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/004/ERR3253084/ERR3253084_subreads.fastq.gz | seqtk seq -A - > data/ERR3253084_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/005/ERR3253085/ERR3253085_subreads.fastq.gz | seqtk seq -A - > data/ERR3253085_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/007/ERR3253087/ERR3253087_subreads.fastq.gz | seqtk seq -A - > data/ERR3253087_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/008/ERR3253088/ERR3253088_subreads.fastq.gz | seqtk seq -A - > data/ERR3253088_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/009/ERR3253089/ERR3253089_subreads.fastq.gz | seqtk seq -A - > data/ERR3253089_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/001/ERR3253091/ERR3253091_subreads.fastq.gz | seqtk seq -A - > data/ERR3253091_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/003/ERR3253093/ERR3253093_subreads.fastq.gz | seqtk seq -A - > data/ERR3253093_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/004/ERR3253094/ERR3253094_subreads.fastq.gz | seqtk seq -A - > data/ERR3253094_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/005/ERR3253095/ERR3253095_subreads.fastq.gz | seqtk seq -A - > data/ERR3253095_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/006/ERR3253096/ERR3253096_subreads.fastq.gz | seqtk seq -A - > data/ERR3253096_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/007/ERR3253097/ERR3253097_subreads.fastq.gz | seqtk seq -A - > data/ERR3253097_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/008/ERR3253098/ERR3253098_subreads.fastq.gz | seqtk seq -A - > data/ERR3253098_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/009/ERR3253099/ERR3253099_subreads.fastq.gz | seqtk seq -A - > data/ERR3253099_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/000/ERR3253100/ERR3253100_subreads.fastq.gz | seqtk seq -A - > data/ERR3253100_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/002/ERR3253102/ERR3253102_subreads.fastq.gz | seqtk seq -A - > data/ERR3253102_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/003/ERR3253103/ERR3253103_subreads.fastq.gz | seqtk seq -A - > data/ERR3253103_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/004/ERR3253104/ERR3253104_subreads.fastq.gz | seqtk seq -A - > data/ERR3253104_pb.fasta

# Some dataset need to be merged
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR350/004/ERR3500074/ERR3500074_subreads.fastq.gz | seqtk seq -A - > data/ERR3500074_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR350/005/ERR3500075/ERR3500075_subreads.fastq.gz | seqtk seq -A - >> data/ERR3500074_pb.fasta

curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR265/005/ERR2651535/ERR2651535_subreads.fastq.gz | seqtk seq -A - > data/ERR2651535_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR265/007/ERR2651537/ERR2651537_subreads.fastq.gz | seqtk seq -A - >> data/ERR2651535_pb.fasta

curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/007/ERR3253077/ERR3253077_subreads.fastq.gz | seqtk seq -A - > data/ERR3253077_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/002/ERR3253092/ERR3253092_subreads.fastq.gz | seqtk seq -A - >> data/ERR3253077_pb.fasta

curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/000/ERR3253090/ERR3253090_subreads.fastq.gz | seqtk seq -A - > data/ERR3253090_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/005/ERR3253105/ERR3253105_subreads.fastq.gz | seqtk seq -A - >> data/ERR3253090_pb.fasta

curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/006/ERR3253086/ERR3253086_subreads.fastq.gz | seqtk seq -A - > data/ERR3253086_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/ERR325/001/ERR3253101/ERR3253101_subreads.fastq.gz | seqtk seq -A - >> data/ERR3253086_pb.fasta

## PRJNA422511

### ONT
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/005/SRR8494915/SRR8494915_1.fastq.gz | seqtk seq -A - > data/SRR8494915_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/006/SRR8494916/SRR8494916_1.fastq.gz | seqtk seq -A - > data/SRR8494916_ont.fasta 
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/007/SRR8494917/SRR8494917_1.fastq.gz | seqtk seq -A - > data/SRR8494917_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/008/SRR8494918/SRR8494918_1.fastq.gz | seqtk seq -A - > data/SRR8494918_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/009/SRR8494919/SRR8494919_1.fastq.gz | seqtk seq -A - > data/SRR8494919_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/000/SRR8494920/SRR8494920_1.fastq.gz | seqtk seq -A - > data/SRR8494920_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/001/SRR8494921/SRR8494921_1.fastq.gz | seqtk seq -A - > data/SRR8494921_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/002/SRR8494922/SRR8494922_1.fastq.gz | seqtk seq -A - > data/SRR8494922_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/003/SRR8494923/SRR8494923_1.fastq.gz | seqtk seq -A - > data/SRR8494923_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/004/SRR8494924/SRR8494924_1.fastq.gz | seqtk seq -A - > data/SRR8494924_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/005/SRR8494935/SRR8494935_1.fastq.gz | seqtk seq -A - > data/SRR8494935_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/006/SRR8494936/SRR8494936_1.fastq.gz | seqtk seq -A - > data/SRR8494936_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/007/SRR8494937/SRR8494937_1.fastq.gz | seqtk seq -A - > data/SRR8494937_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/008/SRR8494938/SRR8494938_1.fastq.gz | seqtk seq -A - > data/SRR8494938_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/009/SRR8494939/SRR8494939_1.fastq.gz | seqtk seq -A - > data/SRR8494939_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/000/SRR8494940/SRR8494940_1.fastq.gz | seqtk seq -A - > data/SRR8494940_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/001/SRR8494941/SRR8494941_1.fastq.gz | seqtk seq -A - > data/SRR8494941_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/002/SRR8494942/SRR8494942_1.fastq.gz | seqtk seq -A - > data/SRR8494942_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/003/SRR8494943/SRR8494943_1.fastq.gz | seqtk seq -A - > data/SRR8494943_ont.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/004/SRR8494944/SRR8494944_1.fastq.gz | seqtk seq -A - > data/SRR8494944_ont.fasta
ln -s SRR8494940_ont.fasta data/real_reads_ont.fasta

### Pacbio

curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/005/SRR8494905/SRR8494905_subreads.fastq.gz | seqtk seq -A - > data/SRR8494905_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/006/SRR8494906/SRR8494906_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494906_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/007/SRR8494907/SRR8494907_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494907_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/008/SRR8494908/SRR8494908_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494908_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/009/SRR8494909/SRR8494909_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494909_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/000/SRR8494910/SRR8494910_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494910_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/001/SRR8494911/SRR8494911_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494911_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/002/SRR8494912/SRR8494912_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494912_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/003/SRR8494913/SRR8494913_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494913_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/004/SRR8494914/SRR8494914_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494914_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/005/SRR8494925/SRR8494925_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494925_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/006/SRR8494926/SRR8494926_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494926_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/007/SRR8494927/SRR8494927_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494927_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/008/SRR8494928/SRR8494928_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494928_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/009/SRR8494929/SRR8494929_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494929_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/000/SRR8494930/SRR8494930_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494930_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/001/SRR8494931/SRR8494931_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494931_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/002/SRR8494932/SRR8494932_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494932_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/003/SRR8494933/SRR8494933_subreads.fastq.gz | seqtk seq -A - >	data/SRR8494933_pb.fasta
curl ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/004/SRR8494934/SRR8494934_subreads.fastq.gz | seqtk seq -A - > data/SRR8494934_pb.fasta
ln -s SRR8494911_pb.fasta data/real_reads_pb.fasta


