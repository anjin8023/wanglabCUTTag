usage(){
        echo "
        Example: $0 treatment control
        "
}

treat=$1
control=$2


bowtie2_index=/home/liukuai/Documents/reference/bowtie2/mm10
picard=/home/liukuai/Documents/sofware/picard/build/libs/picard.jar
tssbed=/home/liukuai/Documents/reference/ucsc_mm10_refseq/ucsc.mm10.refseq.tss.bed
# tssbed=/home/anjin/Documents/project/reference/ucsc_mm10_refseq/ucsc.mm10.cpg.island.txt
workdir=`pwd`
outdir=$workdir/macs2
blacklist=/home/liukuai/Documents/reference/blacklist_region/mm10-blacklist.withchr.v2.bed


mkdir logdir
mkdir aligndir
mkdir sorted
mkdir dup.marked dedup
mkdir bigwig
mkdir bigwigNorm
mkdir $outdir
mkdir plot 
mkdir diffPeaks
mkdir blacklist_filtered

process_align(){
        base=$1
        date | sed 's/$/   XXX           map to genome   started.../' | sed "s/XXX/$name/" >> logdir/"$base".log;
        (bowtie2  --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 2 -x $bowtie2_index -1 fastq/"$base"_1.clean.fq.gz -2 fastq/"$base"_2.clean.fq.gz) 2> logdir/"$base".bowtie2 | samtools view -bS - > aligndir/"$base"_aligned_reads.bam
        ## Filtering unmapped fragments
        samtools view -bh -f 3 -F 4 -F 8 aligndir/"$base"_aligned_reads.bam > sorted/"$base".step1.bam
        ##sorting BAM
        java -jar $picard SortSam -I sorted/"$base".step1.bam -O sorted/"$base".bam -SORT_ORDER coordinate -VALIDATION_STRINGENCY SILENT
        rm -rf sorted/"$base".step1.bam 
        ##Marking duplicates
        java -jar $picard MarkDuplicates \
                -INPUT sorted/"$base".bam -OUTPUT dup.marked/"$base".bam -VALIDATION_STRINGENCY SILENT -METRICS_FILE dup.marked/metrics."$base".txt
        ##Removing duplicates
        samtools view -bh -F 1024 dup.marked/"$base".bam  > dedup/"$base".bam
        ##Remove blacklist region
        bedtools intersect  -v -abam dedup/"$base".bam -b $blacklist > blacklist_filtered/"$base".bam
        ##Creating bam index files
        # samtools index sorted/"$base".bam
        # samtools index dup.marked/"$base".bam
        samtools index dedup/"$base".bam
        samtools index blacklist_filtered/"$base".bam
        #### Convert into bigwig file format
        bamCoverage -b blacklist_filtered/"$base".bam -o bigwig/"$base".bw
        bamCoverage -b blacklist_filtered/"$base".bam --normalizeUsing RPKM -o bigwigNorm/"$base".bw 
        # rm -rf sorted/"$base".bam  dup.marked/"$base".bam  aligndir/"$base".bam  
}




macs2_call_peaks(){
    base=$1
    date | sed 's/$/   XXX         Peak calling using MACS2.../' | sed "s/XXX/$base/" >> logdir/"$base".log;
    macs2 callpeak -t blacklist_filtered/"$base".bam  -g mm -f BAMPE -n $base  --outdir $outdir -q 0.05 -B --SPMR --keep-dup all 2> logdir/"$base".macs2
    ####broad peak calls
    # macs2 callpeak -t blacklist_filtered/"$base".bam  -g mm -f BAMPE -n $base  --outdir $outdirbroad --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all 2> logdir/"$base".broad.all.frag.macs2
    }



visualization(){
    computeMatrix scale-regions -S bigwigNorm/"$treat".bw  bigwigNorm/"$control".bw -R $tssbed --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o plot/"$treat".mat.gz -p 2
    plotProfile -m plot/"$treat".mat.gz -out plot/"$treat".png  --perGroup  --plotTitle ""  --samplesLabel "$treat" "$control"  -T "$treat"  -z ""  --startLabel ""  --endLabel ""  --colors red darkblue

}



differential_peaks(){
    mkdir -p  diffPeaks/$treat
    manorm --pe --rf bam --pf macs2 --p1 $outdir/"$control"_peaks.xls --p2 $outdir/"$treat"_peaks.xls --r1 blacklist_filtered/"$control".bam --r2 blacklist_filtered/"$treat".bam --n1 "$control" --n2 "$treat" -o diffPeaks/$treat
}




process_align $treat &
process_align $control &

wait

macs2_call_peaks $treat &
macs2_call_peaks $control &

wait

visualization &
# differential_peaks &


# bash /home/liukuai/Documents/pipeline/cut_tag/cut_and_tag_pipeline_without_trim.sh  H3K4-KO_FKDL202609935-1a H3K4-WT_FKDL202609934-1a &

