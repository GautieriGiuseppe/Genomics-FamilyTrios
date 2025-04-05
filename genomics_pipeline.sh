# list of assigned cases
cases=("case747" "case685" "case714" "case629" "case723" "case584" "case688" "case696" "case590" "case730")

# path to useful file
REFERENCE="/home/BCG2025_genomics_exam/universe.fasta"
BOWTIE_INDEX="/home/BCG2025_cacciaR/all_var/uni"
EXONS_BED="/home/BCG2025_genomics_exam/exons16Padded_sorted.bed"
VEP_CACHE="/home/BCG2025_genomics_exam/vep_cache"

# index creation if not already present
if [ ! -f "${BOWTIE_INDEX}.1.bt2" ]; then
    echo "Building Bowtie2 index..."
    bowtie2-build "$REFERENCE" "$BOWTIE_INDEX"
fi

# loop for every case
for case in "${cases[@]}"; do
    echo "Processing $case..."

    mkdir -p QC_reports_"${case}"
    mkdir -p VCF_outputs_"${case}"

    mother="/home/BCG2025_genomics_exam/${case}_mother.fq.gz"
    father="/home/BCG2025_genomics_exam/${case}_father.fq.gz"
    child="/home/BCG2025_genomics_exam/${case}_child.fq.gz"

    # FastQC control
    fastqc "$mother" "$father" "$child" -o QC_reports_"${case}"/

    # Bowtie2 aligment
    bowtie2 -x "$BOWTIE_INDEX" -U "$mother" -S "${case}_mother.sam" --rg-id "mother" --rg "SM:mother"
    bowtie2 -x "$BOWTIE_INDEX" -U "$father" -S "${case}_father.sam" --rg-id "father" --rg "SM:father"
    bowtie2 -x "$BOWTIE_INDEX" -U "$child" -S "${case}_child.sam" --rg-id "child" --rg "SM:child"

    # SAM-BAM conversion
    samtools view -bS "${case}_mother.sam" | samtools sort -o "${case}_mother_sorted.bam"
    samtools view -bS "${case}_father.sam" | samtools sort -o "${case}_father_sorted.bam"
    samtools view -bS "${case}_child.sam" | samtools sort -o "${case}_child_sorted.bam"

    # BAM indexing
    samtools index "${case}_mother_sorted.bam"
    samtools index "${case}_father_sorted.bam"
    samtools index "${case}_child_sorted.bam"

    # Qualimap
    qualimap bamqc -bam "${case}_child_sorted.bam" --gff "$EXONS_BED" -outdir "QC_reports_${case}/child"
    qualimap bamqc -bam "${case}_father_sorted.bam" --gff "$EXONS_BED" -outdir "QC_reports_${case}/father"
    qualimap bamqc -bam "${case}_mother_sorted.bam" --gff "$EXONS_BED" -outdir "QC_reports_${case}/mother"
    
    # VCF creation
    freebayes -f "$REFERENCE" "${case}_mother_sorted.bam" "${case}_father_sorted.bam" "${case}_child_sorted.bam" > "VCF_outputs_${case}/${case}_trio.vcf"

    # Generation of the VCF file keeping just the variants included in the target region
    bedtools intersect -a "VCF_outputs_${case}/${case}_trio.vcf" -b "$EXONS_BED" -u > "VCF_outputs_${case}/${case}_filtered.vcf"

    # VEP annotation
    ~/.vep/vep -i "VCF_outputs_${case}/${case}_filtered.vcf" --cache --dir "$VEP_CACHE" --o "VCF_outputs_${case}/${case}_annotated.vcf" --vcf

    # MultiQC per report
    multiqc QC_reports_"${case}"/

    echo "Completed processing for $case"
done

echo "Pipeline finished!"