# all cases
cases=("case747" "case685" "case714" "case629" "case723" "case584" "case688" "case696" "case590" "case730")

# recessive cases
recessive_cases="case688 case747 case685"

# reusable path
REFERENCE="/home/BCG2025_genomics_exam/universe.fasta"
BOWTIE_INDEX="/home/BCG2025_cacciaR/all_var/uni"
EXONS_BED="/home/BCG2025_genomics_exam/exons16Padded_sorted.bed"
VEP_CACHE="/home/BCG2025_genomics_exam/vep_cache"

# indexing building if not present
if [ ! -f "${BOWTIE_INDEX}.1.bt2" ]; then
    echo "Building Bowtie2 index..."
    bowtie2-build "$REFERENCE" "$BOWTIE_INDEX"
fi

# for loop for each case
for case in "${cases[@]}"; do
    echo "Processing $case..."

    mkdir -p QC_reports_"$case"
    mkdir -p VCF_outputs_"$case"
    mkdir -p coverage_"$case"

    mother="/home/BCG2025_genomics_exam/${case}_mother.fq.gz"
    father="/home/BCG2025_genomics_exam/${case}_father.fq.gz"
    child="/home/BCG2025_genomics_exam/${case}_child.fq.gz"

    # FastQC report
    fastqc "$mother" "$father" "$child" -o QC_reports_"$case"/

    # Bowtie2 alignment
    bowtie2 -x "$BOWTIE_INDEX" -U "$mother" -S "${case}_mother.sam" --rg-id "mother" --rg "SM:mother"
    bowtie2 -x "$BOWTIE_INDEX" -U "$father" -S "${case}_father.sam" --rg-id "father" --rg "SM:father"
    bowtie2 -x "$BOWTIE_INDEX" -U "$child"  -S "${case}_child.sam"  --rg-id "child"  --rg "SM:child"

    # SAM → sorted BAM
    samtools view -bS "${case}_mother.sam" | samtools sort -o "${case}_mother_sorted.bam"
    samtools view -bS "${case}_father.sam" | samtools sort -o "${case}_father_sorted.bam"
    samtools view -bS "${case}_child.sam"  | samtools sort -o "${case}_child_sorted.bam"

    # BAM indexing
    samtools index "${case}_mother_sorted.bam"
    samtools index "${case}_father_sorted.bam"
    samtools index "${case}_child_sorted.bam"

    # Qualimap
    qualimap bamqc -bam "${case}_child_sorted.bam" --gff "$EXONS_BED" -outdir "QC_reports_${case}/child"
    qualimap bamqc -bam "${case}_father_sorted.bam" --gff "$EXONS_BED" -outdir "QC_reports_${case}/father"
    qualimap bamqc -bam "${case}_mother_sorted.bam" --gff "$EXONS_BED" -outdir "QC_reports_${case}/mother"

    # Variant calling
    freebayes -f "$REFERENCE" "${case}_mother_sorted.bam" "${case}_father_sorted.bam" "${case}_child_sorted.bam" > "VCF_outputs_${case}/${case}_trio.vcf"

    # Filtering on target region
    bedtools intersect -a "VCF_outputs_${case}/${case}_trio.vcf" -b "$EXONS_BED" -u > "VCF_outputs_${case}/${case}_filtered.vcf"

    # Annotation with VEP
    ~/.vep/vep -i "VCF_outputs_${case}/${case}_filtered.vcf" --cache --dir "$VEP_CACHE" --o "VCF_outputs_${case}/${case}_annotated.vcf" --vcf

    # AR/AD classification
    if [[ " $recessive_cases " =~ " $case " ]]; then
        echo "-> Recessive case detected for $case"
        grep "#" "VCF_outputs_${case}/${case}_annotated.vcf" > "VCF_outputs_${case}/${case}_recessive.vcf"
        grep "1/1.*0/1.*0/1" "VCF_outputs_${case}/${case}_annotated.vcf" >> "VCF_outputs_${case}/${case}_recessive.vcf"
    else
        echo "-> Dominant case detected for $case"
        grep "#" "VCF_outputs_${case}/${case}_trio.vcf" > "VCF_outputs_${case}/${case}_dominant.vcf"
        grep "0/1.*0/0.*0/0" "VCF_outputs_${case}/${case}_filtered.vcf" >> "VCF_outputs_${case}/${case}_dominant.vcf"
    fi

    # MultiQC
    multiqc QC_reports_"$case"/

    # Coverage UCSC
    bedtools genomecov -ibam "${case}_child_sorted.bam"  -bg -trackline -trackopts 'name="child"'  -max 100 > "coverage_${case}/child${case}Cov.bg"
    bedtools genomecov -ibam "${case}_mother_sorted.bam" -bg -trackline -trackopts 'name="mother"' -max 100 > "coverage_${case}/mother${case}Cov.bg"
    bedtools genomecov -ibam "${case}_father_sorted.bam" -bg -trackline -trackopts 'name="father"' -max 100 > "coverage_${case}/father${case}Cov.bg"

    echo "Completed $case"
done

echo "Pipeline COMPLETED for all cases!"
