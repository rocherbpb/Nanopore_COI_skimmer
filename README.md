### Config
THREADS=$NSLOTS
MITO_REF="data/Reference_taxa.fasta"
GENETIC_CODE=5
POLISH=true
MEDAKA_MODEL="r941_min_hac_g507"
OUTDIR="results"
GENBANK_REF="data/Reference.gb"

### Loop over input FASTQ/set output directories
mkdir -p "$OUTDIR"

for fq in data/*.{fastq,fastq.gz}; do
  [ -e "$fq" ] || continue

  sample=$(basename "$fq")
  sample=${sample%%.*}

  echo "Processing $sample ..."

  sample_out="$OUTDIR/$sample"
  mkdir -p "$sample_out"

  ### 0. Filter reads with NanoFilt (quality ≥10, length ≥200)
  echo "Filtering reads with NanoFilt..."
  if [[ "$fq" == *.gz ]]; then
    zcat "$fq" | NanoFilt -q 10 -l 200 > "$sample_out/${sample}_filtered.fastq"
  else
    cat "$fq" | NanoFilt -q 10 -l 200 > "$sample_out/${sample}_filtered.fastq"
  fi

  filtered_reads="$sample_out/${sample}_filtered.fastq"

  ### 1. Map filtered reads to mito reference and extract mapped reads only (flag -F 4 filters unmapped)
  minimap2 -t $THREADS --secondary=no -ax map-ont "$MITO_REF" "$filtered_reads" | \
    samtools view -b -@ $THREADS -F 4 -o "$sample_out/aln.bam"

  samtools fastq "$sample_out/aln.bam" > "$sample_out/mapped_reads.fastq"

  ### 2. Assemble mitogenome from mapped reads (rnabloom2)
  rnabloom -t $THREADS -long "$sample_out/mapped_reads.fastq" -o "$sample_out/rnabloom_mito_out"

  ### 3. Map mapped reads back to rnabloom assembly for polishing input
  minimap2 -t $THREADS -ax map-ont "$sample_out/rnabloom_mito_out/rnabloom.transcripts.fa" "$sample_out/mapped_reads.fastq" > "$sample_out/rnabloom_aln.sam"

  ### 4. Polish assembly with racon using mapped_reads.fastq and rnabloom_aln.sam
  racon -t $THREADS "$sample_out/mapped_reads.fastq" "$sample_out/rnabloom_aln.sam" "$sample_out/rnabloom_mito_out/rnabloom.transcripts.fa" > "$sample_out/racon_round1.fa"

  ### 5. Polish with medaka if flag true
  if [ "$POLISH" = true ]; then
    medaka_consensus -i "$sample_out/mapped_reads.fastq" -d "$sample_out/racon_round1.fa" -o "$sample_out/medaka_out" -t $THREADS -m "$MEDAKA_MODEL"
    final_asm="medaka_out/consensus.fasta"
  else
    final_asm="rnabloom_mito_out/rnabloom.transcripts.fa"
  fi

  ### 6. Run mitofinder from inside sample folder, including GenBank reference
  cd "$sample_out"
  mitofinder -j "$sample" -a "$final_asm" -p $THREADS  -o $GENETIC_CODE -r "../../$GENBANK_REF"
  cd - > /dev/null

  ### 7. Extract COX1 sequences from mitofinder results using seqkit grep
  seqkit grep -r -i -p "cox1|coxi|cytochrome c oxidase subunit 1" \
    $(find "$sample_out/$sample" -name "*final_genes_NT.fasta") > "$sample_out/cox1_extracted.fasta"

  ### 8. Calculate mean coverage of COX1 sequences using coverm
  coverm contig --reference "$sample_out/cox1_extracted.fasta" --single "$filtered_reads" --mapper minimap2-ont \
    --methods mean --output-file "$sample_out/cox1_coverage.tsv" --threads $THREADS --min-read-percent-identity 80

  echo "$sample done."
done
