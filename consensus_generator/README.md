# Gene Consensus Generator (GCG) Tool

A Python tool for comparing gene predictions from MetaGeneMark and Prodigal tools, classifying genes into different categories based on their overlap patterns.

## Overview

The GCG tool analyzes gene prediction results from two popular gene finding tools:
- **MetaGeneMark**: A gene prediction tool for prokaryotic and eukaryotic genomes
- **Prodigal**: A gene prediction tool specifically designed for prokaryotic genomes

The tool classifies genes into four main categories:
1. **Consensus**: Genes predicted by both tools
2. **Insertio**: Genes completely contained within genes from the other tool
3. **Concomitantiam**: Genes that partially overlap with genes from the other tool
4. **Exclusivus**: Genes unique to one tool with no overlap

## Features

- **Comprehensive Analysis**: Compares gene predictions from MetaGeneMark and Prodigal
- **Multiple Output Formats**: Generates GFF files for each gene category
- **Flexible Output**: Optional generation of different gene categories
- **Error Handling**: Robust error handling and validation
- **Progress Reporting**: Clear status messages and statistics

## Requirements

### Python Dependencies
```bash
pip install pandas
```

### Input Files
- **MetaGeneMark GFF file**: Gene predictions from MetaGeneMark tool
- **Prodigal GFF file**: Gene predictions from Prodigal tool

Both input files must be in GFF3 format and contain CDS features.

## Installation

1. Clone or download the script
2. Ensure Python 3.6+ is installed
3. Install required dependencies:
   ```bash
   pip install pandas
   ```

## Usage

### Basic Usage

Generate consensus and unique gene files:
```bash
python GCG_final.py -m metagenemark.gff -p prodigal.gff -o output_directory
```

### Advanced Usage

Generate all gene categories:
```bash
python GCG_final.py -m metagenemark.gff -p prodigal.gff -o output_directory -i -c -e
```

### Command Line Options

| Option | Long Option | Description | Required |
|--------|-------------|-------------|----------|
| `-m` | `--metagenemark` | MetaGeneMark GFF file | Yes |
| `-p` | `--prodigal` | Prodigal GFF file | Yes |
| `-o` | `--output` | Output directory path | Yes |
| `-i` | `--insertio` | Generate insertio gene files | No |
| `-c` | `--concomitantiam` | Generate concomitantiam gene files | No |
| `-e` | `--exclusivus` | Generate exclusivus gene files | No |

## Output Files

The tool generates the following output files in the specified directory:

### Always Generated
- `{output_name}_consensus_metagenemark.gff`: Consensus genes from MetaGeneMark perspective
- `{output_name}_consensus_prodigal.gff`: Consensus genes from Prodigal perspective
- `{output_name}_metagenemark_notmerged.gff`: Unique MetaGeneMark genes
- `{output_name}_prodigal_notmerged.gff`: Unique Prodigal genes

### Optional Files (with `-i`, `-c`, `-e` flags)
- `{output_name}_insertio.gff`: Inserted genes (completely contained)
- `{output_name}_concomitantiam.gff`: Partially overlapping genes
- `{output_name}_exclusivus.gff`: Exclusive genes (no overlap)

## Gene Classification Logic

### Consensus Genes
Genes that are predicted by both tools with identical start and end positions.

### Insertio Genes
Genes that are completely contained within genes from the other tool:
- MetaGeneMark gene is completely within a Prodigal gene
- Prodigal gene is completely within a MetaGeneMark gene

### Concomitantiam Genes
Genes that partially overlap with genes from the other tool:
- Partial overlap at either start or end
- One gene extends beyond the other at both ends

### Exclusivus Genes
Genes that are unique to one tool and have no overlap with any genes from the other tool.

## Example Workflow

1. **Prepare Input Files**:
   ```bash
   # Run MetaGeneMark
   metagenemark your_genome.fasta
   
   # Run Prodigal
   prodigal -i your_genome.fasta -o prodigal.gff -f gff
   ```

2. **Run GCG Analysis**:
   ```bash
   python GCG_final.py -m metagenemark.gff -p prodigal.gff -o results -i -c -e
   ```

3. **Review Results**:
   - Check consensus genes for high-confidence predictions
   - Analyze insertio genes for potential annotation differences
   - Examine concomitantiam genes for partial overlaps
   - Review exclusivus genes for tool-specific predictions

## Output Statistics

The tool provides statistics for each gene category:
```
Processing files:
  MetaGeneMark: metagenemark.gff
  Prodigal: prodigal.gff
  Output: results

Consensus MetaGeneMark genes: 1500
Consensus Prodigal genes: 1500
Unique MetaGeneMark genes: 200
Unique Prodigal genes: 300
Insertio MetaGeneMark genes: 50
Insertio Prodigal genes: 75
Concomitantiam MetaGeneMark genes: 100
Concomitantiam Prodigal genes: 125
Exclusivus MetaGeneMark genes: 150
Exclusivus Prodigal genes: 225

Analysis complete!
```

## File Format

All output files are in GFF3 format with the following columns:
1. **SeqID**: Sequence identifier
2. **Source**: Tool name (MetaGeneMark or Prodigal)
3. **Feature**: Feature type (CDS)
4. **Start**: Start position
5. **End**: End position
6. **Score**: Prediction score
7. **Strand**: Strand orientation (+ or -)
8. **Frame**: Reading frame
9. **Attributes**: Additional attributes

## Error Handling

The tool includes comprehensive error handling:
- **File Validation**: Checks if input files exist
- **Format Validation**: Validates GFF format
- **Data Type Conversion**: Handles coordinate conversions safely
- **Memory Management**: Efficient DataFrame operations

## Troubleshooting

### Common Issues

1. **"File not found" error**:
   - Ensure input file paths are correct
   - Check file permissions

2. **"Import pandas could not be resolved"**:
   - Install pandas: `pip install pandas`

3. **Empty output files**:
   - Verify input files contain CDS features
   - Check GFF format compliance

4. **Memory errors with large files**:
   - Consider processing smaller genome segments
   - Ensure sufficient system memory

## Contributing

To contribute to this project:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This tool is provided as-is for research and educational purposes.

## Citation

If you use this tool in your research, please cite:
- MetaGeneMark: Zhu et al. (2010) MetaGeneMark: A tool for gene prediction in metagenomic sequences
- Prodigal: Hyatt et al. (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification

## Contact

For questions or issues, please send email to emersonwdanzer@gmail.com.

---

**Note**: This tool is designed for comparative analysis of gene prediction tools and should be used as part of a comprehensive genome annotation pipeline.
