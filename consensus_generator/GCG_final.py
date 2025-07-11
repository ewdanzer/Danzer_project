#!/usr/bin/env python3
"""
Gene Consensus Generator (GCG) Tool

This script compares gene predictions from MetaGeneMark and Prodigal tools,
classifying genes into different categories: consensus, insertio, concomitantiam, and exclusivus.
This is a refactored version, any problem should be reported to emersonwdanzer@gmail.com
"""

import sys
import os
import re
import pandas as pd
from optparse import OptionParser
from typing import List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class GFFRecord:
    """Represents a GFF record with all its fields."""
    seqid: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    frame: str
    attributes: str


class GFFParser:
    """Handles parsing and processing of GFF files."""
    
    GFF_COLUMNS = ['SeqID', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attributes']
    
    @staticmethod
    def parse_gff_file(file_path: str) -> pd.DataFrame:
        """Parse a GFF file and return a DataFrame with CDS features only."""
        try:
            with open(file_path, 'r') as file:
                data = [re.split(r'\s+', line.strip()) for line in file if line.strip() and not line.startswith('#')]
            
            df = pd.DataFrame(data, columns=GFFParser.GFF_COLUMNS)
            return df[df['Feature'] == 'CDS']
        except FileNotFoundError:
            print(f"Error: File {file_path} not found.")
            sys.exit(1)
        except Exception as e:
            print(f"Error parsing file {file_path}: {e}")
            sys.exit(1)


class GeneClassifier:
    """Handles gene classification logic."""
    
    def __init__(self, metagenemark_df: pd.DataFrame, prodigal_df: pd.DataFrame):
        """Initialize with MetaGeneMark and Prodigal DataFrames."""
        self.metagenemark_df = metagenemark_df.copy()
        self.prodigal_df = prodigal_df.copy()
        self._prepare_dataframes()
    
    def _prepare_dataframes(self):
        """Prepare DataFrames for comparison by converting coordinates to integers."""
        for df in [self.metagenemark_df, self.prodigal_df]:
            df['Start'] = pd.to_numeric(df['Start'], errors='coerce')
            df['End'] = pd.to_numeric(df['End'], errors='coerce')
    
    def get_consensus_genes(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Find genes that are predicted by both tools (consensus)."""
        merged = pd.merge(
            self.metagenemark_df, 
            self.prodigal_df, 
            on=['SeqID', 'Start', 'End'], 
            suffixes=('_metagenemark', '_prodigal'), 
            how='inner'
        )
        
        consensus_metagenemark = merged[['SeqID', 'Source_metagenemark', 'Feature_metagenemark', 
                                       'Start', 'End', 'Score_metagenemark', 'Strand_metagenemark', 
                                       'Frame_metagenemark', 'Attributes_metagenemark']]
        consensus_metagenemark.columns = GFFParser.GFF_COLUMNS
        
        consensus_prodigal = merged[['SeqID', 'Source_prodigal', 'Feature_prodigal', 
                                   'Start', 'End', 'Score_prodigal', 'Strand_prodigal', 
                                   'Frame_prodigal', 'Attributes_prodigal']]
        consensus_prodigal.columns = GFFParser.GFF_COLUMNS
        
        return consensus_metagenemark, consensus_prodigal
    
    def get_unique_genes(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Find genes that are predicted by only one tool."""
        merged_outer = pd.merge(
            self.metagenemark_df, 
            self.prodigal_df, 
            on=['SeqID', 'Start', 'End'], 
            suffixes=('_metagenemark', '_prodigal'), 
            how='outer', 
            indicator=True
        )
        
        # MetaGeneMark unique genes
        metagenemark_unique = merged_outer[merged_outer['_merge'] == 'left_only'].drop(columns=['_merge'])
        metagenemark_unique = metagenemark_unique.drop(columns=metagenemark_unique.filter(like='_prodigal').columns)
        
        # Prodigal unique genes
        prodigal_unique = merged_outer[merged_outer['_merge'] == 'right_only'].drop(columns=['_merge'])
        prodigal_unique = prodigal_unique.drop(columns=prodigal_unique.filter(like='_metagenemark').columns)
        
        # Reorder columns for Prodigal to match standard format
        desired_order = ['SeqID', 'Source_prodigal', 'Feature_prodigal', 'Start', 'End', 
                       'Score_prodigal', 'Strand_prodigal', 'Frame_prodigal', 'Attributes_prodigal']
        column_positions = [prodigal_unique.columns.get_loc(col) for col in desired_order]
        prodigal_unique = prodigal_unique.iloc[:, column_positions]
        
        return metagenemark_unique, prodigal_unique
    
    def get_insertio_genes(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Find genes that are completely contained within genes from the other tool."""
        metagenemark_unique, prodigal_unique = self.get_unique_genes()
        
        # MetaGeneMark insertio (contained within Prodigal genes)
        merged_mgm = pd.merge(metagenemark_unique, prodigal_unique, on='SeqID', 
                             suffixes=('_metagenemark', '_prodigal'), how='inner')
        
        insertio_metagenemark = merged_mgm[
            (merged_mgm['Start_metagenemark'] >= merged_mgm['Start_prodigal']) & 
            (merged_mgm['End_metagenemark'] <= merged_mgm['End_prodigal'])
        ]
        
        if not insertio_metagenemark.empty:
            insertio_metagenemark = insertio_metagenemark.drop(columns=insertio_metagenemark.filter(like='_prodigal').columns)
            insertio_metagenemark.columns = GFFParser.GFF_COLUMNS
        
        # Prodigal insertio (contained within MetaGeneMark genes)
        merged_prod = pd.merge(prodigal_unique, metagenemark_unique, on='SeqID', 
                              suffixes=('_prodigal', '_metagenemark'), how='inner')
        
        insertio_prodigal = merged_prod[
            (merged_prod['Start_prodigal'] >= merged_prod['Start_metagenemark']) & 
            (merged_prod['End_prodigal'] <= merged_prod['End_metagenemark'])
        ]
        
        if not insertio_prodigal.empty:
            insertio_prodigal = insertio_prodigal.drop(columns=insertio_prodigal.filter(like='_metagenemark').columns)
            insertio_prodigal.columns = GFFParser.GFF_COLUMNS
        
        return insertio_metagenemark, insertio_prodigal
    
    def get_concomitantiam_genes(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Find genes that partially overlap with genes from the other tool."""
        metagenemark_unique, prodigal_unique = self.get_unique_genes()
        
        # MetaGeneMark concomitantiam (overlapping with Prodigal genes)
        merged_mgm = pd.merge(metagenemark_unique, prodigal_unique, on='SeqID', 
                             suffixes=('_metagenemark', '_prodigal'), how='inner')
        
        concomitantiam_metagenemark = merged_mgm[
            ((merged_mgm['Start_metagenemark'] < merged_mgm['Start_prodigal']) & 
             (merged_mgm['End_metagenemark'] > merged_mgm['Start_prodigal']) & 
             (merged_mgm['End_metagenemark'] < merged_mgm['End_prodigal'])) |
            ((merged_mgm['Start_metagenemark'] > merged_mgm['Start_prodigal']) & 
             (merged_mgm['Start_metagenemark'] < merged_mgm['End_prodigal']) & 
             (merged_mgm['End_metagenemark'] > merged_mgm['End_prodigal'])) |
            ((merged_mgm['Start_metagenemark'] < merged_mgm['Start_prodigal']) & 
             (merged_mgm['End_metagenemark'] > merged_mgm['End_prodigal']))
        ]
        
        if not concomitantiam_metagenemark.empty:
            concomitantiam_metagenemark = concomitantiam_metagenemark.drop(columns=concomitantiam_metagenemark.filter(like='_prodigal').columns)
            concomitantiam_metagenemark.columns = GFFParser.GFF_COLUMNS
        
        # Prodigal concomitantiam (overlapping with MetaGeneMark genes)
        merged_prod = pd.merge(prodigal_unique, metagenemark_unique, on='SeqID', 
                              suffixes=('_prodigal', '_metagenemark'), how='inner')
        
        concomitantiam_prodigal = merged_prod[
            ((merged_prod['Start_prodigal'] < merged_prod['Start_metagenemark']) & 
             (merged_prod['End_prodigal'] > merged_prod['Start_metagenemark']) & 
             (merged_prod['End_prodigal'] < merged_prod['End_metagenemark'])) |
            ((merged_prod['Start_prodigal'] > merged_prod['Start_metagenemark']) & 
             (merged_prod['Start_prodigal'] < merged_prod['End_metagenemark']) & 
             (merged_prod['End_prodigal'] > merged_prod['End_metagenemark'])) |
            ((merged_prod['Start_prodigal'] < merged_prod['Start_metagenemark']) & 
             (merged_prod['End_prodigal'] > merged_prod['End_metagenemark']))
        ]
        
        if not concomitantiam_prodigal.empty:
            concomitantiam_prodigal = concomitantiam_prodigal.drop(columns=concomitantiam_prodigal.filter(like='_metagenemark').columns)
            concomitantiam_prodigal.columns = GFFParser.GFF_COLUMNS
        
        return concomitantiam_metagenemark, concomitantiam_prodigal
    
    def get_exclusivus_genes(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Find genes that are completely unique (no overlap with other tool's genes)."""
        metagenemark_unique, prodigal_unique = self.get_unique_genes()
        insertio_metagenemark, insertio_prodigal = self.get_insertio_genes()
        concomitantiam_metagenemark, concomitantiam_prodigal = self.get_concomitantiam_genes()
        
        # MetaGeneMark exclusivus
        to_subtract_mgm = pd.concat([insertio_metagenemark, concomitantiam_metagenemark], ignore_index=True)
        merged_mgm = pd.merge(metagenemark_unique, to_subtract_mgm, how='outer', indicator=True)
        exclusivus_metagenemark = merged_mgm[merged_mgm['_merge'] == 'left_only'].drop(columns=['_merge'])
        
        # Prodigal exclusivus
        to_subtract_prod = pd.concat([insertio_prodigal, concomitantiam_prodigal], ignore_index=True)
        merged_prod = pd.merge(prodigal_unique, to_subtract_prod, how='outer', indicator=True)
        exclusivus_prodigal = merged_prod[merged_prod['_merge'] == 'left_only'].drop(columns=['_merge'])
        
        return exclusivus_metagenemark, exclusivus_prodigal


class GFFWriter:
    """Handles writing GFF files."""
    
    @staticmethod
    def write_gff_file(file_path: str, dataframe: pd.DataFrame):
        """Write a DataFrame to a GFF file."""
        try:
            with open(file_path, 'w') as file:
                file.write('##gff-version 3\n')
                if not dataframe.empty:
                    dataframe.to_csv(file, sep='\t', header=False, index=False)
        except Exception as e:
            print(f"Error writing file {file_path}: {e}")
    
    @staticmethod
    def write_unique_lines(file_path: str, dataframe: pd.DataFrame):
        """Write unique lines from DataFrame to a GFF file."""
        try:
            with open(file_path, 'w') as file:
                file.write('##gff-version 3\n')
                lines_set = set()
                for _, row in dataframe.iterrows():
                    line = '\t'.join(row.astype(str)) + '\n'
                    if line not in lines_set:
                        file.write(line)
                        lines_set.add(line)
        except Exception as e:
            print(f"Error writing file {file_path}: {e}")


class GCGTool:
    """Main class for the Gene Comparison and Classification tool."""
    
    def __init__(self, metagenemark_file: str, prodigal_file: str, output_path: str):
        """Initialize the tool with input files and output path."""
        self.metagenemark_file = metagenemark_file
        self.prodigal_file = prodigal_file
        self.output_path = output_path
        self.arquivo_name = os.path.basename(os.path.normpath(output_path))
        
        # Create output directory
        os.makedirs(output_path, exist_ok=True)
        
        # Parse input files
        self.metagenemark_df = GFFParser.parse_gff_file(metagenemark_file)
        self.prodigal_df = GFFParser.parse_gff_file(prodigal_file)
        
        # Initialize classifier
        self.classifier = GeneClassifier(self.metagenemark_df, self.prodigal_df)
    
    def generate_consensus_files(self):
        """Generate consensus gene files."""
        consensus_metagenemark, consensus_prodigal = self.classifier.get_consensus_genes()
        
        GFFWriter.write_gff_file(
            os.path.join(self.output_path, f'{self.arquivo_name}_consensus_metagenemark.gff'),
            consensus_metagenemark
        )
        
        GFFWriter.write_gff_file(
            os.path.join(self.output_path, f'{self.arquivo_name}_consensus_prodigal.gff'),
            consensus_prodigal
        )
        
        print(f"Consensus MetaGeneMark genes: {len(consensus_metagenemark)}")
        print(f"Consensus Prodigal genes: {len(consensus_prodigal)}")
    
    def generate_unique_files(self):
        """Generate unique gene files."""
        metagenemark_unique, prodigal_unique = self.classifier.get_unique_genes()
        
        GFFWriter.write_gff_file(
            os.path.join(self.output_path, f'{self.arquivo_name}_metagenemark_notmerged.gff'),
            metagenemark_unique
        )
        
        GFFWriter.write_gff_file(
            os.path.join(self.output_path, f'{self.arquivo_name}_prodigal_notmerged.gff'),
            prodigal_unique
        )
        
        print(f"Unique MetaGeneMark genes: {len(metagenemark_unique)}")
        print(f"Unique Prodigal genes: {len(prodigal_unique)}")
    
    def generate_insertio_files(self):
        """Generate insertio gene files."""
        insertio_metagenemark, insertio_prodigal = self.classifier.get_insertio_genes()
        
        GFFWriter.write_gff_file(
            os.path.join(self.output_path, f'{self.arquivo_name}_insertio.gff'),
            pd.concat([insertio_metagenemark, insertio_prodigal], ignore_index=True)
        )
        
        print(f"Insertio MetaGeneMark genes: {len(insertio_metagenemark)}")
        print(f"Insertio Prodigal genes: {len(insertio_prodigal)}")
    
    def generate_concomitantiam_files(self):
        """Generate concomitantiam gene files."""
        concomitantiam_metagenemark, concomitantiam_prodigal = self.classifier.get_concomitantiam_genes()
        
        GFFWriter.write_gff_file(
            os.path.join(self.output_path, f'{self.arquivo_name}_concomitantiam.gff'),
            pd.concat([concomitantiam_metagenemark, concomitantiam_prodigal], ignore_index=True)
        )
        
        print(f"Concomitantiam MetaGeneMark genes: {len(concomitantiam_metagenemark)}")
        print(f"Concomitantiam Prodigal genes: {len(concomitantiam_prodigal)}")
    
    def generate_exclusivus_files(self):
        """Generate exclusivus gene files."""
        exclusivus_metagenemark, exclusivus_prodigal = self.classifier.get_exclusivus_genes()
        
        GFFWriter.write_gff_file(
            os.path.join(self.output_path, f'{self.arquivo_name}_exclusivus.gff'),
            pd.concat([exclusivus_metagenemark, exclusivus_prodigal], ignore_index=True)
        )
        
        print(f"Exclusivus MetaGeneMark genes: {len(exclusivus_metagenemark)}")
        print(f"Exclusivus Prodigal genes: {len(exclusivus_prodigal)}")
    
    def run(self, generate_insertio: bool = False, generate_concomitantiam: bool = False, 
            generate_exclusivus: bool = False):
        """Run the complete analysis."""
        print(f"Processing files:")
        print(f"  MetaGeneMark: {self.metagenemark_file}")
        print(f"  Prodigal: {self.prodigal_file}")
        print(f"  Output: {self.output_path}")
        print()
        
        # Always generate consensus and unique files
        self.generate_consensus_files()
        self.generate_unique_files()
        
        # Generate optional files based on flags
        if generate_insertio:
            self.generate_insertio_files()
        
        if generate_concomitantiam:
            self.generate_concomitantiam_files()
        
        if generate_exclusivus:
            self.generate_exclusivus_files()
        
        print("\nAnalysis complete!")


def main():
    """Main function to handle command line arguments and run the tool."""
    parser = OptionParser(usage="%prog -m METAGENEMARK_FILE -p PRODIGAL_FILE -o OUTPUT_PATH [options]")
    
    parser.add_option("-m", "--metagenemark", dest="metagenemark_file",
                      help="MetaGeneMark GFF file", metavar="FILE")
    
    parser.add_option("-p", "--prodigal", dest="prodigal_file",
                      help="Prodigal GFF file", metavar="FILE")
    
    parser.add_option("-i", "--insertio", action="store_true", dest="generate_insertio", default=False,
                      help="Generate files containing inserted genes")
    
    parser.add_option("-c", "--concomitantiam", action="store_true", dest="generate_concomitantiam", default=False,
                      help="Generate files containing semi-overlapping genes")
    
    parser.add_option("-e", "--exclusivus", action="store_true", dest="generate_exclusivus", default=False,
                      help="Generate files containing exclusive genes")
    
    parser.add_option("-o", "--output", dest="output_path",
                      help="Output directory path", metavar="DIR")
    
    (options, args) = parser.parse_args()
    
    # Validate required arguments
    if not options.metagenemark_file or not options.prodigal_file or not options.output_path:
        parser.print_help()
        sys.exit(1)
    
    # Check if input files exist
    if not os.path.exists(options.metagenemark_file):
        print(f"Error: MetaGeneMark file '{options.metagenemark_file}' not found.")
        sys.exit(1)
    
    if not os.path.exists(options.prodigal_file):
        print(f"Error: Prodigal file '{options.prodigal_file}' not found.")
        sys.exit(1)
    
    try:
        # Initialize and run the tool
        tool = GCGTool(options.metagenemark_file, options.prodigal_file, options.output_path)
        tool.run(
            generate_insertio=options.generate_insertio,
            generate_concomitantiam=options.generate_concomitantiam,
            generate_exclusivus=options.generate_exclusivus
        )
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
