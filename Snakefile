from snakemake.utils import R
import csv
import glob
import gzip
import itertools
import os
from pathlib import Path
import re

# Define SEPINC as shell environment variable
SNAKERULES = os.environ.get('SNAKERULES')
PROJ_NAME = "genomes"
DEPDIR = "deps"
include:
    SNAKERULES + "/Doc.defs.top"

rule bowtie2_index_GRCh37Sponge:
    input: genome=DATA_PATH_GALE+'/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa',\
      sponge=DATA_PATH_GALE+'/human/Miga_Bioinformatids2015_hgsponge/GRCh37_sponge.fa'
    params: D=BOWTIE2_v2p2p9_INDEXES_PATH +'/GRCh37_sponge',bt2_base='GRCh37_sponge'
    output: BOWTIE2_v2p2p9_INDEXES_PATH +'/GRCh37_sponge/GRCh37_sponge.1.bt2'
    shell: """
        module load bowtie2/2.2.9bin
        cd {params.D}; bowtie2-build {input.genome},{input.sponge} {params.bt2_base}
    """

rule bowtie2_index_GRCh38Sponge:
    input: genome=DATA_PATH_GALE+'/Illumina_iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa',\
      sponge=DATA_PATH_GALE+'/human/Miga_Bioinformatids2015_hgsponge/GRCh38_sponge.fa'
    params: D=BOWTIE2_v2p2p9_INDEXES_PATH +'/GRCh38_sponge',bt2_base='GRCh38_sponge'
    output: BOWTIE2_v2p2p9_INDEXES_PATH +'/GRCh38_sponge/GRCh38_sponge.1.bt2'
    shell: """
        module load bowtie2/2.2.9bin
        cd {params.D}; bowtie2-build {input.genome},{input.sponge} {params.bt2_base}
    """

At_6genomes = {'CS70000':'Col-0_pchr.fasta','6911':'Cvi-0_pchr.fasta','LEHLE':'Ler-0_pchr.fasta',\
               '7063':'Can-0_pchr.fasta','5784':'Ty-1_pchr.fasta','1741':'KBS-Mac-74_pchr.fasta'}
rule bowtie2_index_At6g:
    input: genome=PROJ_DATA_PATH_GALE+'/At_6genomes/{tg_ecotypeid}/Sequence/WholeGenomeFasta/genome.fa'
    params: D=PROJ_DATA_PATH_GALE+'/At_6genomes/{tg_ecotypeid}/Sequence/Bowtie2Index',bt2_base="{tg_ecotypeid}"
    output: PROJ_DATA_PATH_GALE+'/At_6genomes/{tg_ecotypeid}/Sequence/Bowtie2Index/{tg_ecotypeid}.1.bt2'
    shell: """
        module load bowtie2/2.2.9bin
        mkdir -p {params.D}; cd {params.D}; bowtie2-build {input.genome} {params.bt2_base}
    """
rule bowtie2_index_At6g_out:
    input: expand(PROJ_DATA_PATH_GALE+'/At_6genomes/{tg_ecotypeid}/Sequence/Bowtie2Index/{tg_ecotypeid}.1.bt2',\
                  tg_ecotypeid=At_6genomes.keys())
    
rule bowtie2_index_Esa1:
    input: genome=PROJ_DATA_PATH_GALE +'/Esa1.0/Sequence/WholeGenomeFasta/genome_chr.fa'
    params: D=BOWTIE2_v2p2p9_INDEXES_PATH +'/Esa1.0',bt2_base='Esa1.0'
    output: BOWTIE2_v2p2p9_INDEXES_PATH +'/Esa1.0/Esa1.0.1.bt2'
    shell: """
        module load bowtie2/2.2.9bin
        cd {params.D}; bowtie2-build {input.genome} {params.bt2_base}
    """

rule wg_fasta_GRCh37Sponge:
    input: genome=DATA_PATH_GALE+'/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa',\
      sponge=DATA_PATH_GALE+'/human/Miga_Bioinformatids2015_hgsponge/GRCh37_sponge.fa'
    output: PROJ_DATA_PATH_GALE +'/GRCh37_sponge/Sequence/WholeGenomeFasta/genome.fa'
    shell: """
        cat {input.genome} {input.sponge} > {output}
    """

rule wg_fasta_GRCh38Sponge:
    input: genome=DATA_PATH_GALE+'/Illumina_iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa',\
      sponge=DATA_PATH_GALE+'/human/Miga_Bioinformatids2015_hgsponge/GRCh38_sponge.fa'
    output: PROJ_DATA_PATH_GALE +'/GRCh38_sponge/Sequence/WholeGenomeFasta/genome.fa'
    shell: """
        cat {input.genome} {input.sponge} > {output}
    """

rule wg_fasta_At6:
    input: nuc=lambda wc: "/gale/netapp/seq4/illumina_runs/florian_6genome/"+At_6genomes[wc.tg_ecotypeid],\
           pt="/gale/netapp/seq4/illumina_runs/florian_6genome/At_chloroplast_TAIR.fasta",\
           mt="/gale/netapp/seq4/illumina_runs/florian_6genome/At_mitochondria_TAIR.fasta"
    output: PROJ_DATA_PATH_GALE+'/At_6genomes/{tg_ecotypeid}/Sequence/WholeGenomeFasta/genome.fa'
    shell: """
        sed -e "s/^>chr_/>/" {input.nuc} > {output}
        sed "1s/.*/>Pt/" {input.pt} | fold -w 60 >> {output}
        sed "1s/.*/>Mt/" {input.mt} | fold -w 60 >> {output}
    """
rule wg_fasta_At6g_out:
    input: expand(PROJ_DATA_PATH_GALE+'/At_6genomes/{tg_ecotypeid}/Sequence/WholeGenomeFasta/genome.fa',tg_ecotypeid=At_6genomes.keys())
    
rule wg_fasta_chr_GRCh37:
    input: fa=PROJ_DATA_PATH_GALE +'/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa'
    output: fa=PROJ_DATA_PATH_GALE +'/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome_chr.fa'
    shell: """
        sed -e "s/^>\([0-9][0-9]\?\|X\|Y\|MT\)/>chr\\1/" {input.fa} > {output.fa}         
    """
rule wg_fasta_chr_GRCh37Sponge:
    input: fa=PROJ_DATA_PATH_GALE +'/GRCh37_sponge/Sequence/WholeGenomeFasta/genome.fa'
    output: fa=PROJ_DATA_PATH_GALE +'/GRCh37_sponge/Sequence/WholeGenomeFasta/genome_chr.fa'
    shell: """
        sed -e "s/^>\([0-9][0-9]\?\|X\|Y\|MT\)/>chr\\1/" {input.fa} > {output.fa}         
    """
rule wg_fasta_chr_TAIR10:
    input: fa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome.fa'
    output: fa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome_chr.fa'
    shell: """
        sed -e "s/^>\([0-9][0-9]\?\|Mt\|Pt\)/>chr\\1/" {input.fa} > {output.fa}         
    """
rule wg_fasta_chr_At6g:
    input: fa=PROJ_DATA_PATH_GALE+'/At_6genomes/{tg_ecotypeid}/Sequence/WholeGenomeFasta/genome.fa'
    output: fa=PROJ_DATA_PATH_GALE+'/At_6genomes/{tg_ecotypeid}/Sequence/WholeGenomeFasta/genome_chr.fa'
    shell: """
        sed -e "s/^>\([0-9][0-9]\?\|Mt\|Pt\)/>chr\\1/" {input.fa} > {output.fa}       
    """
rule wg_fasta_chr_At6g_out:
    input: expand(PROJ_DATA_PATH_GALE+'/At_6genomes/{tg_ecotypeid}/Sequence/WholeGenomeFasta/genome_chr.fa',tg_ecotypeid=At_6genomes.keys())
    
rule wg_fasta_chr_Esa1:
    input: fa=PROJ_DATA_PATH_GALE +'/Esa1.0/Sequence/WholeGenomeFasta/Esalsugineum_genome_v1.0.fa'
    output: genome_fa=PROJ_DATA_PATH_GALE +'/Esa1.0/Sequence/WholeGenomeFasta/genome.fa',\
            genome_chr_fa=PROJ_DATA_PATH_GALE +'/Esa1.0/Sequence/WholeGenomeFasta/genome_chr.fa'
    shell: """
        sed -e "s/^>ch\([0-9]\+\)/>chr\\1/" {input.fa} > {output.genome_chr_fa}
        ln -s {output.genome_chr_fa} {output.genome_fa}
    """

rule wg_fasta_chr_and_lambda_TAIR10:
    input: gfa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome_chr.fa',lfa=PROJ_DATA_PATH_GALE +'/misc/NC001416.fa'
    output: gfa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WG_lambda_Fasta/genome_chr.fa',lfa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WG_lambda_Fasta/NC001416.fa'
    shell: """
        cp {input.gfa} {output.gfa}
        cp {input.lfa} {output.gfa}
    """
rule wg_fasta_GRCh37_TAIR10_lambda:
    input: GRCh37=DATA_PATH_GALE+'/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa',TAIR10=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome.fa',lamb=PROJ_DATA_PATH_GALE +'/misc/NC001416.fa'
    output: PROJ_DATA_PATH_GALE+'/mixed/GRCh37_TAIR10_Lambda/Sequence/WholeGenomeFasta/genome.fa'
    shell: """
        sed -e "s/^>\(.*\)/>GRCh37_Ensembl_\\1/" {input.GRCh37} > {output}
        sed -e "s/^>\(.*\)/>TAIR10_Ensembl_\\1/" {input.TAIR10} >> {output}
        sed -e "s/^>\(.*\)/>Lambda_NCBI_L/" {input.lamb} >> {output}
    """
rule wg_fasta_mm10_lambda:
    input: mm10=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa',\
           lamb=PROJ_DATA_PATH_GALE +'/misc/NC001416.fa'
    output: PROJ_DATA_PATH_GALE+'/mixed/mm10_lambda/Sequence/WholeGenomeFasta/genome.fa'
    shell: """
        cp {input.mm10} {output}
        sed -e "s/^>\(.*\)/>chrL/" {input.lamb} >> {output}
    """

rule nuc_chr_At6g:
    input: PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/SeqInfo/nuclear_chromosomes.txt'
    output: PROJ_DATA_PATH_GALE +'/At_6genomes/{tg_ecotypeid}/Sequence/SeqInfo/nuclear_chromosomes.txt'
    shell: """
        cp {input} {output}
    """
rule nuc_chr_At6g_out:
    input: expand(PROJ_DATA_PATH_GALE +'/At_6genomes/{tg_ecotypeid}/Sequence/SeqInfo/nuclear_chromosomes.txt',tg_ecotypeid=At_6genomes.keys())
    
rule wg_fasta_chr_AGPv3:
    input: fa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Zea_mays/Ensembl/AGPv3/Sequence/WholeGenomeFasta/genome.fa'
    output: fa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Zea_mays/Ensembl/AGPv3/Sequence/WholeGenomeFasta/genome_chr.fa'
    shell: """
        sed -e "s/^>\([0-9][0-9]\?\|Mt\|Pt\)/>chr\\1/" {input.fa} > {output.fa}         
    """
rule wg_fasta_chr_AGPv4:
    input: fa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Zea_mays/Ensembl/AGPv4/Sequence/WholeGenomeFasta/genome.fa'
    output: fa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Zea_mays/Ensembl/AGPv4/Sequence/WholeGenomeFasta/genome_chr.fa'
    shell: """
        sed -e "s/^>\([0-9][0-9]\?\|Mt\|Pt\)/>chr\\1/" {input.fa} > {output.fa}         
    """
rule nuc_fasta_chr_TAIR10:
    input: fa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/nuclear_chromosomes.fa'
    output: fa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/nuclear_chromosomes_chr.fa'
    shell: """
        sed -e "s/^>\([0-9][0-9]\?\|Mt\|Pt\)/>chr\\1/" {input.fa} > {output.fa}         
    """
rule nuc_faidx_chr_TAIR10:
    input: fa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/nuclear_chromosomes_chr.fa'
    output: fai=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/nuclear_chromosomes_chr.fa.fai'
    shell: """
        module load samtools/1.2
        samtools faidx {input.fa}
    """

rule chr_fasta_GRCh37Sponge:# one fasta per chromosome
    input: fa=PROJ_DATA_PATH_GALE +'/GRCh37_sponge/Sequence/WholeGenomeFasta/genome.fa'
    params: D=PROJ_DATA_PATH_GALE +'/GRCh37_sponge/Sequence/Chromosomes'
    output: PROJ_DATA_PATH_GALE +'/GRCh37_sponge/Sequence/Chromosomes/1.fa'
    run:
        f=open(input.fa,"r");
        opened = False
        for line in f :
            if(line[0] == ">") :
                if(opened) :
                    of.close()
                opened = True
                contig_name = line[1:].rstrip()
                ofname = "%s/%s.fa" % (params.D,contig_name) 
                of=open(ofname, "w")
                print("contig: %s; output file %s"%(contig_name,ofname))
            of.write(line)
        of.close()
rule chr_fasta_GRCh38Sponge:# one fasta per chromosome
    input: fa=PROJ_DATA_PATH_GALE +'/GRCh38_sponge/Sequence/WholeGenomeFasta/genome.fa'
    params: D=PROJ_DATA_PATH_GALE +'/GRCh38_sponge/Sequence/Chromosomes'
    output: PROJ_DATA_PATH_GALE +'/GRCh38_sponge/Sequence/Chromosomes/chr1.fa'
    run:
        f=open(input.fa,"r");
        opened = False
        for line in f :
            if(line[0] == ">") :
                if(opened) :
                    of.close()
                opened = True
                contig_name = line[1:].rstrip().split()[0]
                ofname = "%s/%s.fa" % (params.D,contig_name) 
                of=open(ofname, "w")
                print("contig: %s; output file %s"%(contig_name,ofname))
            of.write(line)
        of.close()
rule chr_fasta_At6g:# one fasta per chromosome
    input: fa=PROJ_DATA_PATH_GALE +'/At_6genomes/{tg_ecotypeid}/Sequence/WholeGenomeFasta/genome.fa'
    params: D=PROJ_DATA_PATH_GALE +'/At_6genomes/{tg_ecotypeid}/Sequence/Chromosomes'
    output: PROJ_DATA_PATH_GALE +'/At_6genomes/{tg_ecotypeid}/Sequence/Chromosomes/1.fa'
    run:
        f=open(input.fa,"r");
        opened = False
        for line in f :
            if(line[0] == ">") :
                if(opened) :
                    of.close()
                opened = True
                contig_name = line[1:].rstrip()
                ofname = "%s/%s.fa" % (params.D,contig_name) 
                of=open(ofname, "w")
                print("contig: %s; output file %s"%(contig_name,ofname))
            of.write(line)
        of.close()
rule chr_fasta_At6g_out:
    input: expand(PROJ_DATA_PATH_GALE+'/At_6genomes/{tg_ecotypeid}/Sequence/Chromosomes/1.fa',tg_ecotypeid=At_6genomes.keys())
rule chr_fasta_Esa1:# one fasta per chromosome
    input: fa=PROJ_DATA_PATH_GALE +'/Esa1.0/Sequence/WholeGenomeFasta/genome_chr.fa'
    params: D=PROJ_DATA_PATH_GALE +'/Esa1.0/Sequence/Chromosomes'
    output: PROJ_DATA_PATH_GALE +'/Esa1.0/Sequence/Chromosomes/chr1-01.fa'
    run:
        f=open(input.fa,"r");
        opened = False
        for line in f :
            if(line[0] == ">") :
                if(opened) :
                    of.close()
                opened = True
                contig_name = line[1:].rstrip().split()[0]
                ofname = "%s/%s.fa" % (params.D,contig_name) 
                of=open(ofname, "w")
                print("contig: %s; output file %s"%(contig_name,ofname))
            of.write(line)
        of.close()

    
rule add_chr_prefix:# one fasta per chromosome, file names start with "chr" to be compatible with GEM requirement
    input: fa="{seqence_dir}/Chromosomes/{chr_num}.fa"
    output: fa="{seqence_dir}/Chromosomes_CHR/chr{chr_num}.fa"
    shell: """
       ln -s {input.fa} {output.fa}
    """
rule add_chr_prefix_GRCh37:
    input: expand(PROJ_DATA_PATH_GALE +'/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Chromosomes_CHR/chr{chr_num}.fa',chr_num=[str(x) for x in range(1,23)]+['X','Y','MT'])
rule add_chr_prefix_GRCh37Sponge:
    input: expand(PROJ_DATA_PATH_GALE +'/GRCh37_sponge/Sequence/Chromosomes_CHR/chr{chr_num}.fa',chr_num=[str(x) for x in range(1,23)]+['X','Y','MT'])
rule add_chr_prefix_AGPv3:
    input: expand(PROJ_DATA_PATH_GALE +'/Illumina_iGenomes/Zea_mays/Ensembl/AGPv3/Sequence/Chromosomes_CHR/chr{chr_num}.fa',chr_num=[str(x) for x in range(1,11)]+['MT','Pt'])
rule add_chr_prefix_AGPv4:
    input: expand(PROJ_DATA_PATH_GALE +'/Illumina_iGenomes/Zea_mays/Ensembl/AGPv4/Sequence/Chromosomes_CHR/chr{chr_num}.fa',chr_num=[str(x) for x in range(1,11)]+['MT','Pt'])
rule add_chr_prefix_TAIR10:
    input: expand(PROJ_DATA_PATH_GALE +'/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Chromosomes_CHR/chr{chr_num}.fa',chr_num=[str(x) for x in range(1,6)]+['Mt','Pt'])
rule add_chr_prefix_At6g:
    input: expand(PROJ_DATA_PATH_GALE +'/At_6genomes/{tg_ecotypeid}/Sequence/Chromosomes_CHR/chr{chr_num}.fa',\
                  tg_ecotypeid=At_6genomes.keys(),\
                  chr_num=[str(x) for x in range(1,6)]+['Mt','Pt'])

rule clean_fa_header:# one fasta per chromosome, clean up the FASTA file headers to keep only sequence names
    input: fa="{seqence_dir}/Chromosomes/chr{chr_num}.fa"
    output: fa="{seqence_dir}/Chromosomes_CHR/chr{chr_num}.fa"
    shell: """
        sed -e "s/^\(>[^[:space:]]*\).*/\\1/" {input.fa} > {output.fa}
    """
rule clean_fa_header_GRCh38:
   input: f.replace('/Chromosomes/','/Chromosomes_CHR/') for f in glob.glob(PROJ_DATA_PATH_GALE +'/Illumina_iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/Chromosomes/chr*.fa')
rule clean_fa_header_GRCh38Sponge:
   input: f.replace('/Chromosomes/','/Chromosomes_CHR/') for f in glob.glob(PROJ_DATA_PATH_GALE +'/GRCh38_sponge/Sequence/Chromosomes/chr*.fa')
rule clean_fa_header_GRCh38Decoy:
   input: f.replace('/Chromosomes/','/Chromosomes_CHR/') for f in glob.glob(PROJ_DATA_PATH_GALE +'/Illumina_iGenomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/Chromosomes/chr*.fa')

rule gatk_ref_TAIR10:
   input: PROJ_DATA_PATH_GALE+'/TAIR10/Sequences/whole_chromosomes/col0_tair10_chr.fasta'
   output: PROJ_DATA_PATH_GALE+'/TAIR10/Sequences/whole_chromosomes/col0_tair10_chr.dict'
   shell: """
        module load shhuang java/jdk-8u102-linux-64; java -jar {PICARD_2p8_JAR} CreateSequenceDictionary R={input} O={output}
   """
   
rule samtools_faidx:
    input: "{genome}.fa"
    output: "{genome}.fa.fai"
    shell: """
        module load samtools/1.2
        samtools faidx {input}
    """
rule samtools_nuc_chr:
    input: fa="{sequence_dir}/WholeGenomeFasta/genome.fa",nuc_chr="{sequence_dir}/SeqInfo/nuclear_chromosomes.txt"
    output: nuc_fa="{sequence_dir}/WholeGenomeFasta/nuclear_chromosomes.fa"
    shell: """
        module load samtools/1.2
        samtools faidx {input.fa} `cat {input.nuc_chr}` > {output.nuc_fa}
    """
rule nuc_chr_faidx:
    input: expand(PROJ_DATA_PATH_GALE+"{genome_dir}/nuclear_chromosomes.fa.fai",\
                  genome_dir=['/GRCh37_sponge/Sequence/WholeGenomeFasta',\
                            '/GRCh38_sponge/Sequence/WholeGenomeFasta',\
                            '/Esa1.0/Sequence/WholeGenomeFasta',\
                            '/Illumina_iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta',\
                            '/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta',\
                            '/Illumina_iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta',\
                            '/Illumina_iGenomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta',\
                            '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta'\
                      ] + expand('/At_6genomes/{tg_ecotypeid}/Sequence/WholeGenomeFasta',tg_ecotypeid=At_6genomes.keys()))

rule samtools_chr1:
    input: fa="{sequence_dir}/WholeGenomeFasta/genome.fa",chr1="{sequence_dir}/SeqInfo/chr1.txt"
    output: chr1_fa="{sequence_dir}/WholeGenomeFasta/chr1.fa"
    shell: """
        module load samtools/1.2
        samtools faidx {input.fa} `cat {input.chr1}` > {output.chr1_fa}
    """
rule chr1_faidx:
    input: expand(PROJ_DATA_PATH_GALE+"{genome_dir}/chr1.fa.fai",\
                  genome_dir=['/GRCh37_sponge/Sequence/WholeGenomeFasta',\
                            '/GRCh38_sponge/Sequence/WholeGenomeFasta',\
                            '/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta',\
                            '/Illumina_iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta',\
                            '/Illumina_iGenomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta'\
                      ])
    
rule get_chrom_sizes:
	input: "{genome}.fa.fai"
	output: "{genome}.chrom.sizes"
	shell:"""
        cut -f1,2 {input} > {output}
    """

rule get_chrom_sizes_out:
    input: expand(PROJ_DATA_PATH_GALE+"{genome_dir}/{g}.chrom.sizes",\
                  g=['genome','nuclear_chromosomes'],\
                  genome_dir=['/GRCh37_sponge/Sequence/WholeGenomeFasta',\
                        '/GRCh38_sponge/Sequence/WholeGenomeFasta',\
                        '/Esa1.0/Sequence/WholeGenomeFasta',\
                        '/Illumina_iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta',\
                        '/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta',\
                        '/Illumina_iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta',\
                        '/Illumina_iGenomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta',\
                        '/Illumina_iGenomes/Zea_mays/Ensembl/AGPv3/Sequence/WholeGenomeFasta',\
                        '/Illumina_iGenomes/Zea_mays/Ensembl/AGPv4/Sequence/WholeGenomeFasta',\
                        '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta',\
                        '/Illumina_iGenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta',\
                        '/mixed/GRCh37_TAIR10_Lambda/Sequence/WholeGenomeFasta'\
                      ]) +
            expand(PROJ_DATA_PATH_GALE+"{genome_dir}/{g}.chrom.sizes",\
                  g=['nuclear_chromosomes_chr'],\
                  genome_dir=['/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta'
                      ]) +
            expand(PROJ_DATA_PATH_GALE+"/At_6genomes/{tg_ecotypeid}/Sequence/WholeGenomeFasta/{g}.chrom.sizes",\
                   g=['genome','nuclear_chromosomes'],tg_ecotypeid=At_6genomes.keys())


rule get_chrom_sizes_bed:
    input: "{genome_dir}/Sequence/WholeGenomeFasta/genome.chrom.sizes"
    output: "{genome_dir}/Annotation/Genes/genome.bed"
    shell: """
         awk -v FS='\\t' -v OFS='\\'t '{{ print $1,0,$2 }}' {input} > {output}
    """
rule get_chrom_sizes_bed_out:
    input: expand(PROJ_DATA_PATH_GALE+"{genome_dir}/Annotation/Genes/genome.bed",\
                  genome_dir=["/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10",\
                              "/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37"])
# makewindows
WSIZE_DICT = {'10kb':'10000','100kb':'100000'}
rule make_windows:
    input: "{genome_dir}/Sequence/WholeGenomeFasta/nuclear_chromosomes.chrom.sizes"
    params: w=lambda wc: WSIZE_DICT[wc.wsize]
    output: tsv="{genome_dir}/Sequence/SeqInfo/nuclear_chromosomes.w{wsize}.tsv",\
            bed="{genome_dir}/Sequence/SeqInfo/nuclear_chromosomes.w{wsize}.bed"
    shell: """
        module load bedtools/2.26.0
        bedtools makewindows -g {input} -w {params.w} > {output.bed}
        echo -e 'chrom\\tstart\\tend' > {output.tsv}
        awk -v FS='\\t' -v OFS='\\t' '{{ print $1,$2+1,$3 }}' {output.bed} >> {output.tsv}
    """
rule make_windows_out:
    input: expand(PROJ_DATA_PATH_GALE+"{genome_dir}/Sequence/SeqInfo/nuclear_chromosomes.w{wsize}.tsv",\
                  wsize=['10kb','100kb'],\
                  genome_dir=['/mixed/GRCh37_TAIR10_Lambda'])

# SnpEff
rule snpeff_build_dapv1:
    input: DAP_DISTR_PATH + '/Ath_DAPv1_SnpEff.tar'
    log: SNPEFF_4p1L_INSTALL_PATH + '/data/athalianaTair10/Ath_DAPv1_SnpEff.build.log'
    output: SNPEFF_4p1L_INSTALL_PATH + '/data/athalianaTair10/Ath_DAPv1_SnpEff.build.done'
    shell: """
        module load java/jdk-8u102-linux-64
        mkdir -p {SNPEFF_4p1L_INSTALL_PATH}/data/athalianaTair10/regulation.bed
        tar -C {SNPEFF_4p1L_INSTALL_PATH}/data/athalianaTair10/regulation.bed -xvf {input}
        java -Xmx20G -jar {SNPEFF_4p1L_JAR} build -v -onlyReg athalianaTair10 > {log} 2>&1
    """
    
# methylomes
rule mod_genome_hg19:
    input: allc_h5=HS_METH_PATH + '/hg19/{sample}/allc/allc_{sample}.h5',\
           genome=PROJ_DATA_PATH_GALE+'/Illumina_iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/nuclear_chromosomes.fa'
    output: mod_fa=HS_METH_PATH + '/hg19/{sample}/mod_fasta/modgenome_{sample}.fa'
    shell: """
        unset PYTHONPATH; module load shhuang anaconda/2.4.1; python {LIB_PYTHON_UTILS_PATH}/get_mod_genome.py {input.genome} {input.allc_h5} {output.mod_fa}
    """
rule mod_genome_hg19_out:
    input: expand(HS_METH_PATH + '/hg19/{sample}/mod_fasta/modgenome_{sample}.fa',sample='h1')

# methylomes
rule mod_genome_at1001:
    input: allc_h5=AT1001_ALLC + '/allc_{sample}.h5',\
           genome=PROJ_DATA_PATH_GALE+'/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/nuclear_chromosomes_chr.fa'
    output: mod_fa=AT1001_MODFAS + '/modgenome_{sample}.fa'
    shell: """
        unset PYTHONPATH; module load shhuang anaconda/2.4.1; python {LIB_PYTHON_UTILS_PATH}/get_mod_genome.py --zero {input.genome} {input.allc_h5} {output.mod_fa}
    """
rule mod_genome_at1001_out:
    input: expand(AT1001_MODFAS + '/modgenome_{sample}.fa',sample=['6909_pub'])

# MEME background
rule meme_markov_wg:
	input: "{sequence_dir}/WholeGenomeFasta/{genome}.fa"
	output: "{sequence_dir}/SeqInfo/{genome}.m{order}.bg"
	shell:"""
        module load meme/4.11.2; fasta-get-markov -m {wildcards.order} {input} {output}
    """

rule meme_markov_wg_out:
    input: expand(PROJ_DATA_PATH_GALE+"{genome_dir}/SeqInfo/{g}.m{o}.bg",\
                  g=['genome','nuclear_chromosomes'],\
                  o=[0,1,2],\
                  genome_dir=['/GRCh37_sponge/Sequence',\
                              '/GRCh38_sponge/Sequence',\
                              '/Illumina_iGenomes/Homo_sapiens/UCSC/hg19/Sequence',\
                              '/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence',\
                              '/Illumina_iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence',\
                              '/Illumina_iGenomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence',\
                              '/Illumina_iGenomes/Zea_mays/Ensembl/AGPv3/Sequence',\
                              '/Illumina_iGenomes/Zea_mays/Ensembl/AGPv4/Sequence',\
                              '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence'\
                      ])

rule meme_markov_methyl:
    input: HS_METH_PATH + "/{genome}/{sample}/mod_fasta/modgenome_{sample}.fa"
    params: alph=MEME_v4p11p2_DNA_MOD1_ALPH
	output: HS_METH_PATH + "/{genome}/{sample}/mod_fasta/modgenome_{sample}.m{order}.bg"
	shell:"""
        module load meme/4.11.2; fasta-get-markov -m {wildcards.order} -alph {params.alph} {input} {output}
    """

rule meme_markov_methyl_out:
    input: expand(DATA_PATH_GALE+"/{methyl_data}.m{o}.bg",\
                  methyl_data="human/methylcseq/hg19/h1/mod_fasta/modgenome_h1",\
                  o=[0,1,2])

rule meme_markov_at1001_methyl:
    input: AT1001_MODFAS + "/modgenome_{sample}.fa"
    params: alph=MEME_v4p11p2_DNA_MOD1_ALPH
	output: AT1001_MODFAS + "/modgenome_{sample}.m{order}.bg"
	shell:"""
        module load meme/4.11.2; fasta-get-markov -m {wildcards.order} -alph {params.alph} {input} {output}
    """

rule meme_markov_at1001_methyl_out:
    input: expand(DATA_PATH_GALE+"/{methyl_data}.m{o}.bg",\
                  methyl_data="1001_genomes/mod_fasta/modgenome_6909_pub",\
                  o=[0,1,2])

# Bismark
rule bismark_genome_GRCh37_TAIR10_lambda:
    input: PROJ_DATA_PATH_GALE+'/mixed/GRCh37_TAIR10_Lambda/Sequence/WholeGenomeFasta/genome.fa'
    params: D=PROJ_DATA_PATH_GALE+'/mixed/GRCh37_TAIR10_Lambda/Sequence/BismarkGenome'
    output: PROJ_DATA_PATH_GALE+'/mixed/GRCh37_TAIR10_Lambda/Sequence/BismarkGenome/genome.fa'
    shell: """
        module load bowtie2/2.2.9bin Bismark/0.18.1
        mkdir -p {params.D}; cp {input} {output}
        bismark_genome_preparation {params.D}
    """

rule bismark_genome_TAIR10:
    input: fa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome.fa'
    params: D=PROJ_DATA_PATH_GALE+'/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/BismarkGenome'
    output: PROJ_DATA_PATH_GALE+'/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/BismarkGenome/genome.fa'
    shell: """
        module load bowtie2/2.2.9bin Bismark/0.18.1
        mkdir -p {params.D}; cp {input} {output}
        bismark_genome_preparation {params.D}
    """

rule bismark_genome_mm10:
    input: fa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa'
    params: D=PROJ_DATA_PATH_GALE+'/Illumina_iGenomes/Mus_musculus/UCSC/mm10/Sequence/BismarkGenome'
    output: PROJ_DATA_PATH_GALE+'/Illumina_iGenomes/Mus_musculus/UCSC/mm10/Sequence/BismarkGenome/genome.fa'
    shell:    """
        module load bowtie2/2.2.9bin Bismark/0.18.1
        mkdir -p {params.D}; cp {input} {output}
        bismark_genome_preparation {params.D}
    """

rule bismark_genome_mm10_lambda:
    input: PROJ_DATA_PATH_GALE+'/mixed/mm10_lambda/Sequence/WholeGenomeFasta/genome.fa'
    params: D=PROJ_DATA_PATH_GALE+'/mixed/mm10_lambda/Sequence/BismarkGenome'
    output: PROJ_DATA_PATH_GALE+'/mixed/mm10_lambda/Sequence/BismarkGenome/genome.fa'
    shell: """
        module load bowtie2/2.2.9bin Bismark/0.18.1
        mkdir -p {params.D}; cp {input} {output}
        bismark_genome_preparation {params.D}
    """
rule bismark_genome_rn6:
    input: fa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa'
    params: D=PROJ_DATA_PATH_GALE+'/Illumina_iGenomes/Rattus_norvegicus/UCSC/rn6/Sequence/BismarkGenome'
    output: PROJ_DATA_PATH_GALE+'/Illumina_iGenomes/Rattus_norvegicus/UCSC/rn6/Sequence/BismarkGenome/genome.fa'
    shell:    """
        module load bowtie2/2.2.9bin Bismark/0.18.1
        mkdir -p {params.D}; cp {input} {output}
        bismark_genome_preparation {params.D}
    """

    
# Methylpy
rule bismark_genome_GRCh37:
    input: fa=PROJ_DATA_PATH_GALE + '/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa'
    params: pref=PROJ_DATA_PATH_GALE+'/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/MethylpyGenome_Bt2/GRCh37_bt2'
    log: PROJ_DATA_PATH_GALE+'/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/MethylpyGenome_Bt2/GRCh37_bt2.log'
    output: PROJ_DATA_PATH_GALE+'/Illumina_iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/MethylpyGenome_Bt2/GRCh37_bt2_f.4.bt2'
    shell: """
        module purge; unset PYTHONPATH; module load shhuang python/2.7.10 methylpy/1.2.9 scipy/0.16.1 numpy/1.10.2  bowtie2/2.2.9bin samtools/1.2 java/jdk-8u102-linux-64 picard/2.8.3 ucsc_utils/20160210 
        methylpy build-reference --input-files {input.fa} --output-prefix {params.pref} --bowtie2 True
    """

# nucleotide frequencies
rule nuc_freq:
    input: fa="{sequence_dir}/WholeGenomeFasta/genome.fa"
    output: tsv="{sequence_dir}/SeqInfo/genome.nt{d}_freq.tsv"
    run:
        shell("module load r/3.3.1")
        R("""
        library(Biostrings)
        library(plyr)
        seq = readDNAStringSet("{input.fa}")
        frq_f = oligonucleotideFrequency(seq,width={wildcards.d},simplify.as="list")
        names(frq_f) = names(seq)
        frq_f_df = ldply(frq_f,function(si) data.frame(nt=names(si),count=si,strand="+"),.id="chrom")
        frq_r = oligonucleotideFrequency(reverseComplement(seq),width={wildcards.d},simplify.as="list")
        names(frq_r) = names(seq)
        frq_r_df = ldply(frq_r,function(si) data.frame(nt=names(si),count=si,strand="-"),.id="chrom")
        frq_df = rbind(frq_f_df,frq_r_df)
        write.table(frq_df,"{output.tsv}",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
    """)
rule nuc_freq_out:
    input: expand(PROJ_DATA_PATH_GALE+"{genome_dir}/SeqInfo/genome.nt{d}_freq.tsv",\
                  d=[1,2,3],\
                  genome_dir=['/Illumina_iGenomes/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence',\
                              '/mixed/GRCh37_TAIR10_Lambda/Sequence']\
                  )
