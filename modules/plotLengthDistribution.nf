#!/usr/bin/env nextflow

process plotLengthDistribution {
    publishDir 'RESULTS/8_plots'
    container 'simp-test:latest'

    input:
    path input_fasta_file
    
    output:
    path "${input_fasta_file.baseName}_length_distribution.png"

    script:
    """
    #!/opt/conda/bin/python
    
    from Bio import SeqIO
    import pandas as pd 
    import matplotlib.pyplot as plt

    records = list(SeqIO.parse("$input_fasta_file", "fasta"))
    lengths = [len(record.seq) for record in records]
    no_seqs = len(records)
    
    df = pd.DataFrame({"Sequence Length": lengths})
    length_value_counts_df = df.value_counts().reset_index()
    length_value_counts_df.columns = ["Transcript Length", "Count"]

    ax = length_value_counts_df.plot.bar(x='Count', y='Transcript Length', rot=0, color="skyblue")
    ax.set_xlabel("Number of Transcripts")
    ax.set_ylabel("Transcript Length")
    ax.set_title("Length Distribution of Transcriptome ($input_fasta_file)")
    plt.tight_layout()
    plt.savefig('${input_fasta_file.baseName}_length_distribution.png')
    """ 
}