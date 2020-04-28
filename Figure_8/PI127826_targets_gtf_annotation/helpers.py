# add blast header to blast result files
def add_blast_header_to_file(blast_file_without_header,blast_file_with_header):
    """takes a blast result file (outformat 6 not customized) and add a header
    """
    blast_header = ["qseqid",
                "subject_id",
                "pct_identity",
                "aln_length",
                "n_of_mismatches",
                "gap_openings",
                "q_start",
                "q_end",
                "s_start",
                "s_end",
                "e_value",
                "bit_score"]
    df = pd.read_csv(blast_file_without_header,sep="\t",header=None)
    df.columns = blast_header
    df.to_csv(blast_file_with_header,sep="\t",header=True,index=False)