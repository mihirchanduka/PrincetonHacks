
import vcf as vcf
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def generate_mutated_genome(base_dir):
    """
    Generates a mutated genome by applying variants from a VCF file to a reference genome.
    """
    # Define file paths
    reference_file = f"{base_dir}/genome/reference.fasta"
    vcf_file = f"{base_dir}/genome/variants.vcf"
    output_file = f"{base_dir}/genome/mutated_genome.fasta"

    # Read the reference genome
    reference = SeqIO.read(reference_file, "fasta")
    mutable_seq = MutableSeq(reference.seq)

    # Read the VCF file and apply variants
    with open(vcf_file, 'r') as vcf_reader:
        reader = vcf.Reader(vcf_reader)
        for record in reader:
            # VCF is 1-based, Biopython is 0-based
            pos = record.POS - 1
            ref = record.REF
            alt = record.ALT[0]

            # Check if the reference base matches
            if mutable_seq[pos] == ref:
                mutable_seq[pos] = str(alt)
            else:
                print(f"Warning: Reference base mismatch at position {record.POS}. Expected {ref}, found {mutable_seq[pos]}.")


    # Create a new SeqRecord with the mutated sequence
    mutated_record = SeqRecord(
        Seq(mutable_seq),
        id=f"{reference.id}_mutated",
        description="Mutated sequence"
    )

    # Write the mutated sequence to a new FASTA file
    SeqIO.write(mutated_record, output_file, "fasta")
    print(f"Mutated genome saved to {output_file}")

if __name__ == '__main__':
    generate_mutated_genome('.')
