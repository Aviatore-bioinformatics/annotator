import os
import subprocess
from lib.Logger import log_information


class Blast:
    def __init__(self, config):
        self.config = config
        self.fasta_output_path = os.path.join(self.config['output_dir'], 'fasta')
        self.fasta_files = []

    def merge_sequences(self, fasta_sequences, output_path):
        log_information('Start merging FASTA sequences')
        merged = ''
        for sequence in fasta_sequences:
            with open(sequence, 'r') as file:
                merged += file.read()

        output_dir = os.path.dirname(output_path)
        if not os.path.exists(os.path.dirname(output_path)):
            os.makedirs(output_dir)

        if not os.path.exists(output_path):
            with open(output_path, 'w') as file:
                file.write(merged)
        else:
            log_information(f'The file "{output_path}" already exists')

    def make_db(self, reference):
        log_information(f"Creating a BLAST database from the file '{reference}'")

        result = subprocess.run(["makeblastdb", "-dbtype", "nucl", "-in", reference, "-parse_seqids"],
                                check=True, capture_output=True, text=True)

        log_file_path = os.path.join(os.path.dirname(reference), os.path.basename(reference).replace('.fasta', '_makeblastdb.log'))
        with open(log_file_path, 'w') as file:
            file.write(result.stdout)
            file.write(result.stderr)

    def make_blast(self, sequence, reference, output_dir=None):
        if output_dir is None:
            result = subprocess.run(["blastn", "-db", reference, "-query", sequence, "-outfmt",
                                     "6 sseqid slen sstart send qstart qend sstrand length"],
                                    capture_output=True, text=True, check=True)
            return result
        else:
            output_path = os.path.join(output_dir, os.path.basename(sequence).replace('.fasta', '_out.txt'))
            subprocess.run(["blastn", "-db", reference, "-query", sequence, "-outfmt", "6 sseqid slen sstart send qstart qend sstrand length", "-out", output_path],
                           check=True)

    def change_fasta_seq_name(self):
        if not os.path.exists(self.fasta_output_path):
            os.mkdir(self.fasta_output_path)

        for index, sequence in enumerate(self.config['sequences']):
            self.fasta_files.append(os.path.join(self.fasta_output_path, os.path.basename(sequence)))
            with open(os.path.join(self.fasta_output_path, os.path.basename(sequence)), 'w') as out:
                with open(sequence, 'r') as file:
                    for line in file:
                        if line.startswith('>'):
                            out.write(f">h{index}\n")
                        else:
                            out.write(line)
