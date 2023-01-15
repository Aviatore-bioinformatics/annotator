import os
from lib.Blast import Blast


class Circos:
    def __init__(self, config):
        self.config = config
        self.B_SUBJECT_START_INDEX = 3
        self.B_SUBJECT_END_INDEX = 4
        self.B_QUERY_START_INDEX = 5
        self.B_QUERY_END_INDEX = 6
        self.B_STRAND_ORIENTATION = 7

    def make_links_for_seq_duplications(self, data):
        links = []
        for line in data:
            line_splitted = line.split('\t')

            chrom_prev = line_splitted[1]
            chrom = line_splitted[0]

            sstart = int(line_splitted[self.B_SUBJECT_START_INDEX])
            send = int(line_splitted[self.B_SUBJECT_END_INDEX])
            qstart = int(line_splitted[self.B_QUERY_START_INDEX])
            qend = int(line_splitted[self.B_QUERY_END_INDEX])
            strand = line_splitted[self.B_STRAND_ORIENTATION]

            if sstart < send:
                start_1 = sstart
                end_1 = send
            else:
                start_1 = send
                end_1 = sstart

            if qstart < qend:
                start_2 = qstart
                end_2 = qend
            else:
                start_2 = qend
                end_2 = qstart

            links.append(f"{chrom_prev}\t{start_1}\t{end_1}\t{chrom}\t{start_2}\t{end_2}\ttwist={'yes' if strand == 'plus' else 'no'}")

        return links

    def get_duplication_links(self):
        blast = Blast(self.config)

        blast.change_fasta_seq_name()
        reference_path = os.path.join(self.config['output_dir'], 'blast', 'contigs_all.fasta')
        blast.merge_sequences(blast.fasta_files, reference_path)
        blast.make_db(reference_path)

        contents = []
        for sequence in blast.fasta_files:
            content = []
            output_path = os.path.join(self.config['output_dir'], 'blast',
                                       os.path.basename(sequence).replace('.fasta', '_out.txt'))
            result = blast.make_blast(sequence, reference_path)

            # with open(output_path, 'w') as file:
            #     for index, line in enumerate(result.stdout.splitlines()):
            #         # Removing the first line containing the self-alignment
            #         if index > 0:
            #             file.write(f'{line}\n')

            for index, line in enumerate(result.stdout.splitlines()):
                # Skipping the first line containing the self-alignment
                if index > 0:
                    content.append(line)

            content_sorted = sorted(content, key=lambda item: int(item.split('\t')[7]), reverse=True)

            contents.append(content_sorted)

        output = []
        for index, item in enumerate(contents):
            for line in item:
                line_splitted = line.split('\t')
                if int(line_splitted[-1]) >= self.config['min_duplication_length']:
                    line_splitted.insert(0, f'h{index}')
                    output.append('\t'.join(line_splitted))

        output_sorted = sorted(output, key=lambda item: item.split('\t')[1])

        duplication_links = self.make_links_for_seq_duplications(output_sorted)

        return duplication_links
