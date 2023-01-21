import os
import shutil
from typing import List
from lib.Blast import Blast
from lib.Helpers import separate_array_values, get_center_value, get_fasta_length


class Circos:
    def __init__(self, config):
        self.config = config
        self.blast = Blast(config)
        self.B_SUBJECT_START_INDEX = 3
        self.B_SUBJECT_END_INDEX = 4
        self.B_QUERY_START_INDEX = 5
        self.B_QUERY_END_INDEX = 6
        self.B_STRAND_ORIENTATION = 7
        self.cds_file_path = os.path.join(self.config['output_dir'], 'cds', os.path.basename(self.config['cds_reference']))
        self.duplication_links = self.get_duplication_links()
        self.highlights = {
            'minus_prot': [],
            'minus_rrna': [],
            'minus_trna': [],
            'plus_prot': [],
            'plus_rrna': [],
            'plus_trna': []
        }
        self.names = {
            'minus': [],
            'plus': []
        }
        self.connectors = {
            'minus': [],
            'plus': []
        }
        self.links = {
            'prot': [],
            'rrna': [],
            'trna': []
        }

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
        self.blast.change_fasta_seq_name()
        reference_path = os.path.join(self.config['output_dir'], 'blast', 'contigs_all.fasta')
        self.blast.merge_sequences(self.blast.fasta_files, reference_path)
        self.blast.make_db(reference_path)

        contents = []
        for sequence in self.blast.fasta_files:
            content = []
            output_path = os.path.join(self.config['output_dir'], 'blast',
                                       os.path.basename(sequence).replace('.fasta', '_out.txt'))
            result = self.blast.make_blast(sequence, reference_path)

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
                    line_splitted.insert(0, f'h{index + 1}')
                    output.append('\t'.join(line_splitted))

        output_sorted = sorted(output, key=lambda item: item.split('\t')[1])

        duplication_links = self.make_links_for_seq_duplications(output_sorted)

        return duplication_links

    def prepare_cds_reference(self):
        if not os.path.exists(os.path.dirname(self.cds_file_path)):
            os.mkdir(os.path.dirname(self.cds_file_path))

        shutil.copy(self.config['cds_reference'], self.cds_file_path)

        self.blast.make_db(self.cds_file_path)

    def set_exon_names(self, content):
        output = []

        for line in content:
            if line[6] == 'plus':
                start = line[2]
                end = line[3]
            else:
                start = line[3]
                end = line[2]

            output.append([f"{line[0]}_{start}-{end}"] + line[1:])

        return output

    def get_highlights(self, strand, gene_type, content, index):
        prefix = ''
        if gene_type == 'tRNA':
            prefix = 'trn'
        elif gene_type == 'rRNA':
            prefix = 'rrn'

        highlights = []
        for line_splitted in content:
            if line_splitted[6] == strand and line_splitted[0].startswith(prefix):
                input_arr = [f"h{str(index + 1)}"] + [line_splitted[4], line_splitted[5]]
                highlights.append('\t'.join(input_arr))

        return highlights

    def get_names(self, collection: List[List[str]], strand: str, seq_index: int, fasta_file_path: str, without_names=False) -> List[str]:
        output = []

        highlight_start_end_section = map(lambda x: x[4:6], collection)
        center_values = get_center_value(list(highlight_start_end_section))
        max_value = get_fasta_length(fasta_file_path)
        center_values_dispersed = separate_array_values(center_values, self.config["min_distance_between_labels"], max_value)

        for index, line_splitted in enumerate(collection):
            if not line_splitted[6] == strand:
                continue

            diff = center_values_dispersed[index] - center_values[index]
            start = int(line_splitted[4]) + diff
            end = int(line_splitted[5]) + diff

            if without_names:
                input_arr = [f"h{str(seq_index + 1)}", str(start), str(end)]
            else:
                input_arr = [f"h{str(seq_index + 1)}", str(start), str(end), line_splitted[0]]
            output.append('\t'.join(input_arr))

        return output

    def run(self):
        self.prepare_cds_reference()

        for index, sequence in enumerate(self.blast.fasta_files):

            result = self.blast.make_blast(sequence, self.cds_file_path)

            content = []
            for line in result.stdout.splitlines():
                content.append(line.split('\t'))

            result_sorted_by_name = sorted(content, key=lambda x: x[0])

            modified_names = self.set_exon_names(result_sorted_by_name)

            self.highlights['minus_prot'].extend(self.get_highlights('minus', 'prot', modified_names, index))
            self.highlights['minus_rrna'].extend(self.get_highlights('minus', 'rRNA', modified_names, index))
            self.highlights['minus_trna'].extend(self.get_highlights('minus', 'tRNA', modified_names, index))
            self.highlights['plus_prot'].extend(self.get_highlights('plus', 'prot', modified_names, index))
            self.highlights['plus_rrna'].extend(self.get_highlights('plus', 'rRNA', modified_names, index))
            self.highlights['plus_trna'].extend(self.get_highlights('plus', 'tRNA', modified_names, index))

            self.names['minus'].extend(self.get_names(modified_names, 'minus', index, sequence))
            self.names['plus'].extend(self.get_names(modified_names, 'plus', index, sequence))

            self.connectors['minus'].extend(
                self.get_names(modified_names, 'minus', index, sequence, without_names=True))
            self.connectors['plus'].extend(
                self.get_names(modified_names, 'plus', index, sequence, without_names=True))

        for i in self.connectors['plus']:
            print(i)


        # for key in self.highlights.keys():
        #     print(f'Key: {key}')
        #     for line in self.highlights[key]:
        #         print(line)

        exit()

