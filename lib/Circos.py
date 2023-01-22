import os
import shutil
from typing import List
from lib.Blast import Blast
from lib.Helpers import separate_array_values, get_center_value, get_fasta_length, change_str_to_int, write_to_file
from lib.Config_creator import create_config


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
        self.fasta_lengths = {}
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
        self.names_all = []
        self.connectors = {
            'minus': [],
            'plus': []
        }
        self.links = {
            'prot': [],
            'rrna': [],
            'trna': []
        }
        self.karyotype = []


    def make_links_for_seq_duplications(self, data):
        links = []
        for line_splitted in data:
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

            result = self.blast.make_blast(sequence, reference_path)

            for index, line in enumerate(result.stdout.splitlines()):
                # Skipping the first line containing the self-alignment
                if index > 0:
                    content.append(line)

            content_sorted = sorted(content, key=lambda item: int(item.split('\t')[7]), reverse=True)

            seq_len = int(content_sorted[0].split('\t')[1])
            seq_index = int(content_sorted[0].split('\t')[0].replace('h', '')) - 1
            self.fasta_lengths[seq_index] = seq_len

            contents.append(content_sorted)

        output = []
        for index, item in enumerate(contents):
            for line in item:
                line_splitted = line.split('\t')
                if int(line_splitted[-1]) >= self.config['min_duplication_length']:
                    line_splitted.insert(0, f'h{index + 1}')
                    output.append(line_splitted)

        output_sorted = sorted(output, key=lambda item: item[1])

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
            if prefix == '' and line_splitted[6] == strand and not line_splitted[0].startswith('trn') and not line_splitted[0].startswith('rrn'):
                input_arr = [f"h{str(index + 1)}"] + [line_splitted[4], line_splitted[5]]
                highlights.append('\t'.join(input_arr))
            elif prefix != '' and line_splitted[6] == strand and line_splitted[0].startswith(prefix):
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

    def get_names_all(self, collection: List[List[str]], seq_index: int) -> List[str]:
        output = []

        for index, line_splitted in enumerate(collection):
            start = int(line_splitted[4])
            end = int(line_splitted[5])

            input_arr = [f"h{str(seq_index + 1)}", str(start), str(end), line_splitted[0]]
            output.append('\t'.join(input_arr))

        return output

    def get_links_for_gene(self, collection: List[List[str]]):
        output = []
        indexes_to_change = [3, 4, 5, 6]
        collection = change_str_to_int(collection, indexes_to_change)

        gene_indexes = []
        for i in range(len(collection)):
            gene_indexes.append(i)

        for i in gene_indexes:
            for index, line_splitted in enumerate(collection):
                if index != i:
                    if collection[i][3] < collection[i][4] and collection[index][3] < collection[index][4]:
                        if (collection[i][3] <= collection[index][3] and collection[i][4] > collection[index][3]) \
                                or (collection[i][3] >= collection[index][3] and collection[i][3] < collection[index][4]):
                            if collection[i][3] < collection[index][3]:
                                start1 = collection[i][5] + (collection[index][3] - collection[i][3])
                                start2 = collection[index][5]
                            else:
                                start1 = collection[i][5]
                                start2 = collection[index][5] + (collection[i][3] - collection[index][3])

                            if collection[i][4] < collection[index][4]:
                                end1 = collection[i][6]
                                end2 = collection[index][6] - (collection[index][4] - collection[i][4])
                            else:
                                end1 = collection[i][6] - (collection[i][4] - collection[index][4])
                                end2 = collection[index][6]

                            output_arr = [collection[i][0], str(start1), str(end1), collection[index][0], str(start2), str(end2), 'twist=yes']
                            output.append('\t'.join(output_arr))
                    elif collection[i][3] > collection[i][4] and collection[index][3] > collection[index][4]:
                        if (collection[i][4] <= collection[index][4] and collection[i][3] > collection[index][4]) \
                                or (collection[i][4] >= collection[index][4] and collection[i][4] < collection[index][3]):
                            if collection[i][3] < collection[index][3]:
                                start1 = collection[i][5]
                                start2 = collection[index][5] + (collection[index][3] - collection[i][3])
                            else:
                                start1 = collection[i][5] + (collection[i][3] - collection[index][3])
                                start2 = collection[index][5]

                            if collection[i][4] < collection[index][4]:
                                end1 = collection[i][6] - (collection[index][4] - collection[i][4])
                                end2 = collection[index][6]
                            else:
                                end1 = collection[i][6]
                                end2 = collection[index][6] - (collection[i][4] - collection[index][4])

                            output_arr = [collection[i][0], str(start1), str(end1), collection[index][0], str(start2), str(end2), 'twist=yes']
                            output.append('\t'.join(output_arr))
                    elif collection[i][3] > collection[i][4] and collection[index][3] < collection[index][4]:
                        if (collection[i][4] < collection[index][4] and collection[i][3] > collection[index][3]) or (
                                collection[i][4] > collection[index][4] and collection[i][3] < collection[index][3]):
                            if collection[i][3] < collection[index][4]:
                                start1 = collection[i][5]
                                end2 = collection[index][6] - (collection[index][4] - collection[i][3])
                            else:
                                start1 = collection[i][5] + (collection[i][3] - collection[index][4])
                                end2 = collection[index][6]

                            if collection[i][4] < collection[index][3]:
                                end1 = collection[i][6] - (collection[index][3] - collection[i][4])
                                start2 = collection[index][5]
                            else:
                                end1 = collection[i][6]
                                start2 = collection[index][5] - (collection[i][4] - collection[index][3])

                            output_arr = [collection[i][0], str(start1), str(end1), collection[index][0], str(start2), str(end2), 'twist=no']
                            output.append('\t'.join(output_arr))
                    elif collection[i][3] < collection[i][4] and collection[index][3] > collection[index][4]:
                        if (collection[i][3] < collection[index][3] and collection[i][4] > collection[index][4]) or (
                                collection[i][3] > collection[index][3] and collection[i][4] < collection[index][4]):
                            if collection[i][3] < collection[index][4]:
                                start1 = collection[i][5] + (collection[index][4] - collection[i][3])
                                end2 = collection[index][6]
                            else:
                                start1 = collection[i][5]
                                end2 = collection[index][6] - (collection[i][3] - collection[index][4])

                            if collection[i][4] < collection[index][3]:
                                end1 = collection[i][6]
                                start2 = collection[index][5] + (collection[index][3] - collection[i][4])
                            else:
                                end1 = collection[i][6] - (collection[i][4] - collection[index][3])
                                start2 = collection[index][5]

                            output_arr = [collection[i][0], str(start1), str(end1), collection[index][0], str(start2), str(end2), 'twist=no']
                            output.append('\t'.join(output_arr))

        return output

    def get_links(self, collection: List[List[str]], gene_type: str):
        prefix = ''
        if gene_type == 'tRNA':
            prefix = 'trn'
        elif gene_type == 'rRNA':
            prefix = 'rrn'

        output = []
        gene_name = ""
        gene_collection = []

        for index, line_splitted in enumerate(collection):
            if prefix == '' and not line_splitted[1].startswith('trn') and not line_splitted[1].startswith('rrn'):
                if gene_name == "" or line_splitted[1] == gene_name:
                    gene_collection.append(line_splitted)
                else:
                    output.extend(self.get_links_for_gene(gene_collection))
                    gene_collection = [line_splitted]

                gene_name = line_splitted[1]
            elif prefix != '' and line_splitted[1].startswith(prefix):
                if gene_name == "" or line_splitted[1] == gene_name:
                    gene_collection.append(line_splitted)
                else:
                    output.extend(self.get_links_for_gene(gene_collection))
                    gene_collection = [line_splitted]
                gene_name = line_splitted[1]

        output.extend(self.get_links_for_gene(gene_collection))

        return output

    def create_circos_input_files(self):
        circos_input_files_dir_path = os.path.join(self.config['output_dir'], 'circos_input_files')
        if not os.path.exists(circos_input_files_dir_path):
            os.mkdir(circos_input_files_dir_path)

        write_to_file(self.highlights['minus_prot'],
                      os.path.join(circos_input_files_dir_path, 'minus_prot_names_highlights'))
        write_to_file(self.highlights['minus_rrna'],
                      os.path.join(circos_input_files_dir_path, 'minus_rRNA_names_highlights'))
        write_to_file(self.highlights['minus_trna'],
                      os.path.join(circos_input_files_dir_path, 'minus_tRNA_names_highlights'))
        write_to_file(self.highlights['plus_prot'],
                      os.path.join(circos_input_files_dir_path, 'plus_prot_names_highlights'))
        write_to_file(self.highlights['plus_rrna'],
                      os.path.join(circos_input_files_dir_path, 'plus_rRNA_names_highlights'))
        write_to_file(self.highlights['plus_trna'],
                      os.path.join(circos_input_files_dir_path, 'plus_tRNA_names_highlights'))

        write_to_file(self.names_all,
                      os.path.join(circos_input_files_dir_path, 'names'))
        write_to_file(self.names['minus'],
                      os.path.join(circos_input_files_dir_path, 'minus_names'))
        write_to_file(self.names['plus'],
                      os.path.join(circos_input_files_dir_path, 'plus_names'))

        write_to_file(self.connectors['minus'],
                      os.path.join(circos_input_files_dir_path, 'minus_names_connectors'))

        write_to_file(self.connectors['plus'],
                      os.path.join(circos_input_files_dir_path, 'plus_names_connectors'))

        write_to_file(self.links['prot'],
                      os.path.join(circos_input_files_dir_path, 'prot_links'))

        write_to_file(self.links['rrna'],
                      os.path.join(circos_input_files_dir_path, 'rRNA_links'))

        write_to_file(self.links['trna'],
                      os.path.join(circos_input_files_dir_path, 'tRNA_links'))

        write_to_file(self.karyotype,
                      os.path.join(circos_input_files_dir_path, 'karyotype.txt'))

        write_to_file(self.duplication_links,
                      os.path.join(circos_input_files_dir_path, 'duplication_links'))

    def create_karyotype(self):
        for index, file in enumerate(self.blast.fasta_files):
            line = f"chr\t-\th{index + 1}\t{self.blast.fasta_seq_names[index]}\t0\t{self.fasta_lengths[index]}\tchr{index + 1}"
            self.karyotype.append(line)

    def run(self):
        self.prepare_cds_reference()
        content_all_seq = []

        for index, sequence in enumerate(self.blast.fasta_files):
            result = self.blast.make_blast(sequence, self.cds_file_path)

            content = []
            for line in result.stdout.splitlines():
                content.append(line.split('\t'))
                content_all_seq.append([f"h{index + 1}"] + line.split('\t'))

            result_sorted_by_name = sorted(content, key=lambda x: x[0])

            modified_names = self.set_exon_names(result_sorted_by_name)

            self.highlights['minus_prot'].extend(self.get_highlights('minus', 'prot', modified_names, index))
            self.highlights['minus_rrna'].extend(self.get_highlights('minus', 'rRNA', modified_names, index))
            self.highlights['minus_trna'].extend(self.get_highlights('minus', 'tRNA', modified_names, index))
            self.highlights['plus_prot'].extend(self.get_highlights('plus', 'prot', modified_names, index))
            self.highlights['plus_rrna'].extend(self.get_highlights('plus', 'rRNA', modified_names, index))
            self.highlights['plus_trna'].extend(self.get_highlights('plus', 'tRNA', modified_names, index))

            self.names_all.extend(self.get_names_all(modified_names, index))
            self.names['minus'].extend(self.get_names(modified_names, 'minus', index, sequence))
            self.names['plus'].extend(self.get_names(modified_names, 'plus', index, sequence))

            self.connectors['minus'].extend(
                self.get_names(modified_names, 'minus', index, sequence, without_names=True))
            self.connectors['plus'].extend(
                self.get_names(modified_names, 'plus', index, sequence, without_names=True))

        content_all_seq_sorted_by_name = sorted(content_all_seq, key=lambda x: x[1])

        self.links['prot'] = self.get_links(content_all_seq_sorted_by_name, 'prot')
        self.links['rrna'] = self.get_links(content_all_seq_sorted_by_name, 'rRNA')
        self.links['trna'] = self.get_links(content_all_seq_sorted_by_name, 'tRNA')

        self.create_karyotype()

        self.create_circos_input_files()

        create_config(self.config)
