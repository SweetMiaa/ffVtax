import csv
from collections import defaultdict
import gzip

class SequenceMatcher:
    def __init__(self, bloom_filter_handler):
        self.bloom_filter_handler = bloom_filter_handler

    def match_sequences(self, query_file, output_file):
        """
        从输入文件中逐条读取序列，匹配参考序列，计算 Jaccard 和 Qcov，然后将结果逐条写入 CSV。
        :param query_file: 输入的 query 序列文件路径 (FASTA 格式)。
        :param output_file: 输出的 CSV 文件路径。
        """
        with open(output_file, mode='w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            # 写入表头
            writer.writerow(['Sequence Name', 'Total Input k-mers', 'Reference', 'Matched k-mers', 'Reference Total k-mers', 'GCA Name', 'Jacc', 'Qcov'])

            # 判断是文件路径还是文件对象
            if isinstance(query_file, str):
                open_func = gzip.open if query_file.endswith('.gz') else open
                with open_func(query_file, 'rt') as file:
                    for sequence_name, sequence_data in self._read_fasta_file(file):
                        self._match_sequence(sequence_name, sequence_data, writer)
            else:
                for sequence_name, sequence_data in self._read_fasta_file(query_file):
                    self._match_sequence(sequence_name, sequence_data, writer)

    def _read_fasta_file(self, file):
        sequence_name = None
        sequence_data = []

        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_name is not None:
                    yield sequence_name, ''.join(sequence_data)
                sequence_name = line[1:]
                sequence_data = []
            else:
                sequence_data.append(line)

        if sequence_name is not None:
            yield sequence_name, ''.join(sequence_data)

    def _match_sequence(self, sequence_name, sequence_data, writer):
        """
        计算输入序列的 k-mer 数量，匹配到的参考序列及其匹配的 k-mer 数量，同时计算 Jacc 和 Qcov。
        :param sequence_name: 输入的序列名称。
        :param sequence_data: 输入的序列数据。
        :param writer: CSV writer，用于将匹配结果写入 CSV 文件。
        """
        input_kmers = set(self._extract_kmers(sequence_data))  # 将输入序列的 k-mers 转换为集合
        total_input_kmers = len(input_kmers)  # 计算输入序列的总 k-mers 数量
        reference_kmer_matches = defaultdict(set)  # 存储每个参考序列匹配到的 k-mers

        # 匹配 k-mers 到参考序列
        for kmer in input_kmers:
            if kmer in self.bloom_filter_handler.bloom:
                for ref in self.bloom_filter_handler.kmer_to_reference[kmer]:
                    reference_kmer_matches[ref].add(kmer)

        # 遍历每个匹配的参考序列，计算 Jacc 和 Qcov
        for ref, matched_kmers in reference_kmer_matches.items():
            reference_total_kmers = self.bloom_filter_handler.reference_kmer_count.get(ref, 0)  # 获取参考序列的总 k-mers 数量
            gca_name = self.bloom_filter_handler.reference_gca_map.get(ref, 'Unknown')  # 获取 GCA 名称

            # 将参考序列的 k-mers 转换为集合（假设有对应的存储或提取方法）
            reference_kmers = self.bloom_filter_handler.get_reference_kmers(ref)  # 获取参考序列的 k-mers
            reference_kmers_set = set(reference_kmers)

            # 计算 Jaccard 和 Query Coverage
            if len(input_kmers) > 0 and len(reference_kmers_set) > 0:
                jacc = len(input_kmers & reference_kmers_set) / len(input_kmers | reference_kmers_set)  # Jaccard 指数
                qcov = len(input_kmers & reference_kmers_set) / len(reference_kmers_set)  # Query Coverage
            else:
                jacc = 0
                qcov = 0

            # 打印调试信息
            print(f"Sequence Name: {sequence_name}, Total Input k-mers: {total_input_kmers}, Reference: {ref}, "
                  f"Matched k-mers: {len(matched_kmers)}, Reference Total k-mers: {reference_total_kmers}, Jacc: {jacc}, Qcov: {qcov}")

            # 将结果写入 CSV
            writer.writerow([sequence_name, total_input_kmers, ref, len(matched_kmers), reference_total_kmers, gca_name, f"{jacc:.4f}", f"{qcov:.4f}"])

        # 清理内存
        del input_kmers
        del reference_kmer_matches

    def _extract_kmers(self, sequence):
        """
        从序列中提取 k-mers。
        """
        return [sequence[i:i + self.bloom_filter_handler.kmer_size] for i in range(len(sequence) - self.bloom_filter_handler.kmer_size + 1)]
