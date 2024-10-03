import os
import gzip
from collections import defaultdict
from pybloom_live import BloomFilter
import re
import pickle
import concurrent.futures

class BloomFilterHandler:
    def __init__(self, kmer_size=21, factor=4):
        """
        初始化 Bloom filter 处理器。
        :param kmer_size: k-mer 的长度。
        :param factor: 用于计算所需容量的因子（默认 4）。
        """
        self.kmer_size = kmer_size
        self.factor = factor
        self.error_rate = 0.001
        self.bloom = None
        self.kmer_to_reference = defaultdict(set)
        self.reference_kmer_count = defaultdict(int)
        self.reference_gca_map = {}
        self.reference_sequences = {}

    def estimate_kmers(self, file_path):
        """
        估算数据库文件中的 k-mer 数量。
        :param file_path: 数据库文件的路径。
        :return: 估计的 k-mer 数量。
        """
        kmer_count = 0
        sequence = ""
        open_func = gzip.open if file_path.endswith('.gz') else open

        try:
            with open_func(file_path, 'rt') as db_file:
                print(f"Reading file: {file_path}")  # 输出文件路径
                for line in db_file:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith('>'):
                        # 处理并计算前一个累积的序列
                        kmer_count += self._count_kmers_in_sequence(sequence)
                        sequence = ""
                    else:
                        sequence += line

                # 处理最后一个累积的序列
                kmer_count += self._count_kmers_in_sequence(sequence)
                print(f"Estimated k-mers for {file_path}: {kmer_count}")
        except Exception as e:
            print(f"Error reading file {file_path}: {e}")

        return kmer_count

    def _count_kmers_in_sequence(self, sequence):
        """
        辅助函数：计算给定序列中的 k-mers 数量。
        """
        if len(sequence) >= self.kmer_size:
            return len(sequence) - self.kmer_size + 1
        return 0

    def load_database_multithreaded(self, db_path, num_threads=4):
        """
        加载数据库文件并构建 Bloom filter，支持多线程。
        :param db_path: 数据库文件夹的路径。
        :param num_threads: 并行处理的线程数。
        """
        estimated_kmers = self._estimate_total_kmers_in_db(db_path)

        if estimated_kmers > 0:
            capacity = estimated_kmers * self.factor
            self._initialize_bloom_filter(capacity)
        else:
            raise ValueError("Estimated k-mers count is zero. Check the database files.")

        # 使用线程池并行处理文件
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = []
            for root, _, files in os.walk(db_path):
                for file_name in files:
                    file_path = os.path.join(root, file_name)
                    if file_path.endswith(('.fna', '.fna.gz', '.fasta', '.fasta.gz')):
                        futures.append(executor.submit(self._process_file_kmers, file_path))

            # 等待所有线程完成
            concurrent.futures.wait(futures)

        print(f"Total loaded k-mers across all threads.")

    def _estimate_total_kmers_in_db(self, db_path):
        """
        辅助函数：估算整个数据库中的 k-mers。
        """
        total_kmers = 0
        for root, _, files in os.walk(db_path):
            for file_name in files:
                file_path = os.path.join(root, file_name)
                if file_path.endswith('.fna') or file_path.endswith('.fna.gz'):
                    total_kmers += self.estimate_kmers(file_path)
        return total_kmers

    def _initialize_bloom_filter(self, capacity):
        """
        初始化 Bloom filter。
        :param capacity: Bloom filter 的容量。
        """
        self.bloom = BloomFilter(capacity=capacity, error_rate=self.error_rate)
        print(f"Initializing Bloom filter with capacity {capacity}...")

    def _load_kmers_into_bloom(self, db_path):
        """
        辅助函数：将数据库中的 k-mers 加载到 Bloom filter 中。
        """
        actual_kmer_count = 0

        for root, _, files in os.walk(db_path):
            for file_name in files:
                file_path = os.path.join(root, file_name)
                if (file_path.endswith('.fna') or file_path.endswith('.fna.gz') or
                        file_path.endswith('.fasta') or file_path.endswith('.fasta.gz')):
                    actual_kmer_count += self._process_file_kmers(file_path)

        print(f"Total loaded k-mers: {actual_kmer_count}")

    def _extract_gca_from_filename(self, file_path):
        """
        从文件名中提取 GCA 编号。
        :param file_path: 文件的完整路径。
        :return: GCA 编号，未找到则返回 'Unknown'。
        """
        gca_match = re.search(r'(GCA_\d+\.\d+)', file_path)
        if gca_match:
            return gca_match.group(1)
        return "Unknown"

    def _process_file_kmers(self, file_path):
        """
        辅助函数：从文件中提取 k-mers，并将其添加到 Bloom filter，同时保留 GCA 编号。
        """
        actual_kmer_count = 0
        current_reference = None
        gca_name = self._extract_gca_from_filename(file_path)  # 提取 GCA 编号
        open_func = gzip.open if file_path.endswith('.gz') else open

        sequence_data = []
        with open_func(file_path, 'rt') as db_file:
            for line in db_file:
                line = line.strip()
                if line.startswith('>'):
                # 如果有已处理的参考序列，提取其 k-mers 并加载到 Bloom filter
                    if current_reference and sequence_data:
                        actual_kmer_count += self._add_kmers_to_bloom(''.join(sequence_data), current_reference)
                        sequence_data = []

                    current_reference = self._extract_reference_name(line)  # 提取参考序列名称
                    self.reference_gca_map[current_reference] = gca_name  # 将 GCA 名字与参考序列一起存储
                else:
                    sequence_data.append(line)  # 累积序列数据

        # 处理最后一个参考序列
            if current_reference and sequence_data:
                actual_kmer_count += self._add_kmers_to_bloom(''.join(sequence_data), current_reference)

        return actual_kmer_count

    def _add_kmers_to_bloom(self, sequence, reference_name):
        """
        辅助函数：提取序列中的 k-mers，并将其添加到 Bloom filter 中。
        :param sequence: 输入的序列数据。
        :param reference_name: 参考序列的名称。
        :return: 提取的 k-mers 数量。
        """
        kmers = self.extract_kmers(sequence)

        print(f"Reference Name: {reference_name}, Extracted k-mers: {len(kmers)}")

        for kmer in kmers:
            self.bloom.add(kmer)
            self.kmer_to_reference[kmer].add(reference_name)

        self.reference_kmer_count[reference_name] += len(kmers)  # 更新参考序列的 k-mers 计数
        return len(kmers)  # 返回提取到的 k-mers 数量


    def extract_kmers(self, sequence):
        """
        从序列中提取 k-mers。
        """
        return [sequence[i:i + self.kmer_size] for i in range(len(sequence) - self.kmer_size + 1)]

    def _extract_reference_name(self, header_line):
        """
        辅助函数：从序列头提取参考名称。
        """
        if header_line.startswith('>'):
            return header_line.split()[0][1:]
        else:
            print(f"Warning: Header line format is incorrect: {header_line}")
        return "Unknown"

    def get_reference_sequence(self, reference_name):
        """
        根据参考序列名获取对应的序列数据。
        :param reference_name: 参考序列名称。
        :return: 对应的参考序列数据。
        """
        return self.reference_sequences.get(reference_name, '')

    def get_reference_kmers(self, reference_name):
        """
        根据参考序列名获取参考序列中的所有 k-mers。
        :param reference_name: 参考序列的名称。
        :return: 参考序列的所有 k-mers 列表。
        """
        # 在添加 k-mers 到 Bloom filter 时，也保存每个参考序列的所有 k-mers
        reference_kmers = [kmer for kmer, refs in self.kmer_to_reference.items() if reference_name in refs]
        return reference_kmers

    def save_bloom_filter(self, directory):
        """
        将 Bloom filter 保存到文件中。
        :param directory: 保存的文件路径。
        """
        if not os.path.exists(directory):
            os.makedirs(directory)

        # 构造保存的文件路径
        file_path = os.path.join(directory, 'bloom_filter.pkl')
        data_to_save = {
            'bloom': self.bloom,
            'kmer_to_reference': self.kmer_to_reference,
            'reference_kmer_count': self.reference_kmer_count,
            'reference_gca_map': self.reference_gca_map,
            'reference_sequences': self.reference_sequences
        }

        # 保存 Bloom filter 到文件
        with open(file_path, 'wb') as bloom_file:
            pickle.dump(data_to_save, bloom_file)
        print(f"Bloom filter saved to {file_path}")

    def load_bloom_filter(self, file_path):
        """
        从文件加载 Bloom filter。
        :param file_path: 保存的 Bloom filter 文件路径。
        """
        if os.path.exists(file_path):
            with open(file_path, 'rb') as bloom_file:
                data_loaded = pickle.load(bloom_file)

        # 恢复保存的属性
            self.bloom = data_loaded['bloom']
            self.kmer_to_reference = data_loaded['kmer_to_reference']
            self.reference_kmer_count = data_loaded['reference_kmer_count']
            self.reference_gca_map = data_loaded['reference_gca_map']
            self.reference_sequences = data_loaded['reference_sequences']
            print(f"Bloom filter loaded from {file_path}")
        else:
            print(f"Error: Bloom filter file {file_path} does not exist.")
            raise FileNotFoundError(f"Bloom filter file {file_path} not found.")
