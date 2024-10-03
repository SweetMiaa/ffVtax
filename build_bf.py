import argparse
from bloom_filter_handler import BloomFilterHandler
import concurrent.futures
import os


def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description="Load the database into bloom filter with multithreading support.")
    parser.add_argument('-i', '--input', required=True, help="Path to reference sequences")
    parser.add_argument('-o', '--output', default="db_bloom", help="Output directory for bloom filter")
    parser.add_argument('--kmer_size', type=int, default=21, help="K-mer size (default: 21)")
    parser.add_argument('--error_rate', type=float, default=0.001, help="Bloom filter error rate (default: 0.001)")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads for parallel processing (default: 4)")
    args = parser.parse_args()

    # 初始化 BloomFilterHandler 并加载数据库
    print("Loading reference database with multithreading...")
    bloom_handler = BloomFilterHandler(kmer_size=args.kmer_size, factor=4)
    bloom_handler.load_database_multithreaded(args.input, num_threads=args.threads)
    print("Database loaded successfully.")

    bloom_handler.save_bloom_filter(args.output)


if __name__ == "__main__":
    main()
