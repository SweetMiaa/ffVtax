import argparse
import bloom_filter_handler
from bloom_filter_handler import BloomFilterHandler
from sequence_matcher import SequenceMatcher
import os
import pandas as pd
import gca_to_taxid


def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description="Run viral sequence matching and calculate Jaccard and Qcov.")
    parser.add_argument('-i', '--input', required=True, help="Input query sequence file (FASTA format)")
    parser.add_argument('-d', '--database', default="db_bloom/bloom_filter.pkl", help="Path to the prebuilt bloom filter")
    parser.add_argument('-o', '--output', required=True, help="Output directory for CSV results")
    parser.add_argument('-jacc', type=float, default=0.95, help="Jaccard index threshold (default: 0.95)")
    parser.add_argument('-qcov', type=float, default=0.8, help="Query coverage threshold (default: 0.8)")
    parser.add_argument('-taxid', default="taxid.map", help= "the path to taxid.map file")
    args = parser.parse_args()

    bloom_handler = BloomFilterHandler()
    # 初始化 BloomFilterHandler 并加载数据库
    print("Loading reference database...")
    bloom_handler.load_bloom_filter(args.database)

    print("Database loaded successfully.")

    # 输出文件路径
    output_dir = args.output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    matching_result_file = os.path.join(output_dir, "matching_results.csv")

    # 初始化 SequenceMatcher 并进行序列匹配
    print(f"Matching sequences from {args.input} and calculating Jaccard and Qcov...")
    matcher = SequenceMatcher(bloom_handler)
    matcher.match_sequences(args.input, matching_result_file)
    print(f"Matching results saved to {matching_result_file}")

    gca_to_taxid.add_taxid_to_scored_output(matching_result_file, args.taxid)
    # 过滤匹配结果并保存到 CSV
    filter_and_sort_results(matching_result_file, os.path.join(output_dir, "filtered_results.csv"), args.jacc, args.qcov)
    print(f"Filtered results saved to {os.path.join(output_dir, 'filtered_results.csv')}")



def filter_and_sort_results(input_csv, output_csv, jacc_threshold=0.95, qcov_threshold=0.8):
    """
    读取匹配结果的 CSV 文件，根据 Jaccard 和 Qcov 阈值进行过滤，按优先级排序，并去重。
    最终将过滤和排序后的结果保存到新的 CSV 文件中。

    :param input_csv: 输入的包含匹配结果的 CSV 文件路径
    :param output_csv: 输出的过滤和排序后的 CSV 文件路径
    :param jacc_threshold: Jaccard 阈值，默认值为 0.95
    :param qcov_threshold: Query Coverage 阈值，默认值为 0.8
    """
    # 读取 CSV 文件到 pandas DataFrame
    data = pd.read_csv(input_csv)

    # 过滤出满足 Jaccard 和 Query Coverage 阈值的结果
    filtered_results = data[(data['Jacc'] >= jacc_threshold) & (data['Qcov'] >= qcov_threshold)]

    # 按照 Jaccard 和 Query Coverage 进行排序，优先考虑 Jaccard 其次是 Qcov
    filtered_results = filtered_results.sort_values(by=['Sequence Name', 'Jacc', 'Qcov'], ascending=[True, False, False])

    # 去除重复的 Sequence Name，只保留每个序列的最佳匹配
    filtered_result = filtered_results.drop_duplicates(subset=['Sequence Name'], keep='first')

    # 将过滤和排序后的结果保存到输出文件
    filtered_result.to_csv(output_csv, index=False)
    print(f"Filtered results saved to {output_csv}")


if __name__ == "__main__":
    main()
