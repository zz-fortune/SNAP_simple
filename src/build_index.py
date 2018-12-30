"""
在该模块中建立序列比对需要的索引。

@author: zhang heng
"""

import time
import math


def load_file(filename):
    '''
    从文件中读取数据
    :param filename: 文件名字
    :return: 读取出来的原始数据
    '''
    with open(filename, 'r') as fp:
        data = fp.read()
    return data


def load_ref(filename):
    '''
    对参考序列进行一些处理，去掉描述信息并将数据全部转换为大写
    :param filename: 存储参考序列的文件
    :return: 全部大写的参考序列
    '''
    data = load_file(filename)
    return data[data.index('\n') + 1:].replace('\n', '').upper()


def load_reads(file1, file2, max_num=math.inf):
    '''
    对需要比对的read进行处理，包括将文件中的read抽取出来，并将数据全部转换为大写。

    这里因为是pair-end的read，所以同时将两个文件中的reads分别取出来
    :param file1: 存储read的文件
    :param file2: 存储read的文件
    :param max_num: 读出的reads最大对数，默认是全部读出
    :return: 提取出来的read列表
    '''
    with open(file1, 'r') as fp:
        data1 = fp.read()
    with open(file2, 'r') as fp:
        data2 = fp.read()
    data1 = data1.split('\n')
    data2 = data2.split('\n')
    if max_num == math.inf or max_num > len(data1) / 4:
        num = len(data1) / 4
    elif max_num <= 0:
        num = 0
    else:
        num = max_num
    reads1 = []
    reads2 = []
    for i in range(int(num)):
        reads1.append(data1[i * 4 + 1].upper())
        reads2.append(data2[i * 4 + 1].upper())
    return reads1, reads2


def get_index(filename, k_mer):
    """
    建立索引。遍历一遍参考基因组
    :param filename: 存储参考基因组的文件
    :param k_mer: k-mer的长度
    :return: k-mer 的哈希索引以及原始的参考序列
    """
    ref_sequence = load_ref(filename)  # 完整的参考基因组序列
    kmer_index = dict()  # k-mer的索引表
    for i in range(len(ref_sequence) - k_mer + 1):
        if kmer_index.get(ref_sequence[i:i + k_mer]) is None:
            kmer_index[ref_sequence[i:i + k_mer]] = [i + 1]
        else:
            kmer_index[ref_sequence[i:i + k_mer]].append(i + 1)
    return kmer_index, ref_sequence


if __name__ == '__main__':
    start = time.time()
    get_index('../data/NC_008253.fna', 17)
    # get_index('../data/test.txt', 4)
    end = time.time()
    print(end - start)
