"""
在该模块中建立序列比对需要的索引。

@author: zhang heng
"""

import time


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


def load_reads(filename):
    '''
    对需要比对的read进行处理，包括将文件中的read抽取出来，并将数据全部转换为大写
    :param filename: 存储read的文件
    :return: 提取出来的read列表
    '''
    pass


def sort_starts(starts):
    sorted_start = []
    for v in starts.values():
        sorted_start += v
    return sorted(sorted_start)


def get_unipath_starts(ref, k_mer):
    '''
    这是建立deBGA的需要的索引的时候需要用到的一个辅助数据结构，这里会顺带统计k-mer
    :param filename: 文件
    :return: 统计的k-mer以及各unipath的起始点
    '''
    UT_kmer = dict()  # 索引中的哈希表
    UT_starts = dict()  # 辅助数据结构，用于存储所有unipath的起始位置
    UT_starts[1] = [1]  # 初始化
    flag = True  # 标志变量，True表明当前是一个新的unipath，False表示当前是重复一个unipath
    for i in range(len(ref) - k_mer + 1):
        if UT_kmer.get(ref[i:i + k_mer]) is None and flag:
            UT_kmer[ref[i:i + k_mer]] = i + 1
        elif UT_kmer.get(ref[i:i + k_mer]) is None and flag is False:
            UT_kmer[ref[i:i + k_mer]] = i + 1
            flag = True
            UT_starts[i + 1] = [i + 1]
            next_kmer = UT_kmer.get(ref[i - 1:i + k_mer - 1])
            if UT_kmer.get(ref[next_kmer:next_kmer + k_mer]) == next_kmer + 1 and UT_starts.get(next_kmer + 1) is None:
                UT_starts[next_kmer + 1] = [next_kmer + 1]
        elif UT_kmer.get(ref[i:i + k_mer]) is not None and flag:
            flag = False
            if UT_starts.get(UT_kmer.get(ref[i:i + k_mer])) is not None:
                UT_starts.get(UT_kmer.get(ref[i:i + k_mer])).append(i + 1)
            else:
                UT_starts[UT_kmer.get(ref[i:i + k_mer])] = [UT_kmer.get(ref[i:i + k_mer]), i + 1]
        else:
            if UT_starts.get(UT_kmer.get(ref[i:i + k_mer])) is not None:
                UT_starts.get(UT_kmer.get(ref[i:i + k_mer])).append(i + 1)
    return UT_kmer, UT_starts


def get_index(filename, k_mer):
    ref_sequence = load_ref(filename)
    UT_kmer, UT_starts = get_unipath_starts(ref_sequence, k_mer)
    first_starts = sorted(UT_starts.keys())
    sorted_starts = sort_starts(UT_starts)
    sorted_starts.append(len(ref_sequence) - k_mer + 2)
    UT_seq = []
    UT_pos = []
    for i in range(0, len(first_starts), 1):
        start = first_starts[i]
        start_index = sorted_starts.index(start)
        length = sorted_starts[start_index + 1] - start + k_mer - 1
        UT_seq.append([length, ref_sequence[start - 1:start + length - 1]])
        UT_pos.append(UT_starts.get(start))

    uids = [0] * len(ref_sequence)
    print(len(uids))
    print(len(UT_kmer))
    for i in range(len(UT_pos)):
        for j in range(UT_pos[i][0], UT_pos[i][0] + UT_seq[i][0] - k_mer + 1, 1):
            uids[j] = i + 1
    # i = 0

    for mer in UT_kmer.keys():
        pos = UT_kmer.get(mer)
        uid = uids[pos]
        # while UT_starts.get(pos) is None:
        #     pos -= 1
        offset = pos - UT_pos[uid - 1][0] + 1
        # uid = first_starts.index(pos) + 1
        UT_kmer[mer] = [uid, offset]
        # if i % 1000 == 0:
        #     print(i / 1000)
        # i += 1

    print(UT_kmer)
    print(UT_seq)
    print(UT_pos)


if __name__ == '__main__':
    start = time.time()
    # get_index('../data/NC_008253.fna', 17)
    get_index('../data/test.txt', 4)
    end = time.time()
    print(end - start)
