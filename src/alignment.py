"""
将read和参考基因组进行比对

@author: zhang heng
"""
import copy
import numpy
from src.build_index import *


def unscore_most_hitting(seeds_hitting, scores):
    """
    找到还没有评分（也就是就算编辑距离）的有最多hitting的位置
    :param seeds_hitting: 记录每个位置的hitting次数的表
    :param scores: 记录每个位置得分（编辑距离）的表
    :return: 还没有评分（也就是就算编辑距离）的有最多hitting的位置。如果位置均已评分，返回0
    """
    max_hitting, max_pos = -1, 0  # 分别表示hitting数和对应的位置
    for k in seeds_hitting:
        if scores.get(k) is None and seeds_hitting[k] > max_hitting:
            max_hitting = seeds_hitting[k]
            max_pos = k
    return max_pos


def edit_distance(read, candidate, d_limit):
    """
    计算read和参考基因组上的预测位置的编辑距离。
    具体来说就是计算read与参考基因组的一个子串的编辑距离，这个子串起始位置是利用k-mer预测的起始位置，其长度
    比read稍长一点，这是考虑到允许的最大的indel的长度。最终在编辑距离中，取从与read相同长度的子串到允许的最
    长子串和read的编辑距离的最小值。
    :param read: 待比对的片段
    :param candidate: 参考基因组长的一个子串，预测这个子串和read匹配
    :param d_limit: 允许的最大编辑距离，大于这个距离就视为距离无穷大
    :return: 计算出的编辑距离
    """
    m, n = len(candidate), len(read)
    d = numpy.zeros((m + 1, n + 1))
    # 初始化距离矩阵
    for i in range(m):
        d[i, 0] = i
    for i in range(n):
        d[0, i] = i

    # 使用动态规划算法计算编辑距离
    for i in range(1, m + 1, 1):
        for j in range(1, n + 1, 1):

            # 分为两种情况，一种是接线来待比较的两个字符相等
            if candidate[i - 1] == read[j - 1]:
                if d[i - 1, j - 1] <= d[i, j - 1] + 1 and d[i - 1, j - 1] <= d[i - 1, j] + 1:
                    d[i, j] = d[i - 1, j - 1]
                elif d[i - 1, j] < d[i - 1, j]:
                    d[i, j] = d[i - 1, j] + 1
                else:
                    d[i, j] = d[i, j - 1] + 1

            # 另一种是待比较的两个字符不等
            else:
                if d[i - 1, j - 1] <= d[i, j - 1] and d[i - 1, j - 1] <= d[i - 1, j]:
                    d[i, j] = d[i - 1, j - 1] + 1
                elif d[i - 1, j] < d[i - 1, j]:
                    d[i, j] = d[i - 1, j] + 1
                else:
                    d[i, j] = d[i, j - 1] + 1

    # 判断考虑到indel而添加的字符是否增大了编辑距离，找到将read完全比对的最小的编辑距离
    d_min = d[n, n]
    c = len(candidate) - len(read)
    for i in range(n - c, m + 1, 1):
        if d[i, n] < d_min:
            d_min = d[i, n]
    if d_min > d_limit:
        return math.inf
    else:
        return d_min


def update_distance(d_best, d_second, d_cur, pos_best, pos_second, pos):
    """
    根据最新计算出的一个编辑距离，更新最小的两个编辑距离以及对应的位置
    :param d_best: 原有的最小的编辑距离
    :param d_second: 原有的次小的编辑距离
    :param d_cur: 最新一个计算出的编辑距离
    :param pos_best: 原有的最小的编辑距离对应的位置
    :param pos: 最新一个计算出的编辑距离对应的位置
    :return: 两个编辑距离以及对应的位置
    """
    if d_cur < d_best:
        return copy.copy(d_cur), copy.copy(d_best), copy.copy(pos), copy.copy(pos_best)
    elif d_cur < d_second:
        return copy.copy(d_best), copy.copy(d_cur), copy.copy(pos_best), copy.copy(pos)
    else:
        return d_best, d_second, pos_best, pos_second


def align(reference, kmer_index, read, k_mer, h_max, d_max, c_threshold):
    """
    比对算法的主要部分。将read与参考基因组进行比对，找到read在参考基因组上的可能位置。本算法中只考虑编辑距离
    最小的两个位置。
    :param reference: 参考基因组序列
    :param kmer_index: k-mer的索引表
    :param read: 待比对的片段
    :param k_mer: k-mer的长度
    :param h_max: 每个k-mer允许的最大的匹配数，有些k-mer如“AAAAAAAA”这样的重复序列可能匹配到上千的位置，忽略这些位置
    :param d_max: 允许的read与找到的基因序列的子串之间最大的编辑距离
    :param c_threshold: 置信度，小于这个阈值的位置可以认为已经足够好。此外，第二好的位置如果大于第一的位置超过这个阈值，认为第二好的位置无效
    :return: ead在参考基因组上的可能位置。主要考虑两个，可能没有或者只有一个
    """
    non_overlapping_num = 0  # 用以记录已经考虑过的没有交叠的k-mer的个数
    d_best, d_second = math.inf, math.inf  # 记录最小的两个编辑距离
    pos_best, pos_second = 0, 0  # 记录最小的两个编辑距离对应的位置
    seeds_hitting = dict()  # 记录每个位置的hitting数
    scores = dict()  # 记录每个位置的编辑距离
    for i in range(len(read) - k_mer + 1):  # 依次考虑每个k-mer
        if i % k_mer == 0:  # 这里是依次考虑每个k-mer，所以每k_mer个k-mer会出现一个不交叠的k-mer
            non_overlapping_num += 1
        seed = read[i:i + k_mer]
        locations = (kmer_index.get(seed) or [])  # 在k-mer索引表中查找
        if len(locations) <= h_max:

            # 记录每个当前的k-mer都hit到了哪些位置
            for l in locations:
                pos = l - i
                if pos > 0:
                    seeds_hitting[pos] = (seeds_hitting.get(pos) or 0) + 1
        pos = unscore_most_hitting(seeds_hitting, scores)
        if pos == 0:  # 说明没有未计算编辑距离的位置，直接开始下一轮循环
            continue

        # 更新d_limit的值，更小的d_limit使得计算编辑距离时更快
        if d_best > d_max:
            d_limit = d_max + c_threshold - 1
        elif d_second >= d_best + c_threshold:
            d_limit = d_best + c_threshold - 1
        else:
            d_limit = d_best - 1

        # 计算编辑距离，并更新最小的编辑距离
        if pos - 1 < 0 or pos + len(read) > len(reference):
            continue
        d_cur = edit_distance(read, reference[pos - 1:pos + len(read) + c_threshold], d_limit)
        d_best, d_second, pos_best, pos_second = update_distance(d_best, d_second, d_cur, pos_best, pos_second, pos)
        scores[pos] = d_cur

        # 出现在置信度内的位置，提前结束
        if d_best < c_threshold and d_second < d_best + c_threshold:
            return [pos_best, pos_second]

        # 如果不交叠的k-mer的数量为t，则这个read与参考基因组的距离至少为t，这里也可以提前结束循环
        elif non_overlapping_num >= d_best + c_threshold:
            for k in seeds_hitting.keys():
                if scores.get(k) is None:
                    scores[k] = edit_distance(read, reference[k - 1:k + len(read) + c_threshold], d_limit)
            break

    # 如果次小的距离与最小的距离差距比较大，认为只有一个是有效
    if d_best <= d_max and d_second >= d_best + c_threshold:
        return [pos_best]

    # 如果没有一个位置有命中，说明read和参考基因组存在“AAAAAA”这种串，这个可随意返回一些位置
    elif len(seeds_hitting) == 0:
        for i in range(len(read) - k_mer):
            if kmer_index.get(read[i:i + k_mer]):
                return kmer_index[read[:k_mer]]
        return []

    # 如果最小的距离在允许范围内，就返回那些在允许范围内的位置
    elif d_best <= d_max:
        result = []
        for k in scores.keys():
            if scores.get(k) <= d_max:
                result.append(k)
        return result
    # 否则，说明read在参考基因组上找不到
    else:
        return []


def go_back(pre_graph, min_pos):
    """
    根据动态规划得到的前驱矩阵，重构read和参考基因组的异同。
    :param pre_graph: 动态规划得到的前驱矩阵
    :param min_pos: 回溯开始的位置
    :return: 以sam文件中对应的描述格式输出
    """
    cigar = []
    m, n = pre_graph.shape
    i, j = min_pos, n - 1

    # SAM文件描述会将连续的插入、删除等信息合并
    while i > 0 and j > 0:
        num = 0
        if pre_graph[i, j] == 1:  # 记录连续的插入
            while i > 0 and j > 0 and pre_graph[i, j] == 1:
                num += 1
                j -= 1
            cigar.append(str(num) + 'I')
        elif pre_graph[i, j] == 2:  # 处理连续的match和mismatch
            while i > 0 and j > 0 and pre_graph[i, j] == 2:
                num += 1
                i -= 1
                j -= 1
            cigar.append(str(num) + 'M')
        else:  # 处理连续的删除
            while i > 0 and j > 0 and pre_graph[i, j] == 3:
                num += 1
                i -= 1
            cigar.append(str(num) + 'D')

    # 判断最前端的插入和删除
    if i == 0 and j > 0:
        cigar.append(str(j) + 'I')
    elif j == 0 and i > 0:
        cigar.append(str(i) + 'D')
    cigar.reverse()
    return ''.join(cigar)


def align_graph(read, candidate):
    """
    这里通过动态规划算法计算两个序列建的编辑图，便于之后重构两个序列的差异
    :param read: 待比对片段
    :param candidate: 参考基因组的一个子串，预测该子串与read匹配
    :return: 三个值，第一个是比对结束时参考基因组子串上的位置，第二个类似编辑距离，第三个是前驱矩阵
    """
    m, n = len(candidate), len(read)
    d, pre = numpy.zeros((m + 1, n + 1)), numpy.zeros((m + 1, n + 1))  # 前驱矩阵和距离矩阵

    # 初始化前驱矩阵和距离矩阵
    for i in range(m):
        d[i, 0] = i
        pre[i, 0] = 3
    for i in range(n):
        d[0, i] = i
        pre[0, i] = 1

    # 使用动态规划算法计算编辑距离
    for i in range(1, m + 1, 1):
        for j in range(1, n + 1, 1):

            # 这里利用1表示右移，2表示右下移，3表示下移
            # 分为两种情况，一种是待比对的字符相等
            if candidate[i - 1] == read[j - 1]:
                if d[i - 1, j - 1] <= d[i, j - 1] + 1 and d[i - 1, j - 1] <= d[i - 1, j] + 1:
                    d[i, j] = d[i - 1, j - 1]
                    pre[i, j] = 2
                elif d[i - 1, j] < d[i - 1, j]:
                    d[i, j] = d[i - 1, j] + 1
                    pre[i, j] = 3
                else:
                    d[i, j] = d[i, j - 1] + 1
                    pre[i, j] = 1

            # 另一种是比对的字符不相等
            else:
                if d[i - 1, j - 1] <= d[i, j - 1] and d[i - 1, j - 1] <= d[i - 1, j]:
                    d[i, j] = d[i - 1, j - 1] + 1
                    pre[i, j] = 2
                elif d[i - 1, j] < d[i - 1, j]:
                    d[i, j] = d[i - 1, j] + 1
                    pre[i, j] = 3
                else:
                    d[i, j] = d[i, j - 1] + 1
                    pre[i, j] = 1

    # 找到回溯的起始位置
    min_pos = n
    for i in range(n, m + 1, 1):
        if d[i, n] < d[min_pos, n]:
            min_pos = i
    return min_pos, len(read) - d[min_pos, n], pre


def sum_of_flag(prp1, prp2, seq):
    """
    计算SAN文件中的第二列的数据，这里考虑其中的一部分情况
    :param prp1: pair-end中的其中一个read匹配到的位置
    :param prp2: pair-end中的其中另一个read匹配到的位置
    :param seq: 第一个参数是第几个read的值，1表示read1,2表示read2
    :return: flag值
    """
    flag = 1
    if len(prp1) == 1 and prp2 == 1:
        flag += 2
    if len(prp1) == 0:
        flag += 4
    if len(prp2) == 0:
        flag += 8
    if seq == 1:
        flag += 64
    if seq == 2:
        flag += 128
    return flag


def format_line(col_1, col_2, col_3, col_4, col_5, col_6, col_7, col_8, col_9, col_10, col_11):
    """
    根据提供的信息，构建出SAM文件中的一个结果行
    :param col_1: 第1列信息
    :param col_2: 第2列信息
    :param col_3: 第3列信息
    :param col_4: 第4列信息
    :param col_5: 第5列信息
    :param col_6: 第6列信息
    :param col_7: 第7列信息
    :param col_8: 第8列信息
    :param col_9: 第9列信息
    :param col_10: 第10列信息
    :param col_11: 第11列信息
    :return: 构建行的字符串
    """
    return col_1 + '\t' + str(col_2) + '\t' + col_3 + '\t' + str(col_4) + '\t' + str(
        col_5) + '\t' + col_6 + '\t' + col_7 + '\t' + str(col_8) + '\t' + str(
        col_9) + '\t' + col_10 + '\t' + col_11 + '\n'


def mate_info(prp1, prp2):
    """
    判断pair-end中另一个mate在参考基因组上的位置信息
    :param prp1: 当前read匹配到的位置
    :param prp2: mate read匹配到的位置
    :return: 两个。第一个是染色体的编号。这里只有两种，要么相同，要么没有。第二个是偏移位置，如果没有就为0
    """
    if len(prp1) == 1 and len(prp2) == 1:
        return '=', prp2[0]
    else:
        return '*', 0


def format_to_sam(read1, read2, prps1, prps2, reference, c_threshold):
    """
    将比对的结果以SAM格式存储在文件中
    :param read1: pair-end中的第一部分read
    :param read2: pair-end中的第一部分read
    :param prps1: pair-end中的第一部分read的比对结果
    :param prps2: pair-end中的第一部分read的比对结果
    :param reference: 参考基因组序列
    :param c_threshold: 插入删除允许的阈值
    """
    sam_content = list()  # 存储SAM文件的所有行
    sam_content.append('@HD VN:1.6 SO:coordinate\n')
    sam_content.append('@SQ SN:ref LN:45\n')
    for i in range(len(read1)):
        read_name = 'r' + str(i)
        flag = sum_of_flag(prps1[i], prps2[i], 1)
        if len(prps1[i]) == 1 and len(prps2[i]) == 1:
            length = prps2[i][0] - prps1[i][0]
            if length > 0:
                length += len(read2[i])
            else:
                length -= len(read1[i])
        else:
            length = 0
        if len(prps1[i]) == 0:
            sam_content.append(format_line(read_name, flag, '*', 0, 0, '*', '*', 0, 0, read1[i], '*'))
        else:
            mate_ref, mate_pos = mate_info(prps1[i], prps2[i])
            for p in prps1[i]:
                min_pos, quality, pre_graph = align_graph(read1[i], reference[p - 1:p + len(read1[i]) + c_threshold])
                cigar = go_back(pre_graph, min_pos)
                sam_content.append(
                    format_line(read_name, flag, 'ref', p, quality, cigar, mate_ref, mate_pos, length,
                                read1[i], '*'))
        flag = sum_of_flag(prps2[i], prps1[i], 2)
        if len(prps2[i]) == 0:
            sam_content.append(format_line(read_name, flag, '*', 0, 0, '*', '*', 0, 0, read2[i], '*'))
        else:
            mate_ref, mate_pos = mate_info(prps2[i], prps1[i])
            for p in prps2[i]:
                min_pos, quality, pre_graph = align_graph(read2[i], reference[p:p + len(read2[i]) + c_threshold])
                cigar = go_back(pre_graph, min_pos)
                sam_content.append(
                    format_line(read_name, flag, 'ref', p, quality, cigar, mate_ref, mate_pos, length,
                                read2[i], '*'))
    with open('../data/result.sam', 'w') as fp:
        fp.writelines(sam_content)


if __name__ == '__main__':
    k = 17
    h_max, d_max, c = 500, 10, 4
    start = time.time()
    kmer_index, reference = get_index('../data/NC_008253.fna', k)
    end = time.time()
    print(end - start)
    # kmer_index, reference = get_index('../data/test.txt', k)
    # read = load_reads('../data/read.txt')
    read1, read2 = load_reads('../data/Ecoli_4x.fq1', '../data/Ecoli_4x.fq2', 10000)
    prps1, prps2 = [], []
    # r1 = align(reference, kmer_index, read1[2889], k, h_max, d_max, c)
    # r2 = align(reference, kmer_index, read2[2889], k, h_max, d_max, c)
    # print(read2[2889])
    # print(reference[r2[0]-1:r2[0] + 70])
    for i in range(len(read1)):
        if i % 100 == 0:
            print(i)
        r1 = align(reference, kmer_index, read1[i], k, h_max, d_max, c)
        r2 = align(reference, kmer_index, read2[i], k, h_max, d_max, c)
        prps1.append(r1)
        prps2.append(r2)
    end = time.time()
    print(end - start)
    format_to_sam(read1, read2, prps1, prps2, reference, c)
    end = time.time()
    print(end - start)
