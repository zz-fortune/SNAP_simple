"""
将read和参考基因组进行比对

@author: zhang heng
"""
import copy
import numpy
from src.build_index import *


def seed_generation(seeds, read, kmer, interval, UT_kmer, UT_seq):
    for i in range(0, len(read), interval):

        # 判断该k-mer是否需要扩展
        extend = True
        uid, offset = UT_kmer[read[i:i + kmer]]
        for seed in seeds:
            if seed[3][i] is True and seed[3][i + kmer - 1] is True and seed[1] == uid and seed[2] == offset:
                extend = False
                break
        # 扩展k-mer生成seed
        if extend:
            unitig = UT_seq[uid - 1]
            start, end = i, i
            while start > -1 or end <= len(read) - kmer:
                if start > -1 and read[start:start + kmer] == unitig[
                                                              offset - i + start - 1:offset - i + start - 1 + kmer]:
                    start -= 1
                if end <= len(read) - kmer and read[end:end + kmer] == unitig[
                                                                       offset - i + start - 1:offset - i + end - 1 + kmer]:
                    end += 1
                if offset - i + start == 0 or offset - i + end:
                    pass


def unscore_most_hitting(seeds_hitting, scores):
    max_hitting, max_pos = -1, 0
    for k in seeds_hitting:
        if scores.get(k) is None and seeds_hitting[k] > max_hitting:
            max_hitting = seeds_hitting[k]
            max_pos = k
    return max_pos


def edit_distance(read, candidate, d_limit):
    m, n = len(candidate), len(read)
    d = numpy.zeros((m + 1, n + 1))
    for i in range(m):
        d[i, 0] = i
    for i in range(n):
        d[0, i] = i
    for i in range(1, m + 1, 1):
        for j in range(1, n + 1, 1):
            if candidate[i - 1] == read[j - 1] and d[i - 1, j - 1] <= d[i, j - 1] and d[i - 1, j - 1] <= d[i - 1, j]:
                d[i, j] = d[i - 1, j - 1]
            elif d[i - 1, j - 1] < d[i, j - 1] and d[i - 1, j - 1] < d[i - 1, j]:
                d[i, j] = d[i - 1, j - 1] + 1
            elif d[i, j - 1] < d[i - 1, j]:
                d[i, j] = d[i, j - 1] + 1
            else:
                d[i, j] = d[i - 1, j] + 1
    d_min = d[n, n]
    for i in range(n, m + 1, 1):
        if d[i, n] < d_min:
            d_min = d[i, n]
    if d_min > d_limit:
        return math.inf
    else:
        return d_min


def update_distance(d_best, d_second, d_cur, pos_best, pos):
    if d_cur < d_best:
        return copy.copy(d_cur), copy.copy(d_best), copy.copy(pos), copy.copy(pos_best)
    elif d_cur < d_second:
        return copy.copy(d_best), copy.copy(d_cur), copy.copy(pos_best), copy.copy(pos)


def align(reference, kmer_index, read, k_mer, h_max, d_max, c_threshold):
    non_overlapping_num = 0
    d_best, d_second = math.inf, math.inf
    pos_best, pos_second = 0, 0
    seeds_hitting = dict()
    scores = dict()
    for i in range(len(read) - k_mer + 1):
        if i % k_mer == 0:
            non_overlapping_num += 1
        seed = read[i:i + k_mer]
        locations = (kmer_index.get(seed) or [])
        if len(locations) <= h_max:
            for l in locations:
                pos = l - i
                if pos > 0:
                    seeds_hitting[pos] = (seeds_hitting.get(pos) or 0) + 1
        pos = unscore_most_hitting(seeds_hitting, scores)
        if pos == 0:
            continue
        if d_best > d_max:
            d_limit = d_max + c_threshold - 1
        elif d_second >= d_best + c_threshold:
            d_limit = d_best + c_threshold - 1
        else:
            d_limit = d_best - 1
        d_cur = edit_distance(read, reference[pos - 1:pos + len(read) + c_threshold], d_limit)
        d_best, d_second, pos_best, pos_second = update_distance(d_best, d_second, d_cur, pos_best, pos)
        scores[pos] = d_cur
        if d_best < c_threshold and d_second < d_best + c_threshold:
            return [pos_best, pos_second]
        elif non_overlapping_num >= d_best + c_threshold:
            for k in seeds_hitting.keys():
                if scores.get(k) is None:
                    scores[k] = edit_distance(read, reference[k - 1:k + len(read) + c_threshold], d_limit)
            break
    if d_best <= d_max and d_second >= d_best + c_threshold:
        return [pos_best]
    elif len(seeds_hitting) == 0:
        for i in range(len(read) - k_mer):
            if kmer_index.get(read[i:i + k_mer]):
                return kmer_index[read[:k_mer]]
        return []
    elif d_best <= d_max:
        result = []
        for k in scores.keys():
            if scores.get(k) <= d_max:
                result.append(k)
        return result
    else:
        return []


def go_back(pre, end_pos, re):
    pass


def align_graph(read, reference, prp, c_threshold):
    m, n = len(read) + c_threshold, len(read)
    d = numpy.zeros((m + 1, n + 1))
    pre = numpy.zeros((m + 1, n + 1))
    for i in range(m):
        d[i, 0] = i
        pre[i, 0] = 2
    for i in range(n):
        d[0, i] = i
        pre[0, i] = 1
    for k in range(len(prp)):
        candidate = reference[prp[k]:prp[k] + m]
        for i in range(1, m + 1, 1):
            for j in range(1, n + 1, 1):
                if candidate[i - 1] == read[j - 1] and d[i - 1, j - 1] <= d[i, j - 1] and d[i - 1, j - 1] <= d[
                    i - 1, j]:
                    d[i, j] = d[i - 1, j - 1]
                    pre[i, j] = 3
                elif d[i - 1, j - 1] < d[i, j - 1] and d[i - 1, j - 1] < d[i - 1, j]:
                    d[i, j] = d[i - 1, j - 1] + 1
                    pre[i, j] = 3
                elif d[i, j - 1] < d[i - 1, j]:
                    d[i, j] = d[i, j - 1] + 1
                    pre[i, j] = 1
                else:
                    d[i, j] = d[i - 1, j] + 1
                    pre[i, j] = 2
        min_pos = n
        for i in range(n, m + 1, 1):
            if d[i, n] < d[min_pos, n]:
                min_pos = i
        re = []
        go_back(pre, min_pos, re)


if __name__ == '__main__':
    k = 17
    h_max, d_max, c = 500, 10, 4
    kmer_index, reference = get_index('../data/NC_008253.fna', k)
    # kmer_index, reference = get_index('../data/test.txt', k)
    # read = load_reads('../data/read.txt')
    read = load_reads('../data/Ecoli_4x.fq1', 10)
    r = align(reference, kmer_index, read[0], k, h_max, d_max, c)
    print(read[0])
    print(reference[r[0] - 1:r[0] + len(read[0])])
    print(r)
