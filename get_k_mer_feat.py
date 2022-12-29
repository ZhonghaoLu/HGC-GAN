import torch

def get_basic_g(n, k, letters, m='', k_letter=''):
    if m == '':
        m = n - k
        k_letter = list(letters[:n])
    if n == m + 1: return k_letter
    num_of_perm = []
    for perm in get_basic_g(n - 1, k, letters, m, k_letter):
        temp_letter = [i for i in k_letter]
        num_of_perm += [perm + i for i in temp_letter]
    return num_of_perm

def get_k_mer(xulie):
    k = 3
    a = 'A'
    g = 'G'
    c = 'C'
    t = 'T'
    basic = ['A', 'C', 'G', 'T']
    basic_group = get_basic_g(len(basic), k, basic)
    base_group_count = list(dict(zip(basic_group, [0 for _ in range(len(basic_group))])) for _ in range(len(xulie)))

    for i in range(len(xulie)):
        for j in range(len(xulie[i])-k+1):
            base_group_count[i][xulie[i][j:j+k]] += 1

    feat_f = [list(base_group_count[i].values()) for i in range(len(base_group_count))]
    feat = []
    for i in range(len(feat_f)):
        sum_f = 0
        for j in range(len(basic_group)):
            sum_f += feat_f[i][j]
        feat.append([round(feat_f[i][n] / sum_f, 4) for n in range(len(feat_f[0]))])
    return feat

def get_feature(lncrna_name):

    f = open("./Data/Dataset2/lncRNA_Sequence.txt", 'r', encoding='utf-8')
    contents = f.readlines()
    lncrna_name_seq = []

    for content in contents:
        value = content.split()
        lncrna_name_seq.append(value)
    f.close()

    lncrna_xulie = []
    for name in lncrna_name:
        for xuelie in lncrna_name_seq:
            if xuelie[0] == name:
                lncrna_xulie.append(xuelie[1])
                break

    k_mer = get_k_mer(lncrna_xulie)
    k_mer = torch.tensor(k_mer)

    print("The k_mer feature construction of the node is completed!")
    return k_mer