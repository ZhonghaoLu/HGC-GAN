'''
python3.7
-*- coding: utf-8 -*-

Copyright (C) 2023 Zhonghao Lu, Inc. All Rights Reserved

@Time    : 2023/1/1
@Author  : Zhonghao Lu
@Email   : lzhao6688@163.com
@File    : main.py
@Software: PyCharm
@Description:
Pytorch Implementation of Heterogeneous Graph Convolutional Generative Adversarial Network (HGC-GAN) model in:
Zhonghao Lu et al. Predicting lncRNA-disease associations based on heterogeneous graph convolutional generative adversarial network
'''

import torch
import train
import random
import numpy as np
from Model import model
from sklearn.model_selection import KFold

import get_k_mer_feat
from Data.load_Dataset1 import dataset1
from Data.load_Dataset2 import dataset2
from Data.load_Dataset3 import dataset3
from Data.load_Dataset4 import dataset4
from hetero_graph import get_hetero_graph

import warnings

warnings.filterwarnings("ignore")

if __name__ == '__main__':
    epochs = 150
    pro_ZR = 50
    pro_PM = 50
    alpha = 0.1

    feat_shape = 64
    out_feat = 64

    G_step = 5
    D_step = 2
    batchSize = 32

    G_PATH = './weights/G.pth'
    D_PATH = './weights/D.pth'

    AUC = []

    lncrna_num, mirna_num, disease_num, lncrna_disease, mirna_disease, lncrna_mirna = dataset1()
    # lncrna_num, mirna_num, disease_num, lncrna_disease, mirna_disease, lncrna_mirna, lncRNA_name = dataset2()
    # lncrna_num, mirna_num, disease_num, lncrna_disease, mirna_disease, lncrna_mirna = dataset3()
    # lncrna_num, mirna_num, disease_num, lncrna_disease, mirna_disease, lncrna_mirna = dataset4()

    # lncrna_feat = get_k_mer_feat.get_feature(lncRNA_name)
    lncrna_feat = torch.rand(lncrna_num, feat_shape)
    mirna_feat = torch.rand(mirna_num, feat_shape)
    disease_feat = torch.rand(disease_num, feat_shape)

    print("\nData Details:")
    print("lncRNA:{}  miRNA:{}  Disease:{}".format(lncrna_num, mirna_num, disease_num))
    print("lncRNA-Disease:{}  miRNA-Disease:{}  lncRNA-miRNA:{}".format(len(lncrna_disease[0]), len(mirna_disease[0]),
                                                                        len(lncrna_mirna[0])))
    print("Sparsity of lncRNA-Disease associated data: {:.6f}\n".format(
        len(lncrna_disease[0]) / (lncrna_num * disease_num)))

    edge_data = []
    for i in range(len(lncrna_disease[0])):
        x = []
        x.append(lncrna_disease[0][i])
        x.append(lncrna_disease[1][i])
        edge_data.append(x)
    edge_data = np.array(edge_data)

    negativeSample_edge = []
    for i in range(len(edge_data)):
        row = random.randint(0, lncrna_num - 1)
        col = random.randint(0, disease_num - 1)
        while ([row, col] in edge_data.tolist() or [row, col] in negativeSample_edge):
            row = random.randint(0, lncrna_num - 1)
            col = random.randint(0, disease_num - 1)
        negativeSample_edge.append([row, col])

    kf = KFold(n_splits=10, shuffle=True)
    for train_index, test_index in kf.split(edge_data):
        train_negative = random.sample(negativeSample_edge, int(len(train_index)))
        test_negative = [data_negative for data_negative in negativeSample_edge if data_negative not in train_negative]
        train_index = train_index.tolist()
        train_lncrna_20 = random.sample(train_index, int(len(train_index) * 0.25))
        train_lncrna_60 = [i for i in train_index if i not in train_lncrna_20]
        test_lncrna = test_index.tolist()

        input_net = [[], []]
        for i in train_lncrna_60:
            input_net[0].append(edge_data[i][0])
            input_net[1].append(edge_data[i][1])
        for i in range(len(train_negative)):
            input_net[0].append(train_negative[i][0])
            input_net[1].append(train_negative[i][1])

        true_input_net = [[], []]
        for i in train_index:
            true_input_net[0].append(edge_data[i][0])
            true_input_net[1].append(edge_data[i][1])

        test_input_net = [[], []]
        for i in test_lncrna:
            test_input_net[0].append(edge_data[i][0])
            test_input_net[1].append(edge_data[i][1])
        for i in range(len(test_negative)):
            test_input_net[0].append(test_negative[i][0])
            test_input_net[1].append(test_negative[i][1])

        G = model.generator(disease_num, feat_shape, out_feat)
        D = model.discriminator(disease_num, feat_shape, out_feat)
        Noise_LMDN, Noise_LMDN_h, True_LMDN, True_LMDN_h = get_hetero_graph(input_net, true_input_net, lncrna_mirna,
                                                                            mirna_disease, lncrna_feat, mirna_feat,
                                                                            disease_feat)

        auc = train.main(lncrna_num, disease_num, epochs, pro_ZR, pro_PM, alpha, batchSize,
                         input_net, true_input_net, test_input_net, test_negative,
                         Noise_LMDN, Noise_LMDN_h, True_LMDN, True_LMDN_h,
                         G, D, G_step, D_step, G_PATH, D_PATH)
        AUC.append(auc)
        print('\n10_fold_auc:{}'.format(AUC))
        print('Mean auc:{}\n'.format(sum(AUC) / len(AUC)))
