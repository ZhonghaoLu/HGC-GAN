import copy
import torch
import random
import numpy as np
import torch.nn as nn
from get_adj import adj
from sklearn import metrics
from torch.autograd import Variable
from sklearn.metrics import roc_auc_score, precision_recall_curve


def select_negative_items(Data, num_pm, num_zr, disease_num):
    data = np.array(Data)
    n_items_pm = np.zeros_like(data)
    n_items_zr = np.zeros_like(data)
    for i in range(data.shape[0]):
        p_items = np.where(data[i] != 0)[0]
        all_item_index = random.sample(range(data.shape[1]), disease_num)
        for j in range(p_items.shape[0]):
            all_item_index.remove(list(p_items)[j])
        random.shuffle(all_item_index)
        n_item_index_pm = all_item_index[0: num_pm]
        n_item_index_zr = all_item_index[num_pm: (num_pm + num_zr)]
        n_items_pm[i][n_item_index_pm] = 1
        n_items_zr[i][n_item_index_zr] = 1
    return n_items_pm, n_items_zr


def main(lncrna_num, disease_num, epochCount, pro_ZR, pro_PM, alpha, batchSize,
         input_net, true_input_net, test_input_net, test_negative,
         Noise_LMDN, Noise_LMDN_h, True_LMDN, True_LMDN_h,
         G, D, G_step, D_step, G_PATH, D_PATH
         ):
    best_auc = 0
    t = 0

    regularization = nn.MSELoss()
    d_optimizer = torch.optim.RMSprop(D.parameters(), lr=0.0001)
    g_optimizer = torch.optim.RMSprop(G.parameters(), lr=0.0001)

    Noise_LDA = adj(lncrna_num, disease_num, input_net)
    True_LDA = adj(lncrna_num, disease_num, true_input_net)
    Test_Adj = adj(lncrna_num, disease_num, test_input_net)

    for epoch in range(epochCount):
        for step in range(G_step):
            leftIndex = random.randint(0, lncrna_num - batchSize - 1)
            realData = Variable(copy.deepcopy(True_LDA[leftIndex:leftIndex + batchSize]))
            noiseData = Variable(copy.deepcopy(Noise_LDA[leftIndex:leftIndex + batchSize]))

            e_i = Variable(copy.deepcopy(Noise_LDA[leftIndex:leftIndex + batchSize]))
            n_items_pm, n_items_zr = select_negative_items(noiseData, pro_PM, pro_ZR, disease_num)
            k_i_zp = Variable(torch.tensor(n_items_pm + n_items_zr))
            realData_zp = Variable(torch.ones_like(realData)) * e_i + Variable(torch.zeros_like(realData)) * k_i_zp

            fake_embeding, r_i = G(Noise_LMDN, Noise_LMDN_h, noiseData, batchSize, leftIndex)
            pred_matrix = r_i * (e_i + k_i_zp)
            fakeData_result = D(pred_matrix, fake_embeding)

            g_loss = np.mean((1. - fakeData_result.detach().numpy() + 10e-5)) + alpha * regularization(pred_matrix,
                                                                                                       realData_zp)
            g_optimizer.zero_grad()
            g_loss.backward(retain_graph=True)
            g_optimizer.step()

        for step in range(D_step):
            leftIndex = random.randint(1, lncrna_num - batchSize - 1)
            realData = Variable(copy.deepcopy(True_LDA[leftIndex:leftIndex + batchSize]))
            noise_Data = Variable(copy.deepcopy(Noise_LDA[leftIndex:leftIndex + batchSize]))

            e_i = Variable(copy.deepcopy(Noise_LDA[leftIndex:leftIndex + batchSize]))
            n_items_pm, _ = select_negative_items(noise_Data, pro_PM, pro_ZR, disease_num)
            k_i = Variable(torch.tensor(n_items_pm))

            fake_embeding, r_i = G(Noise_LMDN, Noise_LMDN_h, noise_Data, batchSize, leftIndex)
            pred_matrix = r_i * (e_i + k_i)
            fakeData_result = D(pred_matrix, fake_embeding)

            true_embeding, _ = G(True_LMDN, True_LMDN_h, realData, batchSize, leftIndex)
            realData_result = D(realData, true_embeding)

            d_loss = -np.mean((realData_result.detach().numpy() + 10e-5) +
                              (1. - fakeData_result.detach().numpy() + 10e-5)) + 0 * regularization(pred_matrix,
                                                                                                    realData_zp)
            d_optimizer.zero_grad()
            d_loss.backward(retain_graph=True)
            d_optimizer.step()

            for p in D.parameters(): p.data.clamp_(-0.01, 0.01)

        label = []
        pred = []
        for testUser in range(len(Test_Adj)):
            data = Variable(copy.deepcopy(Noise_LDA[testUser:testUser + 1]))
            _, predData = G(Noise_LMDN, Noise_LMDN_h, data, 1, testUser)
            pred_i = predData[0].tolist()
            test_i = Test_Adj[testUser].tolist()
            for i in range(len(test_i)):
                if test_i[i] == 1 and [testUser, i] in test_negative:
                    label.append(0)
                    pred.append(pred_i[i])
                if test_i[i] == 1 and [testUser, i] not in test_negative:
                    label.append(1)
                    pred.append(pred_i[i])

        auc = roc_auc_score(label, pred)
        p, r, thresholds = precision_recall_curve(label, pred)
        print('Epoch[{}/{}],AUC:{:0.4f},Aupr:{:0.4f}'.format(epoch, epochCount, auc, metrics.auc(r, p)))

        if auc > best_auc:
            best_auc = auc
            torch.save(G, G_PATH)
            torch.save(D, D_PATH)

        if auc < best_auc:
            t -= 1
        else:
            t = 0
        if t <= -10: break

    return best_auc
