import os
import dgl
import torch
from dgl import save_graphs, load_graphs


def get_hetero_graph(input_net, true_input_net, lncrna_mirna, mirna_disease, lncrna_feat, mirna_feat, disease_feat):
    Noise_LMDN_data = {
        ('LncRNA', 'LvsD', 'Disease'): (torch.tensor(input_net[0]), torch.tensor(input_net[1])),
        ('Disease', 'DvsL', 'LncRNA'): (torch.tensor(input_net[1]), torch.tensor(input_net[0])),
        ('LncRNA', 'LvsM', 'MiRNA'): (torch.tensor(lncrna_mirna[0]), torch.tensor(lncrna_mirna[1])),
        ('MiRNA', 'MvsL', 'LncRNA'): (torch.tensor(lncrna_mirna[1]), torch.tensor(lncrna_mirna[0])),
        ('MiRNA', 'MvsD', 'Disease'): (torch.tensor(mirna_disease[0]), torch.tensor(mirna_disease[1])),
        ('Disease', 'DvsM', 'MiRNA'): (torch.tensor(mirna_disease[1]), torch.tensor(mirna_disease[0]))
    }
    Noise_LMDN = dgl.heterograph(Noise_LMDN_data)
    Noise_LMDN.nodes['LncRNA'].data['feat'] = lncrna_feat
    Noise_LMDN.nodes['MiRNA'].data['feat'] = mirna_feat
    Noise_LMDN.nodes['Disease'].data['feat'] = disease_feat
    Noise_LMDN_h = {'LncRNA': Noise_LMDN.nodes['LncRNA'].data['feat'],
                    'MiRNA': Noise_LMDN.nodes['MiRNA'].data['feat'],
                    'Disease': Noise_LMDN.nodes['Disease'].data['feat']}


    True_LMDN_data = {
        ('LncRNA', 'LvsD', 'Disease'): (torch.tensor(true_input_net[0]), torch.tensor(true_input_net[1])),
        ('Disease', 'DvsL', 'LncRNA'): (torch.tensor(true_input_net[1]), torch.tensor(true_input_net[0])),
        ('LncRNA', 'LvsM', 'MiRNA'): (torch.tensor(lncrna_mirna[0]), torch.tensor(lncrna_mirna[1])),
        ('MiRNA', 'MvsL', 'LncRNA'): (torch.tensor(lncrna_mirna[1]), torch.tensor(lncrna_mirna[0])),
        ('MiRNA', 'MvsD', 'Disease'): (torch.tensor(mirna_disease[0]), torch.tensor(mirna_disease[1])),
        ('Disease', 'DvsM', 'MiRNA'): (torch.tensor(mirna_disease[1]), torch.tensor(mirna_disease[0]))
    }
    True_LMDN = dgl.heterograph(True_LMDN_data)
    True_LMDN.nodes['LncRNA'].data['feat'] = lncrna_feat
    True_LMDN.nodes['MiRNA'].data['feat'] = mirna_feat
    True_LMDN.nodes['Disease'].data['feat'] = disease_feat
    True_LMDN_h = {'LncRNA': True_LMDN.nodes['LncRNA'].data['feat'],
                   'MiRNA': True_LMDN.nodes['MiRNA'].data['feat'],
                   'Disease': True_LMDN.nodes['Disease'].data['feat']}

    # Save hetero graph
    # save_path = './save_hetero_graph'
    # Noise_LMDN_mode = 'Noise_LMDN'
    # Noise_LMDN_path = os.path.join(save_path, Noise_LMDN_mode + '.bin')
    # save_graphs(Noise_LMDN_path, Noise_LMDN ,)
    # True_LMDN_mode = 'True_LMDN'
    # True_LMDN_path = os.path.join(save_path, True_LMDN_mode + '.bin')
    # save_graphs(True_LMDN_path, True_LMDN,)
    # Load hetero graph
    # Noise_LMDN, _ = load_graphs(Noise_LMDN_path)
    # True_LMDN, _ = load_graphs(True_LMDN_path)

    return Noise_LMDN, Noise_LMDN_h, True_LMDN, True_LMDN_h
