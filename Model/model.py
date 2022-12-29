import copy
import torch
import torch.nn as nn
import dgl.nn.pytorch as dglnn
from torch.autograd import Variable


class discriminator(nn.Module):
    def __init__(self, itemCount, feat_shape, out_feat_shape):
        super(discriminator, self).__init__()

        self.itemCount = itemCount
        self.feat_shape = feat_shape
        self.out_feat_shape = out_feat_shape

        self.f1 = nn.Linear(self.itemCount + self.out_feat_shape, 1024)
        self.a1 = nn.ReLU(True)
        self.f2 = nn.Linear(1024, 128)
        self.a2 = nn.ReLU(True)
        self.f3 = nn.Linear(128, 16)
        self.a3 = nn.ReLU(True)
        self.f4 = nn.Linear(16, 1)
        self.a4 = nn.Sigmoid()

    def forward(self, Adj, embeding):
        x = torch.cat((Adj, embeding), 1)
        x = self.f1(x)
        x = self.a1(x)
        x = self.f2(x)
        x = self.a2(x)
        x = self.f3(x)
        x = self.a3(x)
        x = self.f4(x)
        return x


class generator(nn.Module):
    def __init__(self, itemCount, feat_shape, out_feat_shape):
        super(generator, self).__init__()

        self.itemCount = itemCount
        self.feat_shape = feat_shape
        self.out_feat_shape = out_feat_shape

        self.HeteroConv1 = dglnn.HeteroGraphConv({
            'LvsD': dglnn.GraphConv(self.feat_shape, 64),
            'DvsL': dglnn.GraphConv(self.feat_shape, 64),
            'LvsM': dglnn.GraphConv(self.feat_shape, 64),
            'MvsL': dglnn.GraphConv(self.feat_shape, 64),
            'DvsM': dglnn.GraphConv(self.feat_shape, 64),
            'MvsD': dglnn.GraphConv(self.feat_shape, 64)
        },
            aggregate='sum')

        self.HeteroConv2 = dglnn.HeteroGraphConv({
            'LvsD': dglnn.GraphConv(64, 64),
            'DvsL': dglnn.GraphConv(64, 64),
            'LvsM': dglnn.GraphConv(64, 64),
            'MvsL': dglnn.GraphConv(64, 64),
            'DvsM': dglnn.GraphConv(64, 64),
            'MvsD': dglnn.GraphConv(64, 64)
        },
            aggregate='sum')

        self.HeteroConv3 = dglnn.HeteroGraphConv({
            'LvsD': dglnn.GraphConv(64, self.out_feat_shape),
            'DvsL': dglnn.GraphConv(64, self.out_feat_shape),
            'LvsM': dglnn.GraphConv(64, self.out_feat_shape),
            'MvsL': dglnn.GraphConv(64, self.out_feat_shape),
            'DvsM': dglnn.GraphConv(64, self.out_feat_shape),
            'MvsD': dglnn.GraphConv(64, self.out_feat_shape)
        },
            aggregate='sum')

        self.f1 = nn.Linear(self.itemCount + self.out_feat_shape, 256)
        self.a1 = nn.ReLU(True)
        self.f2 = nn.Linear(256, 512)
        self.a2 = nn.ReLU(True)
        self.f3 = nn.Linear(512, 1024)
        self.a3 = nn.ReLU(True)
        self.f4 = nn.Linear(1024, itemCount)
        self.a4 = nn.Sigmoid()

    def forward(self, g, h, Adj, size, leftIndex):
        h1 = self.HeteroConv1(g, h)
        h2 = self.HeteroConv2(g, h1)
        h3 = self.HeteroConv3(g, h2)

        lncrna_feat = h3.get('LncRNA')
        lncrna_feat = torch.nn.functional.normalize(lncrna_feat, p=1, dim=1)
        lncrna_feat = torch.tensor(lncrna_feat)

        fake_embeding = Variable(copy.deepcopy(lncrna_feat[leftIndex:leftIndex + size]))

        M = torch.cat([Adj, fake_embeding], 1)

        x = self.f1(M)
        x = self.a1(x)
        x = self.f2(x)
        x = self.a2(x)
        x = self.f3(x)
        x = self.a3(x)
        x = self.f4(x)
        x = self.a4(x)

        return fake_embeding, x
