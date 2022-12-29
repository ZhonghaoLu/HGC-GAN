import torch
import numpy as np

def adj(row, col,data):
    Adj = np.zeros((row, col))
    for i in range(len(data[0])):
        Adj[data[0][i]][data[1][i]] = 1
    Vector = torch.tensor(Adj.astype(np.float32))
    return Vector