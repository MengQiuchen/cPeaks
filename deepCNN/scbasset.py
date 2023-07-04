
import os
import torch
import torch.nn as nn
import numpy as np
import torch.nn.functional as F
import torch.nn as nn
from typing import Callable
import torch.optim as optim
from torch.utils.data import TensorDataset
from torch.utils.data import random_split
from torch.utils.data import DataLoader
from tqdm import tqdm
import matplotlib.pyplot as plt


print('import done!')


# In[3]:
"""
def binary_accuracy_and_recall(output, target):
    # output和target都是大小为(batch_size, num_labels)的张量
    with torch.no_grad():
        # 计算准确率
        pred = (output > 0.5).int()
        correct = pred.eq(target.int()).sum(dim=0)  # 按列相加，得到长度为num_labels的张量
        acc = correct.float() / target.size(0)

        # 计算召回率
        true_positives = (pred & target.int()).sum(dim=0)
        actual_positives = target.int().sum(dim=0)
        recall = true_positives.float() / actual_positives.float()
        
    return acc, recall
"""

def setup_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    #random.seed(seed)
    torch.backends.cudnn.deterministic = True
# 设置随机数种子
setup_seed(2023)

# net

import math
from typing import List
from typing import Callable

def _round(x):
    return int(np.round(x))



class ConvLayer(nn.Module):
    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        kernel_size: int,
        pool_size: int = None,
        batch_norm: bool = True,
        dropout: float = 0.0,
        activation_fn: Callable = nn.GELU(),
    ):
        super().__init__()
        self.conv = nn.Conv1d(
            in_channels=in_channels,
            out_channels=out_channels,
            kernel_size=kernel_size,
            padding=kernel_size // 2,
        )
        self.batch_norm = (
            nn.BatchNorm1d(out_channels) if batch_norm else nn.Identity()
        )
        self.pool = (
            nn.MaxPool1d(pool_size) if pool_size is not None else nn.Identity()
        )
        self.activation_fn = activation_fn
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor):
        x = self.conv(x)
        x = self.batch_norm(x)
        x = self.dropout(x)
        x = self.pool(x)
        x = self.activation_fn(x)
        return x


class DenseLayer(nn.Module):
    def __init__(
        self,
        in_features: int,
        out_features: int,
        use_bias: bool = True,
        batch_norm: bool = True,
        dropout: float = 0.2,
        activation_fn: Callable = nn.GELU(),
    ):
        super().__init__()
        self.dense = nn.Linear(in_features, out_features, bias=use_bias)
        self.batch_norm = (
            nn.BatchNorm1d(out_features) if batch_norm else nn.Identity()
        )
        self.activation_fn = activation_fn
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor):
        x = self.dense(x)
        x = self.batch_norm(x)
        x = self.dropout(x)
        x = self.activation_fn(x)
        return x

    




def _col_round(x):
    frac = x - math.floor(x)
    if frac <= 0.5:
        return math.floor(x)
    return math.ceil(x)


def _get_filter_dim(seq_length: int, pooling_sizes: List[int]):
    filter_dim = seq_length
    for ps in pooling_sizes:
        filter_dim = _col_round(filter_dim / ps)
    return filter_dim


class scBasset(nn.Module):
    """
    PytTorch implementation of scBasset model (Yuan and Kelley, 2022)
    Article link: https://www.nature.com/articles/s41592-022-01562-8
    Original implementation in Keras: https://github.com/calico/scBasset
    Args:
        n_cells: number of cells to predict region accessibility
        n_filters_init: nr of filters for the initial conv layer
        n_repeat_blocks_tower: nr of layers in the convolutional tower
        filters_mult: proportion by which the nr of filters should inrease in the
            convolutional tower
        n_bottleneck_layer: size of the bottleneck layer
        batch_norm: whether to apply batch norm across model layers
        dropout: dropout rate across layers, by default we don't do it for
            convolutional layers but we do it for the dense layers
    """

    def __init__(
        self,
        n_cells: int,
        n_filters_init: int = 288,
        n_repeat_blocks_tower: int = 5,
        filters_mult: float = 1.122,
        n_filters_pre_bottleneck: int = 256,
        n_bottleneck_layer: int = 32,
        batch_norm: bool = True,
        dropout: float = 0.0,
        genomic_seq_length: int = 2000,
    ):
        super().__init__()

        self.stem = ConvLayer(
            in_channels=4,
            out_channels=n_filters_init,
            kernel_size=17,
            pool_size=3,
            dropout=dropout,
            batch_norm=batch_norm,
        )

        tower_layers = []
        curr_n_filters = n_filters_init
        for i in range(n_repeat_blocks_tower):
            tower_layers.append(
                ConvLayer(
                    in_channels=curr_n_filters,
                    out_channels=_round(curr_n_filters * filters_mult),
                    kernel_size=5,
                    pool_size=2,
                    dropout=dropout,
                    batch_norm=batch_norm,
                )
            )
            curr_n_filters = _round(curr_n_filters * filters_mult)
        self.tower = nn.Sequential(*tower_layers)

        self.pre_bottleneck = ConvLayer(
            in_channels=curr_n_filters,
            out_channels=n_filters_pre_bottleneck,
            kernel_size=1,
            dropout=dropout,
            batch_norm=batch_norm,
            pool_size=1,
        )

        # get pooling sizes of the upstream conv layers
        pooling_sizes = [3] + [2] * n_repeat_blocks_tower + [1]
        # get filter dimensionality to account for variable sequence length
        filter_dim = _get_filter_dim(
            seq_length=genomic_seq_length, pooling_sizes=pooling_sizes
        )
        self.bottleneck = DenseLayer(
            in_features=n_filters_pre_bottleneck * filter_dim,
            out_features=n_bottleneck_layer,
            use_bias=True,
            batch_norm=True,
            dropout=0.2,
            activation_fn=nn.Identity(),
        )
        self.final = nn.Linear(n_bottleneck_layer, n_cells)
        self.softmax = nn.Softmax(dim=1)

    def forward(
        self,
        x: torch.Tensor,  # input shape: (batch_size, 4, seq_length)
    ):
        # TODO: add random shift to act as a regularizer on the dataset level
        # TODO: add use reverse complement randomly on the dataset level
        x = self.stem(x)
        x = self.tower(x)
        x = self.pre_bottleneck(x)
        # flatten the input
#         print(x.size())
        x = x.view(x.shape[0], -1)
#         print(x.size())
        x = self.bottleneck(x)
#         print(x.size())
        x = self.final(x)
        x = self.softmax(x)
        return x





# import torch
# import torch.nn as nn

# from layers import ConvLayer, DenseLayer
# from utils import _round


# class scBasset(nn.Module):
#     """
#     PytTorch implementation of scBasset model (Yuan and Kelley, 2022)
#     Article link: https://www.nature.com/articles/s41592-022-01562-8
#     Original implementation in Keras: https://github.com/calico/scBasset

#     Args:
#         n_cells: number of cells to predict region accessibility
#         n_filters_init: nr of filters for the initial conv layer
#         n_repeat_blocks_tower: nr of layers in the convolutional tower
#         filters_mult: proportion by which the nr of filters should inrease in the
#             convolutional tower
#         n_bottleneck_layer: size of the bottleneck layer
#         batch_norm: whether to apply batch norm across model layers
#         dropout: dropout rate across layers, by default we don't do it for
#             convolutional layers but we do it for the dense layers
#     """

#     def __init__(
#         self,
#         n_cells: int,
#         n_filters_init: int = 288,
#         n_repeat_blocks_tower: int = 5,
#         filters_mult: float = 1.122,
#         n_filters_pre_bottleneck: int = 256,
#         n_bottleneck_layer: int = 32,
#         batch_norm: bool = True,
#         dropout: float = 0.0,
#     ):
#         super().__init__()

#         self.stem = ConvLayer(
#             in_channels=4,
#             out_channels=n_filters_init,
#             kernel_size=17,
#             pool_size=3,
#             dropout=dropout,
#             batch_norm=batch_norm,
#         )

#         tower_layers = []
#         curr_n_filters = n_filters_init
#         for i in range(n_repeat_blocks_tower):
#             tower_layers.append(
#                 ConvLayer(
#                     in_channels=curr_n_filters,
#                     out_channels=_round(curr_n_filters * filters_mult),
#                     kernel_size=5,
#                     pool_size=2,
#                     dropout=dropout,
#                     batch_norm=batch_norm,
#                 )
#             )
#             curr_n_filters = _round(curr_n_filters * filters_mult)
#         self.tower = nn.Sequential(*tower_layers)

#         self.pre_bottleneck = ConvLayer(
#             in_channels=curr_n_filters,
#             out_channels=n_filters_pre_bottleneck,
#             kernel_size=1,
#             dropout=dropout,
#             batch_norm=batch_norm,
#             pool_size=2,
#         )
#         self.bottleneck = DenseLayer(
#             in_features=n_filters_pre_bottleneck * 7,
#             out_features=n_bottleneck_layer,
#             use_bias=True,
#             batch_norm=True,
#             dropout=0.2,
#             activation_fn=nn.Identity(),
#         )
#         self.final = nn.Linear(n_bottleneck_layer, n_cells)
#         self.sigmoid = nn.Sigmoid()

#     def forward(
#         self,
#         x: torch.Tensor,  # input shape: (batch_size, 4, seq_length)
#     ):
#         # TODO: add random shift to act as a regularizer on the dataset level
#         # TODO: add use reverse complement randomly on the dataset level
#         x = self.stem(x)
#         x = self.tower(x)
#         x = self.pre_bottleneck(x)
#         # flatten the input
#         x = x.view(x.shape[0], -1)
#         x = self.bottleneck(x)
#         x = self.final(x)
#         return self.sigmoid(x)
