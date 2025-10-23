import os,torch
import os.path
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader,Dataset
import numpy as np
from torch.utils.data import DataLoader
from torch.utils.data.dataset import random_split, Subset
import torch.nn as nn
import torch.nn.functional as F
import os.path
import torch

class onetransformer(nn.Module):
    def __init__(self,vocab_size=33,embed_size = 20, hidden_size=128,num_classes=2,d_model=512,num_layers=1,multitask=False,addseq=True):
        super(onetransformer, self).__init__()
        self.multitask = multitask
        self.addseq = addseq
        if addseq:
            self.embedding = nn.Embedding(vocab_size, embed_size)
            self.esm_linear = nn.Linear(1280,d_model-embed_size)
        else:
            self.esm_linear = nn.Linear(1280,d_model)
        # self.dropout = nn.Dropout(0.1)
        dropout_rate = 0.2
        encoder_layer = nn.TransformerEncoderLayer(d_model=d_model, nhead=8)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)

        self.relu = nn.ReLU()
        # self.fc1 = nn.Linear(d_model, hidden_size, bias=True)
        # self.fc2 = nn.Linear(hidden_size, num_classes, bias=True)

        self.linear1 = nn.Sequential(
            nn.Linear(d_model, hidden_size),
            nn.ReLU(),
            nn.Dropout(dropout_rate)
        )

        self.linear2 = nn.Sequential(
            nn.Linear(hidden_size, 2),
            nn.Sigmoid(),
            nn.Dropout(dropout_rate)
        )
        if multitask:
            self.linear1_othertask = nn.Sequential(
                nn.Linear(d_model, hidden_size),
                nn.ReLU(),
                nn.Dropout(dropout_rate)
            )
    
            self.linear2_othertask = nn.Sequential(
                nn.Linear(hidden_size, 2),
                nn.Sigmoid(),
                nn.Dropout(dropout_rate)
            )
        self._reset_parameters()

    def forward(self, seq_tokens,esm_rep,kla_site):
        if self.addseq:
            seq_embeding = self.embedding(seq_tokens)
            esm_embeding = self.esm_linear(esm_rep)
            out = torch.cat((seq_embeding,esm_embeding),dim=-1) 
        else:
            out = self.esm_linear(esm_rep)

        out = self.transformer_encoder(out) # out.shape=bs,seqlen,512 #kla_site.shape = bs,
        idx = kla_site[:,None,None].repeat(1,1,out.shape[-1])
        all_out = torch.gather(out,1,idx) # all_out.shape = bs,1,512 #取出了K位点的表示
        all_out = all_out.squeeze(1)
        # all_out = self.relu(self.fc1(all_out))
        # all_out = self.fc2(self.dropout((all_out)))
        all_out1 = self.linear1(all_out)
        all_out1 = self.linear2(all_out1)
        # all_out = all_out.squeeze(-1)
        if not self.multitask:
            return all_out1 # out.shape = bs,seqlen
        else:
            all_out2 = self.linear1_othertask(all_out)
            all_out2 = self.linear2_othertask(all_out2)
            return all_out1,all_out2
    def _reset_parameters(self):
        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)