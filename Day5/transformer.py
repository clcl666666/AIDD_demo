# Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
"""
antigen: AF2, aa sequence 
antibody: pssm, aa sequence.
interact
embedding info is not added in this version
"""
from typing import Optional, List

import torch
import torch.nn.functional as F
from torch import nn, Tensor
import numpy as np
import math
import copy
import time
from torch.autograd import Variable
import collections

class Transformer(nn.Module):

    def __init__(self, d_model=128, nhead=8, num_encoder_layers=1,
                 num_decoder_layers=1, dropout=0.1,
                 activation="relu",antigen_in_dimension = 405,antibody_in=1280): #405 # 抗体这的网络要改。
        super().__init__()
        self.seqprocess = nn.Linear(antigen_in_dimension, d_model)
        self.seqprocess2 = nn.Linear(antibody_in, d_model)
        self.copy_ = copy.deepcopy

        self.position = PositionalEncoding(d_model)

        decoder_layer = TransformerDecoderLayer(d_model, nhead, dim_feedforward=64)
        self.decoder = TransformerDecoder(decoder_layer, num_decoder_layers)

        decoder_layer2 = TransformerDecoderLayer2(d_model, nhead, dim_feedforward=64)
        self.decoder2 = TransformerDecoder2(decoder_layer2, num_decoder_layers)

        encoder_layer = TransformerEncoderLayer(d_model,nhead,dim_feedforward=64)
        self.encoder = TransformerEncoder(encoder_layer, num_encoder_layers)
        
        self.d_model = d_model
        self.nhead = nhead

        self.aggragator = AggregateLayer(d_model = 2*d_model)
        self.predictor = GlobalPredictor(d_model = 2*d_model,d_h=128, d_out=10)
        self._reset_parameters()

    def _reset_parameters(self):
        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)

    def forward(self, antigen_seq,antibody_seq1,antibody_seq2):
        # each seq has dimension 40
        
        antigen_seq = self.seqprocess(antigen_seq)  # msa_out.shape = [1, seqlen, 256]
        antibody_seq1 = self.seqprocess2(antibody_seq1)  # msa_out.shape = [1, seqlen, 256]
        antibody_seq2 = self.seqprocess2(antibody_seq2)  # msa_out.shape = [1, seqlen, 256]
        # antigen is query, antibody is key, value.
        antibody_seq1 = self.encoder(antibody_seq1)
        antibody_seq2 = self.encoder(antibody_seq2)

        antibody_seq1 = self.decoder2(antigen_seq,antibody_seq1)
        antibody_seq2 = self.decoder2(antigen_seq,antibody_seq2)

        out1 = self.decoder(antibody_seq1,antigen_seq)
        out2 = self.decoder(antibody_seq2,antigen_seq) # 1,antibody_seqlen,256
        
        out = torch.cat([out1,out2],dim = -1) #  1,antibody_seqlen,2*d_model
        out = self.aggragator(out)
        output = self.predictor(out)

        return output  
class GlobalPredictor(nn.Module):
    def __init__(self, d_model=None, d_h=None, d_out=None, dropout=0.5):
        super(GlobalPredictor, self).__init__()
        self.predict_layer = nn.Sequential(collections.OrderedDict([
            # ('batchnorm', nn.BatchNorm1d(d_model)),
            ('fc1', nn.Linear(d_model, d_h)),
            ('tanh', nn.Tanh()),
            ('dropout', nn.Dropout(dropout)),
            ('fc2', nn.Linear(d_h, d_out))
        ]))

    def forward(self, x):
        x = self.predict_layer(x)
        return x
class AggregateLayer(nn.Module):
    def __init__(self, d_model=None, dropout=0.1):
        super(AggregateLayer, self).__init__()        
        self.attn = nn.Sequential(collections.OrderedDict([
            ('layernorm', nn.LayerNorm(d_model)),
            ('fc', nn.Linear(d_model, 1, bias=False)),
            ('dropout', nn.Dropout(dropout)),
            ('softmax', nn.Softmax(dim=1))
        ]))

    def forward(self, context):
        '''
        Parameters
        ----------
        context: token embedding from encoder (Transformer/LSTM)
                (batch_size, seq_len, embed_dim)
        '''
        weight = self.attn(context)
        # (batch_size, seq_len, embed_dim).T * (batch_size, seq_len, 1) *  ->
        # (batch_size, embed_dim, 1)
        output = torch.bmm(context.transpose(1, 2), weight)
        output = output.squeeze(2)
        return output
class Aggregater(nn.Module):
    def __init__(self,d_model,h):
        super().__init__()
        self.d_model = d_model
        self.attn = nn.Linear(d_model,1)

        self.linear1 = nn.Linear(d_model,h)
        self.linear2 = nn.Linear(h,2)

        self.dropout = nn.Dropout(0.1)
        self.activate = nn.LeakyReLU()
        self.bn = nn.BatchNorm1d(d_model)
    def forward(self,seq):
        # seq.shape = bs,seqlen,d_model
        attn = self.attn(seq)
        # attn.shape = bs,seqlen,1
        seq = torch.bmm(seq.transpose(-1,-2),attn).squeeze(-1)
        # seq.shape = bs,d_model
        out = self.linear2(self.dropout(self.activate(self.linear1(self.bn(seq)))))
        return out
class TransformerDecoder(nn.Module):
    def __init__(self, decoder_layer, num_layers, norm=None):
        super().__init__()
        self.layers = _get_clones(decoder_layer, num_layers)
        self.num_layers = num_layers
    def forward(self, ab_seq, ag_seq):

        for layer in self.layers:
            ag_seq,attn = layer(ab_seq, ag_seq)
        return ag_seq

class TransformerDecoder2(nn.Module):
    def __init__(self, decoder_layer, num_layers, norm=None):
        super().__init__()
        self.layers = _get_clones(decoder_layer, num_layers)
        self.num_layers = num_layers
    def forward(self, ab_seq, ag_seq):

        for layer in self.layers:
            ag_seq,attn = layer(ab_seq, ag_seq)
        return ag_seq

class TransformerEncoder(nn.Module):
    def __init__(self, encoder_layer, num_layers, norm=None):
        super().__init__()
        self.layers = _get_clones(encoder_layer, num_layers)
        self.num_layers = num_layers

    def forward(self, ab_seq):

        for layer in self.layers:
            ab_seq = layer(ab_seq)

        return ab_seq
class MyTransformerEncoderLayer(nn.Module):
    def __init__(self,d_model,dim_feedforword,dropout=0.1,nhead = 8):
        super().__init__()
        self.linear1 = nn.Linear(d_model,dim_feedforword)
        self.linear2 = nn.Linear(dim_feedforword,d_model)
        self.dropout = nn.Dropout(dropout)
        self.multiheadAttn = MyMultiheadAttention(d_model,nhead)
        self.norm = nn.LayerNorm(d_model)
        self.activate = nn.LeakyReLU()
    def forward(self,seq):
        # seq.shape = bs,seqlen,d_model
        attn,out = self.multiheadAttn(seq,seq,seq)
        out = seq + self.dropout(out)
        out = self.norm(out)

        out = out+self.dropout(self.linear2(self.activate(self.linear1(out))))
        return out




class MyMultiheadAttention(nn.Module):
    def __init__(self,d_model=512,nhead=8):
        super().__init__()
        assert d_model%nhead == 0
        self.d_model = d_model
        self.d_model2 = d_model//nhead
        self.nhead = nhead
        self.linears = nn.ModuleList([nn.Linear(d_model,d_model) for i in range(3)])
        self.softmax = nn.Softmax(dim=3)
    def forword(self, query, key, value,mask=None): 
        # query.shape = bs,seqlen,512 1个head
        query,key,value = [self.linears[i][[query,key,value][i]] for i in range(3)]

        query = torch.stack(torch.split(query,split_size_or_sections=self.d_model2,dim=2),dim=0)
        key = torch.stack(torch.split(key,split_size_or_sections=self.d_model2,dim=2),dim=0)
        value = torch.stack(torch.split(value,split_size_or_sections=self.d_model2,dim=2),dim=0)
        # query.shape = 8,bs,seqlen,64
        attn = torch.matmul(query,key.transpose(-1,-2)) #attn.shape = h, bs,seqlen,seqlen
        if mask is not None:
            # 0表示mask, 1表示不mask mask.shape = bs,seqlen.
            bs,seqlen = mask.shape
            mask = mask[None,:,None,:].repeat(self.nhead,1,seqlen,1)
            attn = torch.masked_fill(attn,mask==0,value=-np.inf)

        attn = self.softmax(attn)/math.sqrt(self.d_model)
        out = torch.matmul(attn,value) # out.shape= h,bs,seqlen,d_model2
        out = torch.cat(torch.split(out,split_size_or_sections=1,dim=0),dim=-1).squeeze(0)
        return attn, out

class TransformerEncoderLayer(nn.Module):
## later add a sigmoid layer.
    def __init__(self, d_model, nhead, dim_feedforward, dropout=0.1,
                 activation="relu"):
        super().__init__()
        self.pos = PositionalEncoding(d_model)
        self.attn2 = nn.MultiheadAttention(
            d_model, nhead, dropout=dropout,bias = False,batch_first=True) # test later if use true of false.
        # Implementation of Feedforward model
        self.linear1 = nn.Linear(d_model, dim_feedforward)
        self.linear2 = nn.Linear(dim_feedforward, d_model)

        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        self.dropout1 = nn.Dropout(dropout)
        self.dropout2 = nn.Dropout(dropout)

        self.activation = _get_activation_fn(activation)

    def forward(self, ab_seq):
        ### query is antibody
        k = self.pos(ab_seq)
        seq_out,attention = self.attn2(k, k, value=k)
        src = k + self.dropout1(seq_out)
        src = self.norm1(src)

        src2 = self.linear2(self.dropout1(self.activation(self.linear1(src))))
        out = src + self.dropout2(src2)
        out = self.norm2(out)
        return out

class TransformerDecoderLayer(nn.Module):
## later add a sigmoid layer.
    def __init__(self, d_model, nhead, dim_feedforward, dropout=0.1,
                 activation="relu"):
        super().__init__()
        self.pos = PositionalEncoding(d_model)


        self.attn2 = nn.MultiheadAttention(
            d_model, nhead, dropout=dropout,bias = False,batch_first=True) # test later if use true of false.
        # Implementation of Feedforward model
        self.linear1 = nn.Linear(d_model, dim_feedforward)
        self.linear2 = nn.Linear(dim_feedforward, d_model)

        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        self.dropout1 = nn.Dropout(dropout)
        self.dropout2 = nn.Dropout(dropout)

        self.activation = _get_activation_fn(activation)

    def forward(self, ab_seq, ag_seq):
        ### query is antibody
        k = self.pos(ab_seq)
        q = self.pos(ag_seq)

        seq_out,attention = self.attn2(q, k, value=k)
        src = q + self.dropout1(seq_out)
        src = self.norm1(src)

        src2 = self.linear2(self.dropout1(self.activation(self.linear1(src))))
        out = src + self.dropout2(src2)
        out = self.norm2(out)
        return out,attention


class TransformerDecoderLayer2(nn.Module):
## later add a sigmoid layer.
    def __init__(self, d_model, nhead, dim_feedforward, dropout=0.1,
                 activation="relu"):
        super().__init__()
        self.pos = PositionalEncoding(d_model)

        self.attn2 = nn.MultiheadAttention(
            d_model, nhead, dropout=dropout,bias = False,batch_first=True) # test later if use true of false.
        # Implementation of Feedforward model
        self.linear1 = nn.Linear(d_model, dim_feedforward)
        self.linear2 = nn.Linear(dim_feedforward, d_model)

        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        self.dropout1 = nn.Dropout(dropout)
        self.dropout2 = nn.Dropout(dropout)

        self.activation = _get_activation_fn(activation)

    def forward(self, ab_seq, ag_seq):
        ### query is antibody
        k = self.pos(ab_seq)
        q = self.pos(ag_seq)
        
        seq_out,attention = self.attn2(q, k, value=k)
        src = q + self.dropout1(seq_out)
        src = self.norm1(src)

        src2 = self.linear2(self.dropout1(self.activation(self.linear1(src))))
        out = src + self.dropout2(src2)
        out = self.norm2(out)
        return out,attention

def clones(module, N):
    "Produce N identical layers."
    return nn.ModuleList([copy.deepcopy(module) for _ in range(N)])

def _get_clones(module, N):
    return nn.ModuleList([copy.deepcopy(module) for i in range(N)])


def build_transformer(args):
    return Transformer(
    )

class PositionalEncoding(nn.Module):
    "Implement the PE function."

    def __init__(self, d_model, max_len=2000):
        super(PositionalEncoding, self).__init__()
        # Compute the positional encodings once in log space.
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len).unsqueeze(1).float()
        div_term = torch.exp(torch.arange(
            0, d_model, 2).float() * -(math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        # pe.shape = torch.Size([1, 5000, d_model])
        self.register_buffer('pe', pe)

    def forward(self, x):
        x = x + Variable(self.pe[:, :x.size(1)],
                         requires_grad=False)
        return x

def _get_activation_fn(activation):
    """Return an activation function given a string"""
    if activation == "relu":
        return F.relu
    if activation == "gelu":
        return F.gelu
    if activation == "glu":
        return F.glu
    raise RuntimeError(F"activation should be relu/gelu, not {activation}.")

if __name__ == "__main__":
    antigen_seq = torch.rand(2,10,405)
    antibody_seq1 = torch.rand(2,20,1280)
    antibody_seq2 = torch.rand(2,21,1280)
    model = Transformer()
    out = model(antigen_seq,antibody_seq1,antibody_seq2)