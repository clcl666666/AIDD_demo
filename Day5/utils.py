import sys
import torch
import os 
import pandas as pd
import math
from sklearn.utils import shuffle
import numpy as np
import math
import random

from typing import Sequence, Tuple, List

class Logger(object):
    """Writes both to file and terminal"""
    def __init__(self, savepath, mode='a'):
        self.terminal = sys.stdout
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        self.log = open(os.path.join(savepath, 'logfile.log'), mode)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.log.flush()


def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def randomSeed(random_seed):
    """Given a random seed, this will help reproduce results across runs"""
    if random_seed is not None:
        torch.manual_seed(random_seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(random_seed)

class EarlyStopping(object):
    def __init__(self, 
            patience=100, eval_freq=1, best_score=None, 
            delta=1e-9, higher_better=True):
        self.patience = patience
        self.eval_freq = eval_freq
        self.best_score = best_score
        self.delta = delta
        self.higher_better = higher_better
        self.counter = 0
        self.early_stop = False
    
    def not_improved(self, val_score):
        if np.isnan(val_score):
            return True
        if self.higher_better:
            return val_score < self.best_score + self.delta
        else:
            return val_score > self.best_score - self.delta
    
    def update(self, val_score):
        if self.best_score is None:
            self.best_score = val_score
            is_best = True
        elif self.not_improved(val_score):
            self.counter += self.eval_freq
            if (self.patience is not None) and (self.counter > self.patience):
                self.early_stop = True
            is_best = False
        else:
            self.best_score = val_score
            self.counter = 0
            is_best = True
        return is_best
def seed_everything(seed):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = True
