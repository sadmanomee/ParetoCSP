import random
import numpy as np
from abc import abstractmethod


class DatasetBase:
    def __init__(self,
                 n_train_rate: float = 0.5, n_val_rate: float = 0.125,
                 n_train: int = None, n_val: int = None,
                 is_shuffle: bool = True, rand_seed: int = None):
        structures, targets = self.get_dataset()
        self.structures, self.targets = np.array(structures), np.array(targets)
        data_count = len(targets)
        indexes = list(range(data_count))

        if is_shuffle:
            if rand_seed is not None:
                random.seed(rand_seed)
            random.shuffle(indexes)

        if n_train is None:
            self.n_train = int(data_count * n_train_rate)
        if n_val is None:
            self.n_val = int(data_count * n_val_rate)

    @abstractmethod
    def get_dataset(self):
        pass

    @property
    def training_set(self):
        return self.structures[:self.n_train], self.targets[:self.n_train]

    @property
    def val_set(self):
        return self.structures[self.n_train:self.n_train + self.n_val], self.targets[self.n_train:self.n_train + self.n_val]

    @property
    def test_set(self):
        return self.structures[self.n_train + self.n_val:], self.targets[self.n_train + self.n_val:]
