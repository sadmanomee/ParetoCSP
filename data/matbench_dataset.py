from data.dataset_base import DatasetBase


class MatBenchDataset(DatasetBase):
    def get_dataset(self):
        from matminer.datasets import load_dataset

        df = load_dataset("matbench_mp_e_form")
        structures, targets = df['structure'], df['e_form']

        return structures, targets


