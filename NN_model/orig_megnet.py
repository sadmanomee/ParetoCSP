import numpy as np
from megnet.data.crystal import CrystalGraph
from megnet.data.graph import GaussianDistance
from megnet.models import MEGNetModel, GraphModel


class OrigMEGNet(MEGNetModel):
    def __init__(self):
        r_cutoff = 5
        nfeat_bond = 100
        gaussian_width = 0.5
        gaussian_centers = np.linspace(0, r_cutoff + 1, nfeat_bond)
        graph_converter = CrystalGraph(cutoff=r_cutoff)

        super().__init__(graph_converter=graph_converter,
                         centers=gaussian_centers,
                         width=gaussian_width,
                         # lr=1e-4,
                         metrics=['mae', 'mse'], )

    @classmethod
    def from_file(cls, model_path, graph_converter=None):
        from tensorflow.keras.models import load_model
        from megnet.layers import _CUSTOM_OBJECTS

        if graph_converter is None:
            r_cutoff = 5
            nfeat_bond = 100
            gaussian_width = 0.5
            gaussian_centers = np.linspace(0, r_cutoff + 1, nfeat_bond)
            graph_converter = CrystalGraph(cutoff=r_cutoff)

        model = load_model(model_path, custom_objects=_CUSTOM_OBJECTS)
        return GraphModel(model=model,
                          graph_converter=graph_converter,
                          centers=gaussian_centers,
                          width=gaussian_width,)


def train():
    from data.matbench_dataset import MatBenchDataset

    dataset = MatBenchDataset(rand_seed=100)
    model_dir = 'OrigMEGNet_model'

    gn_model = OrigMEGNet()
    gn_model.train(train_structures=dataset.training_set[0], train_targets=dataset.training_set[1],
                   validation_structures=dataset.val_set[0], validation_targets=dataset.val_set[1],
                   dirname=model_dir,
                   # prev_model=model_dir + '/val_mae_00144_0.037274.hdf5',
                   automatic_correction=True,
                   patience=50,
                   batch_size=200)


def predict():
    import os

    model_path = 'OrigMEGNet_model'
    if os.path.isdir(model_path):
        model_files = os.listdir(model_path)
        model_files = [mf for mf in model_files if '.hdf5' in mf]
        model_files.sort(key=lambda fn: os.path.getmtime(model_path + '/' + fn))
        model_file_path = os.path.join(model_path, model_files[-1])
    elif os.path.isfile(model_path):
        model_file_path = model_path
    else:
        print('Input a GN model path!')
        return
    gn_model = OrigMEGNet.from_file(model_file_path)

    from data.matbench_dataset import MatBenchDataset
    from tqdm import trange

    dataset = MatBenchDataset(rand_seed=100)
    test_structures, test_targets = dataset.test_set
    # pred_targets = model.predict_structures(test_targets).reshape(-1)
    ae = []
    for i in trange(len(test_targets)):
        try:
            # pred_target = pred_targets[i]
            pred_target = gn_model.predict_structure(test_structures[i]).reshape(-1)[0]
            true_target = test_targets[i]
            error = pred_target - true_target
            ae.append(abs(error))
        except:
            # print('bad data, skipped!')
            continue

    mae = sum(ae) / len(ae)
    print('MAE:', mae)


if __name__ == '__main__':
    train()
    predict()
