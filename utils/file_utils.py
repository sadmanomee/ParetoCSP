import os
import pickle
import pandas as pd


def get_program_path():
    program_path = os.path.dirname(os.path.split(__file__)[0])
    return program_path


def check_path(path):
    if not os.path.exists(path):
        os.makedirs(path)


def check_and_rename_path(path):
    path = path.rstrip('/').rstrip('\\')
    path_tmp = path
    i = 1
    while os.path.exists(path_tmp):
        path_tmp = '%s_%d' % (path, i)
        i += 1
    if i > 1:
        os.rename(path, path_tmp)
    os.makedirs(path)


def save_data_csv(path, filename, data=None, header=None):
    if path is None:
        path = os.path.split(filename)[0]
        check_path(path)
        file_path = filename
    else:
        check_path(path)
        file_path = os.path.join(path, filename)
    df = pd.DataFrame(data, index=None, columns=header)
    if header is None:
        df.to_csv(file_path, index=False, header=False)
    else:
        df.to_csv(file_path, index=False)


def read_data_csv(path, filename, header=None, usecols=None, dtype=None):
    if path is None:
        file_path = filename
    else:
        file_path = os.path.join(path, filename)
    return pd.read_csv(file_path, delimiter=',', header=header, index_col=None, usecols=usecols, dtype=dtype)


def save_data_bin(path, filename, data=None):
    if path is None:
        path = os.path.split(filename)[0]
        check_path(path)
        file_path = filename
    else:
        check_path(path)
        file_path = os.path.join(path, filename)
    with open(file_path, 'wb') as f:
        pickle.dump(data, f, protocol=4)


def read_data_bin(path, filename):
    if path is None:
        file_path = filename
    else:
        file_path = os.path.join(path, filename)
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
        # data = pickle.load(f, encoding='iso-8859-1')
    return data
