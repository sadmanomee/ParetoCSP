import os
import re

from utils.file_utils import get_program_path, read_data_bin

elements_info_path = os.path.join(get_program_path(),
                                  'data/elements_info.pkl')
elements_info = read_data_bin(None, elements_info_path)


def compound_split(compound):
    """
        Split the compound into elements and corresponding count
        :param compound:
        :return:
    """
    temp_str = compound.replace(' ', '')
    pattern = re.compile(r'\d+')
    count = re.findall(pattern, temp_str)
    count = [int(x) for x in count]
    elements = re.split(pattern, temp_str)[:len(count)]

    return elements, count
