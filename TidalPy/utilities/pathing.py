from typing import Union, List
import os

def get_all_files_of_type(directory_to_search: str, file_extensions: Union[List[str], str]):
    """ Returns all files with a specified extension(s) in a given directory

    :param directory_to_search:
    :param file_extensions:
    :return:
    """

    files = dict()
    only_files = [f for f in os.listdir(directory_to_search) if os.path.isfile(os.path.join(directory_to_search, f))]
    for filepath in only_files:
        filename, extension = os.path.splitext(filepath)
        if extension in file_extensions:
            files[filename] = filepath

    return files
