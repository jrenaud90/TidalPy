import os
from typing import List, Union


def get_all_files_of_type(directory_to_search: str, file_extensions: Union[List[str], str]):
    """ Returns all files with a specified extension(s) in a given directory

    :param directory_to_search:
    :param file_extensions:
    :return:
    """

    if type(file_extensions) == str:
        file_extensions = list(file_extensions)
    # splitext retains the first period, the user may not have expected that
    for e_i, extension in enumerate(file_extensions):
        if extension[0] != '.':
            file_extensions[e_i] = f'.{extension}'

    files = dict()
    # Only pull out files, not sub-directories
    only_files = [f for f in os.listdir(directory_to_search)
                  if os.path.isfile(os.path.join(directory_to_search, f))]

    for file_with_extension in only_files:
        filename, extension = os.path.splitext(file_with_extension)
        if extension in file_extensions:
            files[filename] = os.path.join(directory_to_search, file_with_extension)

    return files
