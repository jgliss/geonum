import glob
import os

import srtm


def get_local_srtm_files():
    """
    Get list of locally stored SRTM data files

    Returns
    -------
    list
        list of file paths
    """
    fh = srtm.utils.FileHandler()
    return glob.glob(f'{fh.local_cache_dir}/*.hgt')


def delete_local_srtm_files():
    """Delete all locally stored SRTM files"""
    for file in get_local_srtm_files():
        print(f'Deleting SRTM data file at {file}')
        os.remove(file)
