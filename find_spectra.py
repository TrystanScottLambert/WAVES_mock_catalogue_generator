"""
Script to find spectra for all the galaxies in a given mock catalog.
"""

import numpy as np
import pandas as pd
import h5py
from read import read_all_spectra
from write import write_spectra_table_to_parquet


def get_ids(catalog_name: str) -> np.ndarray:
    """
    read in the sky_galaxy_ids from the given array.
    """
    df = pd.read_parquet(catalog_name, engine="pyarrow")
    return np.array(df["id_galaxy_sky"])


def main(mock_cat: str, directory: str, filestub: str, outfile: str) -> None:
    """
    Main script for scoping.
    """
    # quickly get wavelenghts
    with h5py.File(f"{directory}/{filestub}_00.hdf5") as f:
        wavelengths = np.array(f["wavegrid"][:])
    galaxy_ids = get_ids(mock_cat)
    spectra = read_all_spectra(directory, filestub, galaxy_ids)
    write_spectra_table_to_parquet(spectra, wavelengths, outfile)


if __name__ == "__main__":
    pass
