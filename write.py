"""
Module consisting of classes and functions which are used in the writing of the data.
"""

from typing import TextIO
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

HEADER_PREAMBLE_DIR = "header_preamble/header_preamble_v1.head"


def stamp_preamble(file: TextIO) -> None:
    """
    Writes the standard preamble text to an open file object.

    Ex. with open('test.txt', 'w') as f:
            write_preamble(f)
            f.write('something \n')
    """
    with open(HEADER_PREAMBLE_DIR, encoding="UTF-8") as head:
        header_lines = head.readlines()
        for line in header_lines:
            file.write(line)


@dataclass
class CatalogueDetails:
    """
    Stores the data of the details of the mock construction.

    area in square degrees.
    mag_filter is the filter that was chosen for the magnitude cut.
    """

    area: float
    mag_filter: str
    mag_cut: float
    redshift_cut: float
    version: float

    def stamp_details(self, file: TextIO):
        """
        Add the details to the file object when writing.
        """
        file.write("# Lightcone Details: \n")
        file.write(
            f"# area: {self.area} deg2 \n# {self.mag_filter} <= {self.mag_cut} \n"
        )
        file.write(f"# redshift < {self.redshift_cut} \n# VERSION: v{self.version} \n")

    def to_dict(self):
        """Returns a dictionary with all values as strings."""
        orig_dict = self.__dict__
        string_dict = {key: str(value) for key, value in orig_dict.items()}
        return string_dict


def write_to_parquet(
    writeable_dicts: list[dict],
    unit_header: dict,
    cat_details: CatalogueDetails,
    outfile: str,
) -> None:
    """
    Co-adds the dictionaries together and writes a parquet file with the column descriptions in the
    metadata as well as the catalogue details and other meta data.
    """
    # Combine the dictionaries
    final_dict = {
        key: value
        for dictionary in writeable_dicts
        for key, value in dictionary.items()
    }
    # Create the data frame:
    for key, value in final_dict.items():
        print(key, len(value))
    df = pd.DataFrame.from_dict(final_dict)

    # create the meta_data
    pa_list = [pa.field(key, pa.array(value).type) for key, value in final_dict.items()]
    meta_data = unit_header | cat_details.to_dict()
    pa_schema = pa.schema(pa_list, metadata=meta_data)

    # Add together and write to parquet
    table = pa.Table.from_pandas(df, pa_schema)
    pq.write_table(table, outfile)


def write_spectra_table_to_parquet(
    data: np.ndarray[np.ndarray], wavelength: np.ndarray, outfile: str
) -> None:
    """
    Writes a 2D numpy array that represents the spectra table to a parquet.
    The first column is assumed to be 'id_galaxy_sky', and the remaining columns are spectra values.
    """
    df = pd.DataFrame(data, columns=["id_galaxy_sky"] + wavelength.astype(str).tolist())
    df.to_parquet(outfile, index=False)
