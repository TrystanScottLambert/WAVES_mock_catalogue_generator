"""
Module consisting of classes and functions which are used in the writing of the data.
"""

from typing import TextIO
from dataclasses import dataclass

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


def write_catagloue(
    writable_dicts: list[dict],
    unit_header: str,
    cat_details: CatalogueDetails,
    outfile: str,
    delimeter: str = " ",
) -> None:
    """
    Co-adds all the writeable dicts in the order that they were given.
    Adds the headers too.
    """
    final_dict = {key: value for dictionary in writable_dicts for key, value in dictionary.items()}
    with open(outfile, "w", encoding="utf-8") as file:
        # Writing the header
        stamp_preamble(file)
        cat_details.stamp_details(file)
        file.write(unit_header)

        # Writing column names
        file.write(f'{delimeter.join(final_dict.keys())} \n')

        # Writing column data
        for row in zip(*final_dict.values()):
            file.write(f'{delimeter.join(map(str, row))} \n')
