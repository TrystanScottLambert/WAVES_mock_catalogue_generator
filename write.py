"""
Module consisting of classes and functions which are used in the writing of the data.
"""

from typing import TextIO
from dataclasses import dataclass
import numpy as np

HEADER_PREAMBLE = 'header_preamble/header_preamble_v1.head'

def stamp_preamble(file: TextIO) -> None:
    """
    Writes the standard preamble text to an open file object. 

    Ex. with open('test.txt', 'w') as f:
            write_preamble(f)
            f.write('something \n')
    """
    with open(HEADER_PREAMBLE, encoding='UTF-8') as head:
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
        file.write('Lightcone Details: \n')
        file.write(f'area: {self.area} deg2 \n {self.mag_filter} <= {self.mag_cut} \n')
        file.write(f'redshift < {self.redshift_cut} \n VERSION: v{self.version} \n')

@dataclass
class DataProperty:
    """
    Class which insures names, units and data are stored togther correctly. 

    name is the name of the column (e.g. ra)
    unit_str is the unit string that will be included in the header of the file (e.g. [deg] J2000)
    data is the array that was generated. 
    """
    name: str
    unit_str: str
    data: np.ndarray

class TabularData:
    """
    Used to store the data into a tabular format to help with easy writing.
    """
    def __init__(self, data_properties: list[DataProperty], details: CatalogueDetails) -> None:
        """Reading in the data properties and catalogue details"""
        self.data_properties = data_properties
        self.details = details
