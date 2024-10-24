"""
Module consisting of classes and functions which are used in the writing of the data.
"""

from typing import TextIO
from dataclasses import dataclass

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

def write_units_in_header(file: TextIO, properties: list, property_dictionary: dict) -> None:
    """
    Writes the description of the properties in the header of the file using the property dictionary
    which can either be GALAXY_PROPERTIES or GROUP_PROPERTIES.
    """
    for _property in properties:
        file.write(f'# {_property}: {property_dictionary[_property]} \n')
