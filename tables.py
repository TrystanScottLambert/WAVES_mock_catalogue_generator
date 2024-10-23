"""
Functionality for building the different tables from the dictionaries and 
putting them into writable format. 
"""

from dataclasses import dataclass
from abc import ABC, abstractmethod

import numpy as np

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


class CalculatedProperty(ABC):
    """
    Base class for the calculated properties. Ensuring that units and methods are stored correctly.
    """
    @abstractmethod
    def method(self, *args, **kwargs) -> DataProperty:
        """
        The method that is implemented to convert the read in data into the calculated property.
        """
