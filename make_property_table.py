"""
Script to create the .md table with all the Galaxy/Group properties.
"""

from typing import TextIO
from property_dictionaries import GALAXY_PROPERTIES, GROUP_PROPERTIES


def write_property_table(file: TextIO, property_dict: dict) -> None:
    """
    Adds a .md table of the given property dictionary to the given file object.
    """
    file.write('| Property | Description | \n')
    file.write('| --- | --- | \n')
    for key, value in property_dict.items():
        file.write(f'| `{key}` | {value} | \n')

def main():
    """
    Main for scope.
    """
    with open('properties.md', 'w', encoding='utf-8') as file:
        file.write('# Properties \n')
        file.write('## Galaxy Properties \n')
        write_property_table(file, GALAXY_PROPERTIES)
        file.write('## Group Properties \n')
        write_property_table(file, GROUP_PROPERTIES)

if __name__ == '__main__':
    main()
