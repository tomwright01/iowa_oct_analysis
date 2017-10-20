"""
Produce OCT A-scan figures
"""
import re
import os
import pathlib
import logging
from App import OctLayers

logger = logging.getLogger(__name__)

def findFiles(basepath, regex):
    """
    Find files matching regex expression
    For each match returns a tuplet:
        (full_path_to_file, match_object)
    The match object contains the regex matches
    """
    pattern = re.compile(regex)
    file_list = []
    for root, dirs, files in os.walk(basepath):
        for filename in files:
            match = pattern.match(filename)
            if match:
                file_list.append((os.path.join(root, filename), match))
    logger.info('Found:{} files to process.'.format(len(file_list)))
    return(file_list)

def processFile(matchobject):
    oct_file = matchobject[0]
    surface_file = os.path.join(os.path.dirname(matchobject[0]), matchobject[1].group(1), matchobject[1].group(1) + '_Surfaces_Retina.xml')
    oct = OctLayers.OctLayers(filename = surface_file,
                              raw_filename =  oct_file)



if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    basepath = os.path.join(pathlib.Path.home(),'Documents/Sickkids/Stephanie Kletke/Working/Cirrus/')
    file_pattern = '((P\d{7})_Macular Cube .*_(OS|OD)_(sn\d{5})_cube_z).img'
    file_list = findFiles(basepath, file_pattern)
    processFile(file_list[0])
