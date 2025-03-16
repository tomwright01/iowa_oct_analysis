import numpy as np
import re
import xml.etree.ElementTree
import logging

logger = logging.getLogger(__name__)

def readIowaSurfaces(fname):
    """Read the surfaces .xml file generated from OCTExplorer

    Params:
      fname - full path to the .xml file

    Returns:
    {'scan_size' - {'x','y','z'} - scan_size in voxels,
     'voxel_size' - {'x','y','z'} - voxel_size in microns,
     'scan_system' - 'cirrus'|'bioptigen' - Recording system manufacturer,
                         #TODO - find out what the value is for .jpg images and heidelburg
     'eye': 'OD'|'OS'
     'surface_names': [] - list of surfaces identified
     'surface_data': an [nSurfaces,nBscans,nAscans] numpy.maskedarray of the surface data
                     each (Surface,x,y) value is an integer value indicating the depth of the surface
                     pixels are counted from the top of the scan and numbered from 1

    Examples:

    Description:
    OCTExplorer.exe (https://www.iibi.uiowa.edu/content/iowa-reference-algorithms-human-and-murine-oct-retinal-layer-analysis-and-display)
    implements the Iowa Reference Algorithm to segment Ocular Coherence Tomography files from a variety of commercial systems.
    The principle output is an .xml file containing information delineating 10 retinal surfaces.
    This function reads the xml file to extract metadata and returns the surface information as numpy ndarray.
    """
    #define labels for the 10 retinal surfaces
    #surface_labels = ['ILM','RNFL-GCL','GCL-IPL','IPL-INL','INL-OPL',
                      #'OPL-HFL','BMEIS','IS/OSJ','IB_OPR','IB_RPE','OB_RPE']
    surface_labels = {}
    logger.debug('Loading surfaces file:{}'.format(fname))

    xml_root = xml.etree.ElementTree.parse(fname).getroot()

    # OCTexplorer version 3 used a <z> element for suface heights, version 4 uses <y>
    # going to remap everything back to version 3 format
    version = xml_root.find('version').text
    if version.startswith('3.0'):
        value_element = 'z'
        coord_element = 'y'
        
            # first lets extract the scan size information
        scan_size = {'x': int(xml_root.find('./scan_characteristics/size/x').text),
                     'y': int(xml_root.find('./scan_characteristics/size/y').text),
                     'z': int(xml_root.find('./scan_characteristics/size/z').text)}
        voxel_size = {'x': float(xml_root.find('./scan_characteristics/voxel_size/x').text),
                      'y': float(xml_root.find('./scan_characteristics/voxel_size/y').text),
                      'z': float(xml_root.find('./scan_characteristics/voxel_size/z').text)}


    else:
        value_element = 'y'
        coord_element = 'z'
        # first lets extract the scan size information
        scan_size = {'x': int(xml_root.find('./scan_characteristics/size/x').text),
                     'y': int(xml_root.find('./scan_characteristics/size/z').text),
                     'z': int(xml_root.find('./scan_characteristics/size/y').text)}
        voxel_size = {'x': float(xml_root.find('./scan_characteristics/voxel_size/x').text),
                      'y': float(xml_root.find('./scan_characteristics/voxel_size/z').text),
                      'z': float(xml_root.find('./scan_characteristics/voxel_size/y').text)}
        #voxel_size = {k:v*1000 for k,v in voxel_size.items()}

    # compatibility with version 3.8.0
    voxel_units = xml_root.find('./scan_characteristics/voxel_size/unit').text
    if voxel_units == 'mm':
        voxel_size['x'] = voxel_size['x'] * 1000
        voxel_size['y'] = voxel_size['y'] * 1000
        voxel_size['z'] = voxel_size['z'] * 1000

    # try to exract system information
    system = xml_root.find('./scan_characteristics/manufacturer').text.lower()
    if bool(re.search('carl zeiss',system)):
        system = 'cirrus'
    elif bool(re.search('Bioptigen',system)):
        system = 'bioptigen'
    else:
        logger.warn('Unknown system type')
        system = 'unknown'

    # structure to hold layer measurements
    # data in this structure is in pixels and can be used by the centering function
    nlayers = int(xml_root.find('surface_num').text)
    data = np.empty((nlayers,
                     scan_size['y'],
                     scan_size['x']))

    p = re.compile(r'.*\((.*)\)')
    for surface in xml_root.findall('surface'):
        # identify which surface this is and assign an index
        # can't use the label in the xml file as these are not contiguous
        surface_name = surface.find('name').text
        logger.debug('Loading surface:{}'.format(surface_name))
        surface_idx = np.NaN
        # extract the surface label
        match = re.match(p, surface_name)
        if match:
            if not match.group(1) in surface_labels.keys():
                #surface not seen before add the label and description
                surface_labels[match.group(1)] = (match.group(0),len(surface_labels))

            surface_idx = surface_labels[match.group(1)][1]
        else:
            logger.warning('Failed to identify surface:{}'.format(surface_name))
            break

        logger.debug('Surface index:{}'.format(surface_idx))
        # loop through all the bscans
        surface_bscans = surface.findall('bscan')
        for bscan_idx in range(data.shape[1]):
            bscan = surface_bscans[bscan_idx]
            data[surface_idx,bscan_idx,:] = [int(z.text) for z in bscan.findall(value_element)]

    # .xml file may also contain information on where segmentation has failed
    # create a structure to store this information
    undef_mask = np.zeros(data.shape,dtype=bool)
    undef_xml = xml_root.find('undefined_region')

    if undef_xml is not None:
        for ascan in undef_xml.findall('ascan'):
            x = int(ascan.find('x').text)
            y = int(ascan.find(coord_element).text)
            undef_mask[:,y,x] = True

    data = np.ma.MaskedArray(data, mask = undef_mask)

    laterality = xml_root.find('scan_characteristics').find('laterality').text
    if laterality.upper() in ['OD','OS']:
        laterality = laterality.upper()
    else:
        # not defined in the xml file, see if we can extract from the filename
        p = re.compile('(OD|OS)')
        m = re.search(p,fname)
        if m:
            laterality = m.group(0)

    return {'scan_size':scan_size,
            'voxel_size':voxel_size,
            'scan_system':system,
            'eye':laterality,
            'surface_names':surface_labels,
            'surface_data':data}

def readIowaCenter(fname):
    """Load the GridCenter.xml file
    Params:
      fname - full path to the _GridCenter_Iowa.xml file

    Returns:
      (center_x,center_y) - scan center in pixels
    """

    xml_root = xml.etree.ElementTree.parse(fname)
    
    version = xml_root.find('version').text
    if version.startswith('3'):
        value_element = 'z'
        coord_element = 'y'
    else:
        value_element = 'y'
        coord_element = 'z'


    c = xml_root.find('center')
    center_x = int(c.find('x').text)
    center_y = int(c.find(coord_element).text)

    return (center_x,center_y)
