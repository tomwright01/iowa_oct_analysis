import numpy as np
import datetime
import logging

logger = logging.getLogger(__name__)

def readBioptigenOct(fname):
    """Read a biotigen .oct file
    Params:
      fname - full path to .oct file
      
    Returns:
      {'headers':[],
       'frame_headers':[],
       'iamge_data':np.array(nBscans,scan_depth,nAscans),
       'doppler_data':np.array(nBscans,scan_depth,nAscans)}
       
    Examples
      headers, frame_headers, image_data, doppler_data = readBioptigenOct(fname)
      
    Based on the matlab function extractOctFunction.m
    by Brad Bower, Bioptigen Inc.
    """
    # create a dict to store the file header info
    # each item has the datatype and the length, 
    # if length is 0 then the value of data_length is used

    header_defs = {'FRAMEHEADER':[np.uint32,1],
                   'FRAMECOUNT':[np.uint32,1],
                   'LINECOUNT':[np.uint32,1],
                   'LINELENGTH':[np.uint32,1],
                   'SAMPLEFORMAT':[np.uint32,1],
                   'DESCRIPTION':[np.uint8,0],
                   'XMIN':[np.float64,1],
                   'XMAX':[np.float64,1],
                   'XCAPTION':[np.uint8,0],
                   'YMIN':[np.float64,1],
                   'YMAX':[np.float64,1],
                   'YCAPTION':[np.uint8,0],
                   'SCANTYPE':[np.uint32,1],
                   'SCANDEPTH':[np.float64,1],
                   'SCANLENGTH':[np.float64,1],
                   'AZSCANLENGTH':[np.float64,1],
                   'ELSCANLENGTH':[np.float64,1],
                   'OBJECTDISTANCE':[np.float64,1],
                   'SCANANGLE':[np.float64,1],
                   'SCANS':[np.uint32,1],
                   'FRAMES':[np.uint32,1],
                   'FRAMESPERVOLUME':[np.uint32,1],
                   'DOPPLERFLAG':[np.uint32,1],
                   'CONFIG':[np.uint8,0]}
    headers = {}
    
    # create a dict to store the file header info
    # each item has the datatype and the length, 
    # if length is 0 then the value of data_length is used
    frame_defs = {'FRAMEDATA':[np.uint32,1],
                  'FRAMEDATETIME':[np.uint16,0],
                  'FRAMETIMESTAMP':[np.float64,1],
                  'FRAMELINES':[np.uint32,1],
                  'FRAMESAMPLES':[np.uint16,0],
                  'DOPPLERSAMPLES':[np.uint16,0]}
    frame_headers = [] #copy this for each frame
    
    with open(fname, 'rb') as fid:
        magic_number = np.fromfile(fid, np.uint16, 2)
        version_number = np.fromfile(fid, np.uint16, 1)
        key_length = np.fromfile(fid, np.uint32, 1)
        key = ''.join([chr(v) for v in np.fromfile(fid, np.uint8, key_length)])
        data_length = np.fromfile(fid, np.uint32, 1)
        
        if key != 'FRAMEHEADER':
            logger.error('Failed loading frame header for file %s'.format(fname))
            raise RuntimeError('Error loading frame header')
        
        header_flag = False # set to True when all header keys are read
        while not header_flag:
            #note need to correct for the final read of these 3 lines after all headers are read.
            key_length = np.fromfile(fid, np.uint32, 1)
            key = ''.join([chr(v) for v in np.fromfile(fid, np.uint8, key_length)])
            data_length = np.fromfile(fid, np.uint32, 1)
            
            if key in header_defs.keys():
                if header_defs[key][1]:
                    headers[key] = np.fromfile(fid, header_defs[key][0], 1)
                else:
                    headers[key] = np.fromfile(fid, header_defs[key][0], data_length)
            else:
                header_flag = True
                logger.info('All headers read.')
                
        # Test for unsupported scan types
        if headers['SCANTYPE'] == 6:
            logger.warning('Mixed Density ("S") Scans Not Supported')
            raise RuntimeError('Mixed Density ("S") Scans Not Supported')
        
        # Capture File Header Length.
        fid.seek(-4,1) # correct for 4 preceeding 4 byte key length read
        file_header_length = fid.tell()
        logger.info('File Header Length=%d', file_header_length)
        
        # Initalise structure to store Frame Data
        logger.info('Initialising frame data structures')
        frame_headers = [] #list to store the frame header details
        image_data = np.empty((headers['FRAMECOUNT'], headers['LINELENGTH'], headers['LINECOUNT']),
                          dtype=np.uint16)
        
        if headers['DOPPLERFLAG']:
            doppler_data = np.copy(frames)
        else:
            doppler_data = None
        
        
        #Start reading frame data
        logger.info('Starting Reading frame data')
        current_frame_idx = 0 # Counter for frame position

        while current_frame_idx < headers['FRAMECOUNT']:
            key_length = np.fromfile(fid, np.uint32, 1)
            key = ''.join([chr(v) for v in np.fromfile(fid, np.uint8, key_length)])
            data_length = np.fromfile(fid, np.uint32, 1)
            
            frame_flag = False # set to True when the current frame is read
        
            if key == 'FRAMEDATA':
                logger.info('Starting frame:%i', current_frame_idx);
                frame_headers.append({})
        
                while not frame_flag:
                    key_length = np.fromfile(fid, np.uint32, 1)
                    key = ''.join([chr(v) for v in np.fromfile(fid, np.uint8, key_length)])
                    data_length = np.fromfile(fid, np.uint32, 1)
                    
                    #Deal with a couple of special cases
                    if key == 'FRAMEDATETIME':
                        frame_date_time = np.fromfile(fid, np.uint16, data_length/2)
                        date_time = datetime.datetime(frame_date_time[0],
                                                      frame_date_time[1],
                                                      frame_date_time[3],
                                                      frame_date_time[4],
                                                      frame_date_time[5],
                                                      frame_date_time[6],
                                                      frame_date_time[7])
                        frame_headers[current_frame_idx][key]=date_time
                    elif key == 'FRAMESAMPLES':
                        try:
                            # Framelines is not always defined on a per frame basis
                            # missing in some earlier versions
                            nLines = frame_headers[current_frame_idx]['FRAMELINES']
                        except KeyError:
                            nLines = headers['LINECOUNT']
                        image = np.fromfile(fid, np.uint16, nLines * headers['LINELENGTH'])
                        image_data[current_frame_idx,:,:] = image.reshape((headers['LINELENGTH'],
                                                                           nLines))
                    elif key =='DOPPLERSAMPLES':
                        try:
                            # Framelines is not always defined on a per frame basis
                            # missing in some earlier versions
                            nLines = frame_headers[current_frame_idx]['FRAMELINES']
                        except KeyError:
                            nLines = headers['LINECOUNT']
                        image = np.fromfile(fid, np.uint16, nLines * headers['LINELENGTH'])
                        doppler_data[current_frame_idx,:,:] = image.reshape((headers['LINELENGTH'],
                                                                             nLines))
                    elif key in frame_defs.keys():
                        frame_headers[current_frame_idx][key] = \
                            np.fromfile(fid, frame_defs[key][0],frame_defs[key][1])
                    else:
                        # frame complete
                        frame_flag = True
                        fid.seek(-4,1) # correct for keylength
                        current_frame_idx = current_frame_idx + 1
    return {'headers':headers,
            'frame_headers':frame_headers,
            'image_data':image_data,
            'doppler_data':doppler_data}
                    
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    fname = '../Data/Sample/Bioptigen/sample_OD_V_12x12_0_0000036.OCT'
    result = readBioptigenOct(fname)
    