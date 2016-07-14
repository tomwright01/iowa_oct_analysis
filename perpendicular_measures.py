import numpy as np
import logging

logger = logging.getLogger(__name__)

def find_thickness(surface2, start_pixel, vector, step=0.1, start_t=0):
    """Extends the 3D vector until it intersects with the surface
    Arguments:
    surface2 - an ndArray of values (the surface to be measured to)
    start_pixel [x,y,z] coordinates of the starting point of the vector (the pixel to be measured from)
    vector - [x,y,z] defining a 3d vector perpendicular to the starting surface
    step - step size to extend the vector, larger step size runs faster, smaller step size is more accurate
    start_t - an initial value of t (vector length). This allows the function to be initially run with a large step size, then a smaller one for more accurancy
    
    Returns:
    t - the largest vector multiple that doesn't cross the target surface
    
    Caution:
    Currently the routine can get stuck in a loop if the vector points in the wrong direction
    """
    
    found = False
    t = start_t
    
    while not found:
        # the main loop, vector is multiplied by factor t until it crosses the target surface
        t = t + step
        scaled_vec = t * vector
        new_pix = scaled_vec + start_pixel
        
        # First check the new x,y coordinates are valid. If not return np.nan
        if new_pix[0] > surface2.shape[1] - 1  or new_pix[1] > surface2.shape[0] - 1 \
           or new_pix[0] < 0 or new_pix[1] < 0:
            #import pdb; pdb.set_trace()
            return np.nan
       
        # Check if both vector components are 0
        if np.isclose(scaled_vec[0], 0) and np.isclose(scaled_vec[1], 0):
            # no interpolation required
            new_val = surface2[new_pix[1], new_pix[0]]

        # identify the surrounding pixels for interpolation
        low_x = np.floor(new_pix[0])
        high_x = np.ceil(new_pix[0])
        low_y = np.floor(new_pix[1])
        high_y = np.ceil(new_pix[1])

        # now check if the new pixels fall exactly on a real pixel (ie. new_pix is an integer)
        if not (np.isclose(new_pix[0], int(new_pix[0])) or np.isclose(new_pix[1], int(new_pix[1]))):
            # bilinear interpolation
            points = [(low_x,low_y,surface2[low_y,low_x]),
                      (low_x,high_y,surface2[low_y,high_x]),
                      (high_x,low_y,surface2[high_y,low_x]),
                      (high_x,high_y,surface2[high_y,high_x])]
            new_val = bilinear_interpolation(new_pix[0], new_pix[1], points)
            
        elif np.isclose(new_pix[0], int(new_pix[0])):
            # only x pixel exists, linear interp on y
            new_val = np.interp(new_pix[1], [low_y,high_y],[surface2[low_y,new_pix[0]],
                                                            surface2[high_y,new_pix[0]]])
        elif np.isclose(new_pix[1], int(new_pix[1])):
            # only y vector component is 0, linear interp on x
            new_val = np.interp(new_pix[0], [low_x,high_x],[surface2[new_pix[1],low_x],
                                                            surface2[new_pix[1],high_x]])
        else:
            # no interpolation required
            new_val = surface2[start_pixel[1], start_pixel[0]]
                        
            
        if new_pix[2] <= new_val:
            vec_length = np.linalg.norm(t * vector)
            return t - step
        
def bilinear_interpolation(x, y, points):
    '''Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,
        ...                        [(10, 4, 100),
        ...                         (20, 4, 200),
        ...                         (10, 6, 150),
        ...                         (20, 6, 300)])
        165.0

    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    points = sorted(points)               # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        pdb.pm()
        raise ValueError('points do not form a rectangle')
        
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        raise ValueError('(x, y) not within the rectangle')

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
           ) / ((x2 - x1) * (y2 - y1) + 0.0)

def generateNormalVectors(surface, size=1, voxel_size=[1,1,1]):
    """Generate normal vectors for a surface
    A vector is generated for each valid pixel,
    The normal vector is the vector that is perpendicular to the surface
    Arguments:
    surface - an ndArray, the source surface
    size - The size of sample to use for each voxel, 
        i.e. if size=2 the area -2:3 will be used to determine the surface slope
    voxel_size - [x,y,z] conversion factors from pixels to microns
    
    Returns:
    an ndArray of the same shape as surface where each point is an (x,y,z) vector
       Invalid points (ones where the sample area would be out of bounds)
       have the value np.nan
    """
    
    # start and end points for each vector across a point
    points = np.zeros((4,2),
                      dtype=[('x', 'i4'),('y', 'i4'),('z','f8')])
    points[:,0]['x'] = [i * calc_size for i in [-1, 0, 1, 1]]
    points[:,0]['y'] = [i * calc_size for i in [1, 1, 1, 0]]
    points[:,1]['x'] = [i * calc_size for i in [1, 0, -1, -1]]
    points[:,1]['y'] = [i * calc_size for i in [-1, -1, -1, 0]]
    
    # not going to calulate vectors for pixels on the edges of the surface
    valid_size = (my_surface.shape[0]-(calc_size * 2)) * (my_surface.shape[1]-(calc_size * 2))

    # the points array for each valid pixel
    points = np.repeat(np.expand_dims(points,2),valid_size,2)
    
    coords_x = range(calc_size,my_surface.shape[1]-calc_size)
    coords_y = range(calc_size,my_surface.shape[0]-calc_size)
    
    xv, yv = np.meshgrid(coords_x,coords_y, indexing='ij')
    
    points[:,0,:]['x'] = points[:,0,:]['x'] + xv.flatten()
    points[:,1,:]['x'] = points[:,1,:]['x'] + xv.flatten()
    points[:,0,:]['y'] = points[:,0,:]['y'] + yv.flatten()
    points[:,1,:]['y'] = points[:,1,:]['y'] + yv.flatten()
    
    # extract the z coordinates
    x = points['x']
    y = points['y']
    points['z'] = my_surface[y,x]  
    
    # points is now a 3d array [a, b, c]
    # each datapoint is a coordinate (x,y,z)
    # a = 4 - The number of vectors crossing each pixel
    # b = 2 - Point coordinates of the end of each vector (0=start, 1=end)
    # c = N number of valid pixels in the surface    
    
    # Calculate the vectors
    vectors = np.empty((4,points.shape[2]),
                       dtype=[('x','f8'),('y','f8'),('z','f8')])
    vectors['x'] = points[:,1,:]['x'] - points[:,0,:]['x']
    vectors['y'] = points[:,1,:]['y'] - points[:,0,:]['y']
    vectors['z'] = points[:,1,:]['z'] - points[:,0,:]['z']
    
    # Define the vector pairs
    pairs = np.array([(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)], dtype=np.int)
    
    # define some structures to hold the data
    normal_vector = np.empty(valid_size, dtype=(np.float,3))
    vector_pairs = np.empty((6,2), dtype=(np.float,3))
    normal_vectors = np.empty(my_surface.shape, dtype=(np.float,3))
    normal_vectors[:] = np.nan
    
    # convert voxels to microns
    # this is only rally needed if the x,y,z dimensions are not equal
    vectors["x"] = vectors["x"] * voxel_size[0]
    vectors["y"] = vectors["y"] * voxel_size[1]
    vectors["z"] = vectors["z"] * voxel_size[2]    
    
    old_pcomp = 0.0 # variable to store progress
    for i in range(valid_size):
        # long running process so calculate progress
        pcomp = (float(i) / valid_size)*100
        if int(pcomp) > old_pcomp:
            old_pcomp = pcomp
            logging.info('percent: {}'.format(int(pcomp)))
            
        # need to convert from structured arrays (np.void) to ndarray
        vector_pairs[:, 0, :] = np.matrix([vectors[pairs[:,0], i]['x'],
                                           vectors[pairs[:,0], i]['y'],
                                           vectors[pairs[:,0], i]['z']]).transpose()
        vector_pairs[:, 1, :] = np.matrix([vectors[pairs[:,1], i]['x'],
                                           vectors[pairs[:,1], i]['y'],
                                           vectors[pairs[:,1], i]['z']]).transpose()
        
        # convert the vector pairs from coordinates to microns
        #import pdb; pdb.set_trace()
        #vector_pairs[:, 0, :] = vector_pairs[:, 0, :] * voxel_sizes
        #vector_pairs[:, 1, :] = vector_pairs[:, 1, :] * voxel_sizes
        
        # calculate a normal vector for each vector pair
        potential_vectors = np.cross(vector_pairs[:, 0],
                                  vector_pairs[:, 1])
        
        #import pdb; pdb.set_trace()
        # average the 6 normal vectors to find the mean value for the pixel
        normal_vector = potential_vectors.mean(axis=0) 
        
        # Normalise the vector to unit length
        normal_vector = normal_vector / np.linalg.norm(normal_vector)
        
        # Convert the vectors back into pixel coordinates
        normal_vector = normal_vector / voxel_sizes
        
        # The normal vector could be in either direction
        # TODO: need to add a check for direction
        if normal_vector[2] > 0:
            normal_vector = 0 - normal_vector
        #normal_vector[1] = 0    
        # create an array to hold the output thickness values
        #thicknesses = np.empty((my_surface.shape[0]-2,my_surface.shape[1]-2),
                               #dtype=(np.float,3))
    
        # note output from unravel is [y,x]
        start_pixel = np.unravel_index(i, (my_surface.shape[0]-(2 * calc_size),my_surface.shape[1]-(2 * calc_size)))
        
        # append the origin z coordinate to the start pixel
        surface_pixel = [val + calc_size for val in start_pixel]
        surface_pixel = surface_pixel[::-1]
        surface_pixel.append(my_surface[surface_pixel[1], surface_pixel[0]])
        
        # Store the normal vector incase I want it later
        normal_vectors[surface_pixel[1], surface_pixel[0]] = normal_vector
        #import pdb;pdb.set_trace()
    