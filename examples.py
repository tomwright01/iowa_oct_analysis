import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import App.OctLayers as oct

def layerMaskExample(obj):
    """Demonstrates how to get a mask for an OCT layer and apply it to the OCT images"""
    layer_mask = obj.getOctLayerMask(3,4,True)
    myoct = np.copy(obj.octdata)
    myoct[layer_mask] = 255
    plt.imshow(myoct[0,:,:])
    plt.show()
    
def etdrsThicknessExample(obj):
    """Demonstrates how to get layer thickness values for the 9 etdrs regions"""
    thickness = obj.getEtdrsThickness(3,4,True)
    print(thickness)

def etdrsIntensityExample(obj):
    """Demonstrates how to get layer intensity values for the 9 etdrs regions"""
    intensities = obj.getEtdrsIntensity(3,4,True)
    print(intensities)

if __name__=='__main__':
    plt.gray()
    
    fname = 'Data/Sample/Cirrus/IMGExportFiles/Sample1/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41329_cube_z.img'
    fname_layers = 'Data/Sample/Cirrus/IMGExportFiles/Sample1/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z_Surfaces_Iowa.xml'
    fname_centers = 'Data/Sample/Cirrus/IMGExportFiles/Sample1/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z_GridCenter_Iowa.xml'
    
    data_oct = oct.OctLayers(filename=fname_layers,
                             center_filename=fname_centers,
                             raw_filename=fname)
    data_oct.centerData()
    
    # layerMaskExample(data_oct)
    # etdrsThicknessExample(data_oct)
    etdrsIntensityExample(data_oct)