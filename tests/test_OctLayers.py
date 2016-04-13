import pytest
import App.OctLayers as Oct

@pytest.fixture(scope="module",
                params=['cirrus','bioptigen'])
def oct(request):
    if request.param == 'cirrus':
        fname = 'Data/Sample/Cirrus/Sample1/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z.img'
        fname_layers = 'Data/Sample/Cirrus/Sample1/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z_Surfaces_Iowa.xml'
        fname_centers = 'Data/Sample/Cirrus/Sample1/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z_GridCenter_Iowa.xml'    
    elif request.param == 'bioptigen':
        fname = 'Data/Sample/Bioptigen/sample_OD_V_12x12_0_0000036.OCT'
        fname_layers = 'Data/Sample/Bioptigen/sample_OD_V_12x12_0_0000036/sample_OD_V_12x12_0_0000036_Surfaces_Iowa.xml'
        fname_centers = 'Data/Sample/Bioptigen/sample_OD_V_12x12_0_0000036/sample_OD_V_12x12_0_0000036_GridCenter_Iowa.xml'            

    return Oct.OctLayers(filename=fname_layers,
                        center_filename=fname_centers,
                        raw_filename=fname)

def test_load(oct):
    assert oct.data is not None
    if oct.system == 'cirrus':
        assert oct.center_x == 253
        assert oct.center_y == 66
    elif oct.system == 'bioptigen':
        assert oct.center_x == 500
        assert oct.center_y == 50
    else:
        print 'Invalid system'