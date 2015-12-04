import logging
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import re
import fnmatch
import App.OctLayers as Oct


def getLayerAsZScore(patient,controls,surface1,surface2):
    layer = patient.getThickness(surface1, surface2)
    layer_mean = controls.getLayerValues(surface1, surface2, 'mean')
    layer_sd = controls.getLayerValues(surface1, surface2, 'stdev')
    layer_zScore = (layer - layer_mean) / layer_sd
    return layer_zScore

def load_control_data(path, laterality):
    controls = Oct.OctCollection(folder = path,
                                 laterality = laterality)
    return controls
    
def plotLayerZscores(patient,controls,surface1,surface2,filename):
    colors = [(0,0,0),(1,0,0),(0,1,0)]
    cm = mpl.colors.ListedColormap(colors)    
    thickness = getLayerAsZScore(patient, controls, surface1, surface2)
    thickness_map = np.zeros(thickness.shape,dtype=np.int)
    thickness_map[thickness.data >= 2] = 1
    thickness_map[thickness.data < 2] = 2
    thickness_map[thickness.mask] = 0
    fig = plt.figure(figsize=(3,3),dpi=300)
    ax = fig.add_subplot(111)
    ax.imshow(thickness_map, aspect='auto',cmap=cm, interpolation='none')
    plt.savefig(filename)
    plt.close(fig)
    
if __name__=='__main__':
    logging.basicConfig(level=logging.DEBUG)
    controls = load_control_data(path = 'Data/Sample/',
                                 laterality='OD')

    p=re.compile('^P(\d{7}).*(\d{1,2}-\d{1,2}-\d{4}).*(OD|OS)_(sn\d+)')
    p2=re.compile('^(.+)_Surfaces_Iowa.xml')
    
    surface_files = []
    center_files = []    
    for dirname, _, filenames in os.walk('Data/Patient/'):
        for name in fnmatch.filter(filenames,'*Surfaces_Iowa.xml'):
            surface_files.append((dirname,name))
        for name in fnmatch.filter(filenames,'*GridCenter_Iowa.xml'):
            center_files.append((dirname,name))     
    
    for surface_file in surface_files:
        basename = re.match(p2,surface_file[1]).group(1)
        try:
            i = [i for i,v in enumerate(center_files) if re.search(basename,v[1])][0]
        except IndexError:
            i = None
            
        if i is not None:
            patient = Oct.OctLayers(filename = os.path.join(surface_file[0], surface_file[1]),
                                    center_filename = os.path.join(center_files[i][0],center_files[i][1]))
            patient.centerData()
        else:
            logger.warning('GridCenter file not found for surface file:{}'.format(surface_file[1]))
            patient = OctLayers(filename = os.path.join(surface_file[0], surface_file[1]))
        patient.setOrient('OD')
        
        file_info = re.match(p,surface_file[1])
        layer_set = [(0,10,'Total'),
                     (0,1,'RNFL'),
                     (0,3,'RNFL+GCL'),
                     (0,4,'Inner'),
                     (5,9,'Outer')]
        
        #open a text file for writing
        tex_basename = '{}-{}-{}-{}'.format(file_info.group(1),
                                            file_info.group(2),
                                            file_info.group(3),
                                            file_info.group(4))
        #create a folder to hold the figures
        tex_figure_folder = os.path.join('Output/',tex_basename)
        if not os.path.exists(tex_figure_folder ):
            os.mkdir(tex_figure_folder)
        
        #open the tex file for writing
        tex_filename = os.path.join('Output/',tex_basename + '.tex')
        tex_file = open(tex_filename, 'w')
        
        tex_file.write('\\documentclass{article}\n')
        tex_file.write('\\usepackage{graphicx}\n')
        tex_file.write('\\begin{document}\n')
        tex_file.write('\\title{{{}}}\n'.format(tex_basename))
        tex_file.write('\\maketitle\n')
        for layers in layer_set:
            # generate the figure
            fig_filename = os.path.join(tex_figure_folder,layers[2]+'.png')
            tex_fig_filename = os.path.join(tex_basename,layers[2]+'.png')
            plotLayerZscores(patient, controls, layers[0], layers[1], fig_filename)
            tex_file.write('\\begin{figure}\n')
            tex_file.write('\t\\centering\n')
            tex_file.write('\t\\includegraphics[width=3.0in]{{{}}}\n'.format(tex_fig_filename))
            tex_file.write('\t\\caption{{layer {} to {} - {}}}\n'.format(*layers))
            tex_file.write('\\end{figure}\n')
        tex_file.write('\\end{document}\n')
        tex_file.close()
    