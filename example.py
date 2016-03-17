import logging
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import re
import datetime
import fnmatch
import dill
import App.OctLayers as Oct
from PIL import Image

def getLayerAsZScore(patient,controls,surface1,surface2):
    layer = patient.getThickness(surface1, surface2)
    layer_mean = controls.getLayerValues(surface1, surface2, 'mean')
    layer_sd = controls.getLayerValues(surface1, surface2, 'stdev')
    layer_zScore = (layer - layer_mean) / layer_sd
    layer_zScore.mask = layer.mask
    return layer_zScore

def getLayerRange(controls,surface1,surface2):
    layer_mean = controls.getLayerValues(surface1, surface2, 'mean')
    layer_sd = controls.getLayerValues(surface1, surface2, 'stdev')
    
    layer_vals = layer_mean + (4 * layer_sd)
    return (layer_vals.min(), layer_vals.max())

def plotControlLayer(controls,range,surface1,surface2,filename):
    layer = controls.getLayerValues(surface1, surface2, 'mean')
    layer[0,0] = range[0]
    layer[-1,0] = range[1]
    
    # scale values to be between 0 and 1
    layer = (layer - layer.min()) / layer.max()
    
    fig = plt.figure(figsize=(3,3),dpi=300)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=1, bottom=0, left=0, right=1)
    ax.imshow(layer,
              aspect='auto',
              cmap=plt.get_cmap('jet'),
              interpolation='none')
    plt.savefig(filename,dpi=300)
    plt.close(fig)    
    
def load_control_data(path, laterality, nmax=None):
    #pkl_filename = 'Data/controls_obj.pkl'
    
    #if os.path.isfile(pkl_filename):
    #    with open(pkl_filename, 'rb') as input:
    #        controls = dill.load(input)
    #else:
    controls = Oct.OctCollection(folder = path,
                                 laterality = laterality,
                                 nmax=nmax)
    # extract the patient id from the filename and append to the OCTLayer object
    p=re.compile('.*([0-9]{3})_[OD|OS]')
    for data in controls.data:
        data.serialNo = re.match(p, data.filename).group(1)
    #    with open(pkl_filename, 'wb') as output:
    #        dill.dump(controls, output)

    return controls

def plotLayer(patient, controls, surface1, surface2, filename):
    layerRange = getLayerRange(controls, surface1, surface2)
    layer = patient.getThickness(surface1, surface2)
    
    # going to force two pixels to range max and min
    layer[0,0] = layerRange[0]
    layer[-1,0] = layerRange[1]
    layerRange = (layer.min(), layer.max())
    # scale values to be between 0 and 1
    #layer = (layer - layer.min()) / layer.max()
    
    fig = plt.figure(figsize=(3,3),dpi=300)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=1, bottom=0, left=0, right=1)
    im = ax.imshow(layer,
                   aspect='auto',
                   cmap=plt.get_cmap('jet'),
                   interpolation='none')
    fig.colorbar(im)
    plt.savefig(filename,dpi=300)
    plt.close(fig)
    return layerRange
    
def plotScale(range, filename):
    gradient = np.linspace(0, 1, 255)
    gradient = np.vstack((gradient,gradient)).transpose()
    #gradient = np.flipud(gradient)
    
    scale = range[1] - range[0]
    gradient = gradient * scale
    gradient = gradient + range[0]
    
    fig = plt.figure(figsize=(0.75,3),dpi=300)
    ax = fig.add_subplot(111)
    ax.axes.get_xaxis().set_visible(False)
    fig.subplots_adjust(top=1, bottom=0, left=0.5, right=1)
    im = ax.imshow(gradient,
                   aspect='auto',
                   cmap=plt.get_cmap('jet'),
                   interpolation='none')
    fig.colorbar(im)
    plt.savefig(filename,dpi=300)
    plt.close(fig)    


def plotLayerZscores(patient,controls,surface1,surface2,filename):
    colors = [(0,0,0),(1,0,0),(0,1,0)]
    cm = mpl.colors.ListedColormap(colors)    
    thickness = getLayerAsZScore(patient, controls, surface1, surface2)
    thickness_map = np.zeros(thickness.shape,dtype=np.int)
    thickness_map[thickness.data >= 2] = 1
    thickness_map[thickness.data < 2] = 2
    thickness_map[thickness.mask] = 0
    thickness_map[0,0] = 0
    thickness_map[0,1] = 1
    thickness_map[0,2] = 2
    fig = plt.figure(figsize=(3,3),dpi=300)
    ax = fig.add_subplot(111)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    fig.subplots_adjust(top=1, bottom=0, left=00, right=1)
    ax.imshow(thickness_map, aspect='auto',cmap=cm, interpolation='none')
    plt.savefig(filename,dpi=300)
    plt.close(fig)
    
def assemble_images(basepath):
    im_scale = Image.open(basepath + '_scale.png')
    im_layer = Image.open(basepath + '_layer.png')
    im_inset = Image.open(basepath + '_inset.png')
    im_control = Image.open(basepath + '_control.png')

    im_inset.thumbnail((150,150), Image.ANTIALIAS)
    im_control.thumbnail((150,150), Image.ANTIALIAS)
    
    newWidth = im_scale.size[0] + im_layer.size[0]
    newWidth = im_layer.size[0]
    result = Image.new("RGB", (newWidth,im_scale.size[1]))
    #result.paste(im_scale,(0,0,im_scale.size[0],im_scale.size[1]))
    result.paste(im_layer,(0,0,im_layer.size[0],im_layer.size[1]))
    result.paste(im_inset, (0, 0, im_inset.size[0], im_inset.size[1]))
    result.paste(im_control, (0,
                              result.size[1] - im_control.size[1],
                              im_control.size[0],
                              result.size[1]))
    result.save(basepath + '.png')
    
if __name__=='__main__':
    logging.basicConfig(level=logging.DEBUG)

    # define the layers of interest
    #    [topSurface,bottomSurface,Name,applyMask]
    #    surfaces are numbered from 0 (not 1 as in segmentation software)
    layer_set = [(0,10,'Total',False),
                 (0,1,'RNFL',True),
                 (1,2,'GCL',True),
                 (2,3,'IPL',True),
                 (1,3,'GCL + IPL',True),
                 (1,4,'GCL + IPL + INL',True),
                 (0,3,'RNFL + GCL + IPL',True),
                 (0,4,'RNFL + GCL + IPL + INL',True),
                 (5,8,'Outer',True),
                 (6,8,'PhR',True),
                 (3,4,'INL',True),
                 (4,5,'OPL',True),
                 (5,6,'ONL',True),
                 (1,5,'Inner Complex',True),
                 (5,8,'Outer Complex',True),
                 (4,8,'Outer Complex + OPL',True),
                 (1,10,'Total GCL-Bruchs',True),
                 (5,10,'Outer Retina',True)]
    
    # load the control data
    controls = load_control_data(path = 'Data/Control/',
                                 laterality='OD',
                                 nmax=None)

    data_etdrs = []
    for layer in layer_set:
        tmp_etdrs = controls.getEtdrs(layer[0], layer[1], True)
        tmp_etdrs['layer'] = layer[2]
        tmp_etdrs['group'] = 0  # controls
        tmp_etdrs['id'] = 0
        data_etdrs.append(tmp_etdrs)
    data_etdrs = pd.concat(data_etdrs)
    data_etdrs = pd.melt(data_etdrs, id_vars=['layer','group','id'])
    
    p=re.compile('^P(\d{7}).*(\d{1,2}-\d{1,2}-\d{4}).*(OD|OS)_(sn\d+)')
    p2=re.compile('^(.+)_Surfaces_Iowa.xml')
    
    surface_files = []
    center_files = []    
    # find the files containing the patient data
    for dirname, _, filenames in os.walk('Data/Patient/'):
        for name in fnmatch.filter(filenames,'*Surfaces_Iowa.xml'):
            m = re.match(p,name)
            # convert date format
            dot = datetime.datetime.strptime(m.group(2),'%m-%d-%Y').strftime('%Y-%m-%d')
            surface_files.append({'dirname':dirname,
                                  'filename':name,
                                  'subject':m.group(1),
                                  'dot':dot,
                                  'eye':m.group(3),
                                  'serial':m.group(4)})
        # check if a GridCenter file has been created
        for name in fnmatch.filter(filenames,'*GridCenter_Iowa.xml'):
            center_files.append((dirname,name))     
            
    # sort by subject and test date
    surface_files = sorted(surface_files, key=lambda k: k['subject'] + k['dot'] + k['eye'])
    
    last_subject = None
    # loop through the patient data
    patient_data = []
    for surface_file in surface_files:
        logging.debug('Serial:{}'.format(surface_file['serial']))

        if not last_subject == surface_file['subject']:
            last_eye = None
            last_subject = surface_file['subject']
            # Chapter for each subject
            # tex_file.write('\\chapter{{Subject: {}}}\n'.format(surface_file['subject']))
            
        basename = re.match(p2,surface_file['filename']).group(1)
        try:
            i = [i for i,v in enumerate(center_files) if re.search(basename,v[1])][0]
        except IndexError:
            i = None
        # load the patient data and apply the center file if it exists
        if i is not None:
            patient = Oct.OctLayers(filename = os.path.join(surface_file['dirname'], surface_file['filename']),
                                    center_filename = os.path.join(center_files[i][0],center_files[i][1]))
            patient.centerData()
        else:
            logger.warning('GridCenter file not found for surface file:{}'.format(surface_file['filename']))
            patient = Oct.OctLayers(filename = os.path.join(surface_file['dirname'], surface_file['filename']))
        patient.serialNo = surface_file['serial']
        
        # generate the etdrs data from the patient
        # first we need to set the patient orientation to OD
        orig_orient = surface_file['eye']
        patient.setOrient('OD')
        patient_etdrs = []
        for layer in layer_set:
            tmp_etdrs = patient.getEtdrsThickness(layer[0],layer[1],layer[3])
            tmp_etdrs['layer'] = layer[2]
            tmp_etdrs['group'] = 1 # patients
            tmp_etdrs['id'] = patient.serialNo
            patient_etdrs.append(pd.DataFrame([tmp_etdrs])) #convert from a series to a df with one row
        patient_etdrs = pd.concat(patient_etdrs)
        patient_etdrs = pd.melt(patient_etdrs, id_vars=['layer','group','id'])
        data_etdrs = pd.concat([data_etdrs,patient_etdrs])
        data_etdrs['test']=0
        # reset the orientation
        patient.setOrient(surface_file['eye'])
        # controls.setOrient(surface_file['eye'])
        patient_data.append({'data':patient,
                             'file':surface_file})
        
    # sort the data into test1, test2 pairs
    testPairs = []
    subjects = set([file['subject']+file['eye'] for file in surface_files])
    for subjecteye in subjects:
        tests = [test for test in surface_files if test['subject']+test['eye'] == subjecteye]   
        if len(tests) > 2:
            logger.warning('More than two tests identified for subject:'.format(subjecteye))
        tests = sorted(tests, key=lambda k: k['dot'])
        test1 = tests[0]
        test2 = tests[-1]
        testPairs.append((test1,test2))
        data_etdrs.loc[data_etdrs['id']==test1['serial'],'test']=1
        data_etdrs.loc[data_etdrs['id']==test2['serial'],'test']=2
        data_etdrs.loc[data_etdrs['id']==test1['serial'],'subject']=test1['subject']
        data_etdrs.loc[data_etdrs['id']==test2['serial'],'subject']=test2['subject']
        
    # save the etdrs data for processing in R
    data_etdrs.to_csv('Data/etdrs_data.csv')
    
    # open a tex file for writing
    tex_filename = os.path.join('Output/','analysis' + '.tex')
    tex_file = open(tex_filename, 'w')    
    # write the document header
    tex_file.write('\\documentclass{report}\n')
    tex_file.write('\\usepackage{graphicx}\n')
    tex_file.write('\\usepackage[space]{grffile}\n') # fix for graphics files with spaces
    tex_file.write('\\usepackage{titlesec}\n')
    tex_file.write('\\usepackage{pdflscape}\n')
    #tex_file.write('\\usepackage[landscape]{geometry}\n')
    tex_file.write('\\newcommand{\\sectionbreak}{\\clearpage}\n')
    tex_file.write('\\newcommand{\\subsectionbreak}{\\clearpage}\n')
    tex_file.write('\\title{{{}}}\n'.format('OCT Thickness Analysis'))
    tex_file.write('\\author{Tom Wright}\n')
    tex_file.write('\\date{{{}}}\n'.format(datetime.datetime.now().strftime("%Y-%m-%d")))
    
    tex_file.write('\\begin{document}\n')
    tex_file.write('\\maketitle\n')



    for tests in testPairs:
        data_1 = [patient for patient in patient_data if patient['data'].serialNo == tests[0]['serial']]
        data_2 = [patient for patient in patient_data if patient['data'].serialNo == tests[1]['serial']]
        p=(data_1[0],data_2[0])
        
        controls.setOrient(p[0]['file']['eye'])
        
        tex_basename = '{}-{}'.format(p[0]['file']['subject'],
                                      p[0]['file']['eye'])
        
        #tex_file.write('\t\\section{{Test date:{}}}\n'.format(surface_file['dot']))
        # subsection for each test
        tex_file.write('\t\t\\section{{Patient:{}, Eye:{}}}\n'.format(p[0]['file']['subject'],
                                                                      p[0]['file']['eye']))
        #create a folder to hold the figures
        tex_figure_folder = os.path.join('Output/',tex_basename)
        tex_file.write('\t\t\t\\graphicspath{{{{{}/}}}}\n'.format(tex_basename))
        
        if not os.path.exists(tex_figure_folder ):
            os.mkdir(tex_figure_folder)
            
        # loop through the layers of interest
        for i,layers in enumerate(layer_set):
            for test in [0,1]:
                # generate the figure
                fig_basename = os.path.join(tex_figure_folder,'test' + '_' + str(test) + '_' + layers[2])
                fig_inset_filename = fig_basename + '_inset.png'
                fig_layer_filename = fig_basename + '_layer.png'
                fig_scale_filename = fig_basename + '_scale.png'
                fig_control_filename = fig_basename + '_control.png'
                tex_fig_filename = os.path.join(tex_basename,layers[2]+'.png')
                plotLayerZscores(p[test]['data'], controls, layers[0], layers[1], fig_inset_filename)
                r = plotLayer(p[test]['data'], controls, layers[0], layers[1], fig_layer_filename)
                plotScale(r, fig_scale_filename)
                plotControlLayer(controls, r, layers[0], layers[1], fig_control_filename)
                
                assemble_images(fig_basename)
            
            tex_file.write('\t\t\t\t\\begin{figure}\n')
            tex_file.write('\t\t\t\t\t\\centering\n')
            #if '+' in layers[2]:
                #tex_file.write('\t\t\t\t\t\\includegraphics[width=3.5in]{{"test_0_{}".png}}\n'.format(layers[2]).replace(' ','\space '))
                #tex_file.write('\t\t\t\t\t\\includegraphics[width=3.5in]{{"test_1_{}".png}}\n'.format(layers[2]).replace(' ','\space '))
            #else:
                #tex_file.write('\t\t\t\t\t\\includegraphics[width=3.5in]{{test_0_{}.png}}\n'.format(layers[2]))
                #tex_file.write('\t\t\t\t\t\\includegraphics[width=3.5in]{{test_1_{}.png}}\n'.format(layers[2]))
            tex_file.write('\t\t\t\t\t\\includegraphics[width=3.5in]{{test_0_{}.png}}\n'.format(layers[2]).replace(' ','\space '))
            tex_file.write('\t\t\t\t\t\\includegraphics[width=3.5in]{{test_1_{}.png}}\n'.format(layers[2]).replace(' ','\space '))
            tex_file.write('\t\t\t\t\t\\caption{{layer {} to {}}}\n'.format(layers[0],layers[1]))
            tex_file.write('\t\t\t\t\t\t{}\n'.format(layers[2]).replace('+','$+$'))
            #tex_file.write('\t\t\t\t\t\\hspace{0.5cm}\n')
            tex_file.write('\t\t\t\t\\end{figure}\n')
            tex_file.write('\t\t\t\\clearpage\n')

    # finish and close the tex doc
    tex_file.write('\\end{document}\n')
    tex_file.close()
    