# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 09:38:11 2016

@author: Antoine
"""

import csv
from Tkinter import Tk
import tkFileDialog
import svg_toolbox as SVGT
from geometry1 import IDT2svg
from numpy import deg2rad,array,mean,vstack
#from tkFileDialog import askopenfilename
import scipy.io as sio
import gdsCAD as gds
import shapely.geometry as shapely_geom
import shapely.affinity as shapely_affinity
from shapely.ops import cascaded_union

import matplotlib.pyplot as plt
from descartes import PolygonPatch
from figures import BLUE, RED, BLACK
from shapely.validation import explain_validity
import tkMessageBox
#from tqdm import tqdm


def Create_group():
    Tk().withdraw()
    IDT_group = {'IDT': []}
    IDT_group_dir = tkFileDialog.askdirectory(title = 'Project_filename')
    IDT_group['IDT_group_dir'] = IDT_group_dir
    return IDT_group

def Import_IDT_parameters(IDT_group):
    IDT_group_dir = IDT_group['IDT_group_dir']
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    csv_filename = tkFileDialog.askopenfilename(title = 'IDT design file ?', defaultextension = 'csv',initialdir = IDT_group_dir)
    slowness_database = {}   
    with open(csv_filename) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            parameters = Read_IDT_param_from_row(row)
            IDT_data = {'parameters' : parameters}
            IDT_data['parameters']['reticule_filename']=IDT_group_dir + '/' +parameters['reticule_filename'] + '.svg'     
            key =parameters['slowness_substrate']            
            if key not in slowness_database.keys():
                try:
                    slowness_filename = IDT_group_dir + '/' + key + '.mat'
                    slowness_database[key] = Import_slowness_from_matlab_NoGui(slowness_filename)
                except IOError:
                    slowness_database[key] = Import_slowness_from_matlab(IDT_group,query = "File %s not found" %key)
            IDT_data['parameters']['slowness_substrate'] = slowness_database[key]
                    #CONTINUE FROM HERE (exception handling file not found)
            IDT_group['IDT'].append(IDT_data)
    
def Import_slowness_from_matlab(IDT_group, query = 'Matlab substrate properties'):
    IDT_group_dir = IDT_group['IDT_group_dir']
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    matlab_filename = tkFileDialog.askopenfilename(title=query, defaultextension = 'mat',initialdir = IDT_group_dir);
    matlab_results = sio.loadmat(matlab_filename)
    matlab_results['filename'] = matlab_filename
    matlab_results['psi']=matlab_results['psi'].squeeze()
    matlab_results['a']=matlab_results['a'].squeeze()
    if mean(matlab_results['slowness'])>0.01:
        matlab_results['s_Ray']=1e-3*matlab_results['slowness'].squeeze()#mm/us
    else:
        matlab_results['s_Ray']=matlab_results['slowness'].squeeze()
    #IDT_group['substrate_properties'] = matlab_results  
    return matlab_results

def Import_slowness_from_matlab_NoGui(matlab_filename):
    matlab_results = sio.loadmat(matlab_filename)
    matlab_results['filename'] = matlab_filename
    matlab_results['psi']=matlab_results['psi'].squeeze()
    matlab_results['a']=matlab_results['a'].squeeze()
    if mean(matlab_results['slowness'])>0.01:
        matlab_results['s_Ray']=1e-3*matlab_results['slowness'].squeeze()#mm/us
    else:
        matlab_results['s_Ray']=matlab_results['slowness'].squeeze()
    return matlab_results        
            
def Read_IDT_param_from_row(row):
    parameters = {'slowness':[],'z':[],'layers':{},'electrodes_angle':[]}
    for key in row.keys():
        if 'zlayer' in key:
            ind = key.find('_')
            layer_name = key[(ind+1):]
            if layer_name not in parameters['layers']:
                parameters['layers'][layer_name] = {'z':-1.0,'slowness':-1.0}
            parameters['layers'][layer_name]['z'] = float(row[key])
        elif 'slowness' in key:
            ind = key.find('_')
            layer_name = key[(ind+1):]
            if layer_name not in parameters['layers']:
                parameters['layers'][layer_name] = {'z':-1.0,'slowness':-1.0}
            try:
                parameters['layers'][layer_name]['slowness'] = float(row[key])
            except ValueError:
                parameters['layers'][layer_name]['slowness'] = row[key] #the filename of the mat-file of the slowness
        elif 'electrodes_angle' in key:
            parameters['electrodes_angle'].append(deg2rad(float(row[key])))
        else:
            try:
                parameters[key] = float(row[key]) 
            except ValueError:
                parameters[key] = row[key]
    #now we need to sort the layers based on their z
    layer_seq = []
    for layer_name in parameters['layers']:
        layer_seq.append((parameters['layers'][layer_name]['z'],parameters['layers'][layer_name]['slowness']))
    layer_seq.sort(key=lambda foo: foo[0])
    for layer in layer_seq:   
        parameters['z'].append(layer[0])
        parameters['slowness'].append(layer[1])
    parameters['z'] = array( parameters['z'])
    parameters['slowness_substrate'] = str(parameters['slowness'][0])
    parameters['slowness'] = array(parameters['slowness'][1:])
    return parameters
                 
def IDTgroup2svg(IDT_group):
    IDT_group_dir = IDT_group['IDT_group_dir']
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    svg_filename = tkFileDialog.asksaveasfilename(title='Wafer routing filename',defaultextension = 'svg',initialdir = IDT_group_dir);
    svg_dwg = SVGT.Init_svg(svg_filename)
    for IDT_data in IDT_group['IDT']:    
        IDT2svg(IDT_data,svg_dwg)
    svg_dwg.save();    
    IDT_group['svg_route']=svg_filename     

def importSVGroute(IDT_group):
    IDT_group_dir = IDT_group['IDT_group_dir']
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    svg_filename = tkFileDialog.askopenfilename(title='Wafer routing filename',defaultextension = 'svg',initialdir = IDT_group_dir);
    all_points = SVGT.read_colored_path_from_svg(svg_filename)
    route =[[],[]]
    for points in all_points['red']:
        polygon = shapely_geom.Polygon(points.T)
        polygon_validity = explain_validity(polygon)
        if polygon_validity=='Valid Geometry':
            route[0].append(polygon)
        else:
            tkMessageBox.showwarning('Error in svg import', polygon_validity)
    for points in all_points['blue']:
        polygon = shapely_geom.Polygon(points.T)
        polygon_validity = explain_validity(polygon)
        if polygon_validity=='Valid Geometry':
            route[1].append(polygon)
        else:
            tkMessageBox.showwarning('Error in svg import', polygon_validity)
    route[0] = shapely_geom.MultiPolygon(route[0])
    route[1] = shapely_geom.MultiPolygon(route[1])
    #outbox = route[0].bounds
    #dx = outbox[2]-outbox[0]
    #dy = outbox[3]-outbox[1]
    #x0 = outbox[0]+dx/2
    #y0 = outbox[1]+dy/2
    x0 = 4000
    y0 = 4000
    factor = 1e-5;
    route[0] = shapely_affinity.translate(route[0], xoff=-x0, yoff=-y0) 
    route[0] = shapely_affinity.scale(route[0], xfact = factor, yfact= factor, origin=(0,0,0))
    route[1] = shapely_affinity.translate(route[1], xoff=-x0, yoff=-y0) 
    route[1] = shapely_affinity.scale(route[1], xfact = factor, yfact= factor, origin=(0,0,0))
    IDT_group['route'] = route 

def combineIDT_group(IDT_group):
    print('combining IDTs, this may take several minutes')
    plot = False
    route0 = []
    route1 = []
    for polygon in IDT_group['route'][0]:
        route0.append(polygon)
    for polygon in IDT_group['route'][1]:   
        route1.append(polygon)
    for IDT_data in IDT_group['IDT']:
        removal0 = IDT_data['removal']['removal0']
        removal1 = IDT_data['removal']['removal1']
        removal_center = IDT_data['removal']['removal_center']
        IDT = IDT_data['IDT']
        for ind_route0 in range(len(route0)):
            if IDT_data['convex_hull'].intersects(IDT_group['route'][0][ind_route0]):
                branch0 = route0[ind_route0]
                
                if plot:
                    ax = plt.subplot(111)
                    for polygon in [branch0]+IDT[0]:
                         patch = PolygonPatch(polygon, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2)
                         ax.add_patch(patch)
                    for polygon in route1:
                         patch = PolygonPatch(polygon, fc=RED, ec=RED, alpha=0.5, zorder=2)
                         ax.add_patch(patch)
                    plt_bnds = branch0.bounds
                    ax.axis([plt_bnds[0],plt_bnds[2], plt_bnds[1], plt_bnds[3]])
                    ax.set_aspect(1)
                    plt.show()  
                
                #for polygon in IDT[0]:
                #    branch0 = branch0.union(polygon)
                branch0 = cascaded_union([branch0]+IDT[0])
                branch0 = branch0.difference(removal1)
                if 'Hole_elect0' in IDT_data['removal'].keys():
                    branch0 = branch0.difference(IDT_data['removal']['Hole_elect0'])
                if type(branch0) is shapely_geom.multipolygon.MultiPolygon:
                    #the difference is made on the wrong side of the IDT. find a way to switch the polarities on time
                    IDT = IDT[::-1] #switch polarities
                    removal0 = IDT_data['removal']['removal1']
                    removal1 = IDT_data['removal']['removal0']
                    branch0 = route0[ind_route0]
                    branch0 = cascaded_union([branch0]+IDT[0])
                    branch0 = branch0.difference(removal1)
                    if 'Hole_elect1' in IDT_data['removal'].keys():
                        branch0 = branch0.difference(IDT_data['removal']['Hole_elect1'])
                branch0 = branch0.difference(removal_center)
                #print "fusion finished, the new polygon has %d vertices" %len(branch0.exterior)
                route0[ind_route0]=branch0
        for ind_route1 in range(len(route1)):
            if IDT_data['convex_hull'].intersects(IDT_group['route'][1][ind_route1]):
                branch1 = route1[ind_route1] 
                #for polygon in IDT[1]:
                #    branch1 = branch1.union(polygon)
                
                if False:
                    print(explain_validity(branch1))
                    #plt_bnds = IDT_data['convex_hull'].bounds
                    plt_bnds = branch1.bounds
                    ax = plt.subplot(111)
                    for polygon in [branch1]+IDT[1]:
                         patch = PolygonPatch(polygon, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2)
                         ax.add_patch(patch)
                    ax.axis([plt_bnds[0],plt_bnds[2], plt_bnds[1], plt_bnds[3]])
                    ax.set_aspect(1)
                    plt.show() 

                branch1 = cascaded_union([branch1]+IDT[1])
                branch1 = branch1.difference(removal0)
                if 'Hole_elect1' in IDT_data['removal'].keys():
                    branch1 = branch1.difference(IDT_data['removal']['Hole_elect1'])
                if type(branch1) is shapely_geom.multipolygon.MultiPolygon:
                    #the difference is made on the wrong side of the IDT. find a way to switch the polarities on time
                    IDT = IDT[::-1] #switch polarities
                    removal0 = IDT_data['removal']['removal1']
                    removal1 = IDT_data['removal']['removal0']
                    branch1 = route1[ind_route1]
                    branch1 = cascaded_union([branch1]+IDT[1])
                    branch1 = branch1.difference(removal0)
                    if 'Hole_elect0' in IDT_data['removal'].keys():
                        branch1 = branch1.difference(IDT_data['removal']['Hole_elect0'])
                branch1 = branch1.difference(removal_center)
                route1[ind_route1]=branch1
    IDT_group['final_IDT'] = route0+route1        
    #for polygon in route:        
    #    IDT_group['final_IDT'].append(polygon)
        
def IDT_group2gds(IDT_group):
    IDT_group_dir = IDT_group['IDT_group_dir']
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    gds_filename = tkFileDialog.asksaveasfilename(title='Wafer gds filename',defaultextension = 'gds',initialdir = IDT_group_dir);
    layout = gds.core.Layout('LIBRARY')
    cell = gds.core.Cell('Main')
    for polygon in IDT_group['final_IDT']:
        points_exterior =  1e6*array(polygon.exterior.coords) #best practice: dimensions in um        
        points_interior_has_to_be_created  = True       
        for interior_ring in polygon.interiors:
            if points_interior_has_to_be_created:
                points_interior = 1e6*array(interior_ring.coords)
                points_interior_has_to_be_created = False
            else:
                points_interior = vstack((points_interior,1e6*array(interior_ring.coords)))       
        if not polygon.interiors:
            boundaries = gds.core.Boundary(points_exterior) 
        else:
            boundaries = gds.core.Boundary(vstack((points_interior,points_exterior))) 
        cell.add(boundaries)
    for IDT_data in IDT_group['IDT']:
        reticule = IDT_data['reticule']
        for polygon in reticule:
            points = 1e6*array(polygon.exterior.coords)
            boundaries = gds.core.Boundary(points)
            cell.add(boundaries)
    layout.add(cell)
    layout.save(gds_filename)        
    
def IDT_group2svg(IDT_group):
    IDT_group_dir = IDT_group['IDT_group_dir']
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    svg_filename = tkFileDialog.asksaveasfilename(title='SVG export filename',defaultextension = 'svg',initialdir = IDT_group_dir);
    plt_bnds = [0,8000,0,8000]
    x0 = 40e-3
    y0 = 40e-3
    factor = 1e5
    fig, ax = plt.subplots()    
    for polygon in IDT_group['final_IDT']:
        polygon = shapely_affinity.translate(polygon, xoff=x0, yoff=y0) 
        polygon = shapely_affinity.scale(polygon, xfact = factor, yfact= factor, origin=(0,0,0))
        patch = PolygonPatch(polygon, fc=BLACK, ec=None, alpha=1.0, zorder=2)
        ax.add_patch(patch)
    for IDT_data in IDT_group['IDT']:
        reticule = IDT_data['reticule']
        for polygon in reticule:
            polygon = shapely_affinity.translate(polygon, xoff=x0, yoff=y0) 
            polygon = shapely_affinity.scale(polygon, xfact = factor, yfact= factor, origin=(0,0,0))
            patch = PolygonPatch(polygon, fc=BLACK, ec='none', alpha=1.0, zorder=2)
            ax.add_patch(patch)
    ax.axis(plt_bnds)
    ax.set_aspect(1)
    #plt.show() 
    fig.savefig(svg_filename, format='svg', dpi=1200)