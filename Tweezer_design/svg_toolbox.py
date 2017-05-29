# -*- coding: utf-8 -*-
"""
Created on Sat Apr 02 18:49:30 2016
version 2.0
@author: Antoine
 """


import numpy as np
from lxml import etree
import svgwrite 
#import Tkinter
import tkMessageBox

def string2points(string2convert):
    #str.find() finds the first occurence of a str in another one, and returns -1 if no such string is found
    Points = [];    
    for i in range(len(string2convert)): 
        string2convert[i] = string2convert[i].replace(',',' ');
        string2convert[i] = string2convert[i].split(' ');#results in a list of list of strings
        x = [];
        y = [];
        current_iter = 0;  
        Continue = False;
        for j in range(0,len(string2convert[i])):
            if string2convert[i][j]=='M':
                if Continue:
                    Points.append(np.array([x,y]))
                    x = [];
                    y = [];
                    current_symbol = 'M';
                else:
                    current_symbol = 'L';
            elif string2convert[i][j]=='m':
                if Continue:
                    Points.append(np.array([x,y]))
                    x = [];
                    y = [];
                    current_symbol = 'm';
                else:
                    current_symbol = 'l';
            elif string2convert[i][j]=='l':
                current_symbol = 'l';
            elif string2convert[i][j]=='L':
                current_symbol = 'L';
            elif string2convert[i][j]=='C':
                current_symbol = 'C';
            elif string2convert[i][j]=='c':
                current_symbol = 'c';
            elif string2convert[i][j]=='z':
                Points.append(np.array([x,y]));
                if j==(len(string2convert[i])-1):
                    break
                else:
                    if string2convert[i][j+1] in ('>','/','"',''):
                        break
                    else:
                        x = [];
                        y = [];
                        current_iter = 0;
                        Continue = False;
            else:
                if current_iter ==0:
                    current_iter = 1;
                    if (current_symbol == 'M') | (current_symbol == 'L'):
                        x.append(float(string2convert[i][j]))
                    elif current_symbol == 'l':
                        if Continue:
                            x.append(float(string2convert[i][j])+x[-1]);
                        else:
                            x.append(float(string2convert[i][j]));
                    elif (current_symbol == 'C')|(current_symbol == 'c'):
                        current_iter = 0;
                    else:
                        if Continue:
                            x.append(float(string2convert[i][j])+Points[-1][0,-1]);
                        else:
                            x.append(float(string2convert[i][j]));
                elif current_iter==1:
                    current_iter = 0;
                    if (current_symbol == 'M') | (current_symbol == 'L'):
                        y.append(float(string2convert[i][j]))
                    elif current_symbol == 'l':
                        if Continue:
                            y.append(float(string2convert[i][j])+y[-1]); 
                        else:
                            y.append(float(string2convert[i][j]));
                    elif (current_symbol == 'C')|(current_symbol == 'c'):
                        current_iter = 1;
                        Continue = False
                    else:
                        if Continue:
                            y.append(float(string2convert[i][j])+Points[-1][1,-1]);
                        else:
                            y.append(float(string2convert[i][j]));   
                    Continue = True;
                elif current_iter == 5:
                     current_iter = 6;
                     if (current_symbol == 'C'):
                        x.append(float(string2convert[i][j]))
                     else:
                        x.append(float(string2convert[i][j])+x[-1])
                elif current_iter == 6:
                    current_iter = 0;
                    if (current_symbol == 'C'):
                       y.append(float(string2convert[i][j]))
                    else:
                        y.append(float(string2convert[i][j])+y[-1])
                    Continue = True
                else:
                    current_iter = current_iter +1;
    return Points

def read_path_from_svg(svg_filename):
   # re_split = re.compile('\s+|,')
    try:
        tree = etree.parse(open(svg_filename))
    except (OSError, IOError):
        tkMessageBox.showwarning(title = 'Error during reticule import',message = "the file %s was not found" %svg_filename)   
    path_list = [];
    for element in tree.iter():
        if element.tag.split("}")[1] == "path":
            #re_split.split(element.get("d"))
            path_list.append(element.get('d'))
            #path_list.append(re_split)  
    Points = string2points(path_list)
    return Points

def read_colored_path_from_svg(svg_filename):
   # re_split = re.compile('\s+|,')
    tree = etree.parse(open(svg_filename))
    path_list = {}
    for element in tree.iter():
        if element.tag.split("}")[1] == "path":
            #re_split.split(element.get("d"))
            path_color = element.get('style')
            if type(path_color) == str:
                ind = path_color.find('fill:')
                path_color = path_color[(ind+5):(ind+12)]
                if path_color == '#ff0000':
                    path_color = 'red'
                elif path_color == '#00ff00':
                    path_color = 'green'
                elif path_color == '#0000ff':
                    path_color = 'blue'
                else:
                    tkMessageBox.showwarning('Some colors were not rekognized', 'Only blue, red and green traces are imported')
            else:
                path_color = element.get('fill')
            if path_color not in path_list.keys():
                path_list[path_color]=[]
            path_list[path_color].append(element.get('d'))
            #path_list.append(re_split)
    Points = {}
    for path_color in path_list:        
        Points[path_color] = string2points(path_list[path_color])
    return Points    

def Init_svg(output_name):
    svg_Lx = 8000;
    svg_Ly = 8000;
    dwg = svgwrite.Drawing(output_name,(svg_Lx,svg_Ly),debug=True)
    #dwg.add(dwg.rect(insert=(0,0),size = ('100%','100%'),rx = None,ry = None,fill ='white'))
    return dwg
    
def Array2svg(svg_dwg,Array,simplify_optn,stroke2='red',fill2='red'):  
    offx = 4000
    offy = 4000
    diameter = 220 #*10um
    step_min = diameter*np.pi/20
    start = [Array[0,0]+offx,Array[0,1]+ offy] ;
    s1 = 'M {0[0]} {0[1]}'.format(start)
    p = svg_dwg.path(d=s1,stroke_width=1,stroke = stroke2,fill = fill2);
    points_seq = [];
    assert Array.shape[1] is 2, "Input Array %r should be an Nx2 array" %Array 
    assert np.amax(np.sqrt(Array[:,0]**2+Array[:,1]**2))< 5000, "Input Array %r elements are expected to be in 10*um" %Array
    for i in range(1,Array.shape[0]):
        if simplify_optn is 'simplify':
            if len(points_seq)<2:
                current_step = np.sqrt((Array[0,0]-Array[i,0])**2 + (Array[0,1]-Array[i,1])**2)
            else:
                current_step = np.sqrt((points_seq[-2] - (Array[i,0]+offx))**2 + (points_seq[-1] - (Array[i,1]+offx))**2)           
            #print(current_step)            
            if current_step > step_min:            
                points_seq.append(Array[i,0]+offx)
                points_seq.append(Array[i,1]+offy)
        else:
            points_seq.append(Array[i,0]+offx)
            points_seq.append(Array[i,1]+offy)
    points_seq.append('z') 
    #print(points_seq)
    p.push('L',points_seq)
    svg_dwg.add(p)
    return svg_dwg 