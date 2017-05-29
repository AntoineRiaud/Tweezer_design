# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 12:25:02 2016

@author[code]: Antoine Riaud
@coauthors: Michael Baudoin, Jean-Louis Thomas, Olivier Bou Matar
IEMN - Institut d'Electronique, Microelectronique et Nanotechnologies
INSP - Institut des Nanosciences de Paris 
"""

#function summary
# Default_properties[Debug] emulates the properties of a piezo material

#import goes here
import math
import numpy
import shapely.geometry as shapely_geom
import shapely.affinity as shapely_affinity
import shapely.ops as shapely_ops
import matplotlib.pyplot as plt
from filtered_der import periodic_derivative
from descartes import PolygonPatch
from figures import BLUE, RED, GREY
from copy import deepcopy
import scipy.io as sio
import gdsCAD as gds
import svg_toolbox as SVGT
import tkMessageBox



#material properties import
def Default_properties(anisotropy_factor):
    plot=False
    n_psi = 100
    i = range(n_psi)
    psi = numpy.array(i)*math.pi/n_psi
    s_Ray = 1/(3500*(1+anisotropy_factor*numpy.sin(2*psi)))
    a = numpy.ones(n_psi)
    a[25:]=-1
    if plot:
        ax = plt.subplot(111, projection='polar')
        ax.plot(psi.tolist(), s_Ray.tolist(), color='r', linewidth=3)
        ax.set_rmax(max(s_Ray))
        ax.grid(True)
        
        ax.set_title("Slowness curve", va='bottom')
        plt.show()
    return({'psi': psi, 's_Ray': s_Ray, 'a': a})

def cart2pol(xy):
    x = xy[:,0]
    y = xy[:,1]
    rho = numpy.sqrt(x**2 + y**2)
    phi = numpy.arctan2(y, x)
    return(rho, phi)

def pol2cart(phi,rho):
    x = rho * numpy.cos(phi)
    y = rho * numpy.sin(phi)
    xy = numpy.vstack((x,y)).T
    return(xy)

def Degeneration_correction(z,slowness,s_Ray,psi):
    # z and slowness are numpy arrays
    # the 0th dimension is for the propagation angle psi
    # the 1st dimension is for the material layer
    plot=False
    slowness = numpy.atleast_2d(slowness)
    z = numpy.atleast_2d(z)
    #s_Ray2 = numpy.atleast_2d(s_Ray)
    if z.size <=1:
        mu_0 = 0*psi
    else:
        if slowness.shape[0]==len(psi):      
            sz = numpy.sqrt(slowness**2 - numpy.tile(s_Ray**2,(slowness.shape[1],1)).T)
            dz = numpy.diff(z)
            dz = numpy.tile(dz,(psi.size,1))
            mu_0 = numpy.sum(dz*sz,1)
        else:
            sz = numpy.sqrt(numpy.tile(slowness**2,(s_Ray.size,1)) - numpy.tile(s_Ray**2,(slowness.size,1)).T)
            mu_0 = numpy.sum(numpy.diff(z,1,1)*sz,1)
    if plot:
        ax = plt.subplot(111, projection='polar')
        ax.plot(psi.tolist(), mu_0.tolist(), color='r', linewidth=3)
        ax.set_rmax(max(mu_0))
        ax.grid(True)
        ax.set_title("Correction", va='bottom')
        plt.show()
    return(mu_0)

def Mesh_Theta(target_size,omega,mu_0,s_Ray,N_turns,l,electrodes_angle):
    freq = omega/(2*numpy.pi)
    lambda_approx = 1/(freq*numpy.mean(s_Ray))
    R0 = target_size + omega*numpy.mean(mu_0)*lambda_approx/(2*numpy.pi)
    if l==0:
        Theta = []
        phi_0 = 2*numpy.pi*R0/lambda_approx
        for n in range(N_turns):
            Theta_max = 2*numpy.pi
            R_approx = R0+n*lambda_approx
            dTheta = lambda_approx/(6*R_approx)
            Theta.append(numpy.arange(electrodes_angle[0],Theta_max+electrodes_angle[0],dTheta))
    else:
        Theta_max = 2*numpy.pi*numpy.ceil(N_turns/numpy.absolute(l))
        Theta = [electrodes_angle[0]]
        while Theta[-1]<(electrodes_angle[0]+Theta_max):
            R_approx = R0+numpy.abs(l)*Theta[-1]*lambda_approx/(2*numpy.pi)
            dTheta = lambda_approx/(6*R_approx)
            Theta.append(Theta[-1]+dTheta)
        if l<0:
            Theta = Theta[::-1]
            phi_0 = numpy.abs(l)*(Theta[0]+2*numpy.pi*R0/lambda_approx)
        else:
            phi_0 = (2*numpy.pi*R0/lambda_approx)
    return {'Theta':Theta,'phi_0':phi_0}

def thickness_factor(electrode_type):
    return {
        'IDT':1.0,
        'splitIDT':0.5
    }[electrode_type]

def IDT_Master_Curve(psi,Theta,mu_0,s_Ray,l,omega,phi_0,a):
    plot = False
    save = False  
    
    sprime_Ray = periodic_derivative(psi,s_Ray,16)
    cos_beta = s_Ray/numpy.sqrt(s_Ray**2+sprime_Ray**2)
    sin_beta = sprime_Ray/numpy.sqrt(s_Ray**2+sprime_Ray**2)
    beta = numpy.arctan2(sin_beta,cos_beta)
    psibar_Theta = Theta+numpy.interp(Theta,psi,beta,period = 2*numpy.pi)
    Num1 = numpy.interp(Theta,psi,-mu_0*omega,period=2*numpy.pi)
    Num2 = l*numpy.array(Theta)
    Num3r = numpy.interp(psibar_Theta,psi,+1.0*numpy.real(a),period=2*numpy.pi)
    Num3i = numpy.interp(psibar_Theta,psi,+1.0*numpy.imag(a),period=2*numpy.pi)
    if numpy.sum(Num3r**2)>10*numpy.sum(Num3i**2):
        Num3 = (1-numpy.sign(Num3r))*numpy.pi/2
    else:       
        Num3 = numpy.angle(Num3r+1j*Num3i)
    #Num3 = Num3 - 2*numpy.pi*round(numpy.mean(Num3)/(2*numpy.pi))
    
    Den1 = numpy.interp(psibar_Theta,psi,omega*s_Ray,period=2*numpy.pi)
    Den2 = numpy.interp(Theta,psi,cos_beta,period = 2*numpy.pi)
    #phi_0 = 120.91
    R = (phi_0 +Num1 + Num2 +Num3)/(Den1 * Den2)
    if plot:
        Thetaplot = numpy.linspace(Theta[0],Theta[-1],num=5000)
        if l==0:
            Rplot =  numpy.interp(Thetaplot,Theta,R)
        else:
            Rplot =  numpy.interp(Thetaplot,Theta[::-1],R[::-1])
            Den1Den2plot =  numpy.interp(Thetaplot,Theta[::-1],(Den1*Den2)[::-1])
            phi_0Num1Num2plot =  numpy.interp(Thetaplot,Theta[::-1],(phi_0+Num1+Num2+Theta)[::-1])
            Num3plot =  numpy.interp(Thetaplot,Theta[::-1],(Num3)[::-1])
            betaplot = numpy.interp(Thetaplot,Theta[::-1],(psibar_Theta-Theta)[::-1])
        #ax = plt.subplot(111, projection='polar')
        ax = plt.subplot(221)
        ax.plot(Thetaplot.tolist(), (Den1Den2plot).tolist(), color='r', linewidth=3)
        ax = plt.subplot(222)
        ax.plot(Thetaplot.tolist(), (betaplot).tolist(), color='r', linewidth=3)
        ax = plt.subplot(223)
        ax.plot(Thetaplot.tolist(), (Num3plot).tolist(), color='r', linewidth=3)
        ax = plt.subplot(224)
        ax.plot(Thetaplot.tolist(), (phi_0Num1Num2plot).tolist(), color='r', linewidth=3)
        #ax.set_rmax(max(Rplot))
        ax.grid(True)
        ax.set_title("R", va='bottom')
        plt.show()
        
    if save:
        sio.savemat('OutputR',{'R':R,'Theta':Theta,'phi_0':phi_0,'Num1':Num1,'Num2':Num2,'Num3':Num3,'Den1':Den1,'Den2':Den2})
        
    return R

def Electrode_pattern(electrode_type,phi_0):
    return{
    'IDT': 
        [phi_0+2.*numpy.pi*numpy.array([0]),phi_0+2.*numpy.pi*numpy.array([0.5])],
    'splitIDT': 
        [phi_0+2.*numpy.pi*numpy.array([0,0.25]),phi_0+2.*numpy.pi*numpy.array([0.5,0.75])],
    }[electrode_type]

def MasterCurve2Electrodes(electrode_type,phi_0,psi,Theta,mu_0,s_Ray,l,omega,a,lambda_approx):
    #plot = False    
    Phi_0 = Electrode_pattern(electrode_type,phi_0)
    #N = Phi_0[0].size+Phi_0[1].size #number of electrodes over lambda
    #IDT = [[],[]]
    Rlist = [[],[]]
    Thetalist = [[],[]]
    for polarity in range(2):
        for i in range(len(Phi_0[polarity])):
            if l==0:
                n_turns = len(Theta)
                for n in range(n_turns):
                    R = IDT_Master_Curve(psi,Theta[n],mu_0,s_Ray,0,omega,Phi_0[polarity][i]+n*2*numpy.pi,a)                     
                    Rlist[polarity].append(R) 
                    Thetalist[polarity].append(Theta[n]) 
                    #xy = pol2cart(Theta[n],R)
                    #Rline = shapely_geom.asLineString(xy)
                    #IDT[polarity].append(Rline.buffer(lambda_approx/(4*N))) #for a coverage of 50%
            else:
                if l<0:
                    for n in range(abs(l)):
                        R = IDT_Master_Curve(psi,Theta,mu_0,s_Ray,l,omega,Phi_0[polarity][i]-n*2*numpy.pi,a) 
                        Rlist[polarity].append(R) 
                        Thetalist[polarity].append(Theta) 
                        #xy = pol2cart(Theta,R)
                        #Rline = shapely_geom.asLineString(xy)
                        #IDT[polarity].append(Rline.buffer(lambda_approx/(4*N))) #for a coverage of 50%
                else:
                    for n in range(l):
                        R = IDT_Master_Curve(psi,Theta,mu_0,s_Ray,l,omega,Phi_0[polarity][i]+n*2*numpy.pi,a) 
                        Rlist[polarity].append(R) 
                        Thetalist[polarity].append(Theta) 
                        #xy = pol2cart(Theta,R)
                        #Rline = shapely_geom.asLineString(xy)
                        #IDT[polarity].append(Rline.buffer(lambda_approx/(4*N))) #for a coverage of 50%
                
    
    IDT_data = {'R':Rlist,'Theta':Thetalist} #'IDT': IDT, 
    return IDT_data
        
def Place_Electrodes(l,phi_0,lambda_approx,electrodes_angle,electrodes_thickness,removal_thickness,IDT_data,electrode_type,electrode_rotate_angle=0):
    plot = False
    Rlist = IDT_data['R']
    Thetalist = IDT_data['Theta']
    Rmin = []
    Rmax = []
    Rline = [[],[]]
    for polarity in range(len(Rlist)):
        for i in range(len(Rlist[polarity])):
            R = Rlist[polarity][i]
            Theta = Thetalist[polarity][i]
            Rmin = min([R.min(),Rmin])
            Rmax = max([R.max(),Rmax])
            xy = pol2cart(Theta,R)
            if l==0:
                #Rline[polarity].append(shapely_geom.LinearRing(shapely_geom.asLineString(xy)))
                Rline[polarity].append(shapely_geom.asLinearRing(xy))
            else:
                Rline[polarity].append(shapely_geom.asLineString(xy))
    Rmax = 1.5*Rmax
    xyelect0 = pol2cart([electrodes_angle[0],electrodes_angle[0]],[Rmin,Rmax])
    xyelect1 = pol2cart([electrodes_angle[1],electrodes_angle[1]],[Rmin,Rmax])
    elect0 = shapely_geom.LineString([(xyelect0[0,0],xyelect0[0,1]),(xyelect0[1,0],xyelect0[1,1])])
    elect1 = shapely_geom.LineString([(xyelect1[0,0],xyelect1[0,1]),(xyelect1[1,0],xyelect1[1,1])])
    elect0 = shapely_affinity.rotate(elect0, electrode_rotate_angle, origin=(xyelect0[0,0],xyelect0[0,1]), use_radians=False) 
    elect1 = shapely_affinity.rotate(elect1, electrode_rotate_angle, origin=(xyelect1[0,0],xyelect1[0,1]), use_radians=False)
    removal0 = deepcopy(elect0)
    removal1 = deepcopy(elect1)
    
    removal_center = shapely_geom.Point(0.,0.).buffer(Rmin-electrodes_thickness*thickness_factor(electrode_type)/2)
    
    Elect0_polyg = elect0.buffer(electrodes_thickness/2) 
    Elect1_polyg = elect1.buffer(electrodes_thickness/2)
    
    IDT_data['removal'] = {'removal0':removal0.buffer(removal_thickness/2),
                            'removal1':removal1.buffer(removal_thickness/2),
                            'removal_center':removal_center}  
                            
    if False:#electrode_rotate_angle != 0. : #draw hollow electrode to avoid forcing opposite voltage to the wave
        Hole_elect0 = elect0.buffer(electrodes_thickness/4) #up to 50% of the power branch can be hollow
        Hole_elect1 = elect1.buffer(electrodes_thickness/4)
        Phi_0 = Electrode_pattern(electrode_type,phi_0)
        N = Phi_0[0].size+Phi_0[1].size #number of electrodes over lambda
        Hole_elect0 = Hole_elect0.intersection(shapely_geom.MultiLineString(Rline[1]).buffer(lambda_approx/(4*N)))
        Hole_elect1 = Hole_elect1.intersection(shapely_geom.MultiLineString(Rline[0]).buffer(lambda_approx/(4*N)))
        IDT_data['removal']['Hole_elect0'] = Hole_elect0
        IDT_data['removal']['Hole_elect1'] = Hole_elect1
        #Elect0_polyg  = Elect0_polyg.difference(Hole_elect0)  
        #Elect1_polyg  = Elect1_polyg.difference(Hole_elect1)
  
         
    Rline = clean_edges(Rline,elect0,elect1)     
    
    IDT_data['electrodes'] = {'elect0':Elect0_polyg,
                            'elect1':Elect1_polyg}                        
    IDT_data['Rline'] = Rline
    
    if plot:
        ax = plt.subplot(111)
        patch1 = PolygonPatch(elect0.buffer(electrodes_thickness/2), fc=RED, ec=RED, alpha=0.5, zorder=2)
        ax.add_patch(patch1)
        patch2 = PolygonPatch(elect1.buffer(electrodes_thickness/2), fc=BLUE, ec=BLUE, alpha=0.5, zorder=2)
        ax.add_patch(patch2)    
        ax.set_title('a) dilation, cap_style=3')
        xrange = [-0.005, +0.005]
        yrange = [-0.005, +0.005]
        ax.axis([xrange[0],xrange[1], yrange[0], yrange[1]])
        ax.set_aspect(1)
        plt.show()
    #return IDT_data



def clean_edges(Rlinelist,elect0,elect1):
    plot = False
    Electrodes = shapely_geom.MultiLineString([elect0,elect1])
    Rline_cleaned = [[],[]]
    for polarity in range(len(Rlinelist)):
        #Rline = shapely_geom.collection.GeometryCollection(Rlinelist[polarity])   
        for n in range(len(Rlinelist[polarity])):
            Rline =  Rlinelist[polarity][n]
            Rsegments = Rline.difference(Electrodes)
            check_linear_ring = Rsegments.geoms[0].coords[0]==Rsegments.geoms[-1].coords[-1]
            check_no_intersection = not shapely_geom.Point(Rsegments.geoms[0].coords[0]).buffer(1e-9).intersects(Electrodes)
            if check_linear_ring and check_no_intersection: #there seems to be a bug in the linear ring intersection of shapely
                Rsegments_geoms_0 = shapely_ops.linemerge([Rsegments.geoms[0],Rsegments.geoms[-1]])
                Rsegments_temp = [Rsegments_geoms_0]                
                for i in range(len(Rsegments.geoms)-2):
                    Rsegments_temp.append(Rsegments[i+1])
                Rsegments = shapely_geom.MultiLineString(Rsegments_temp)
            for segment in Rsegments.geoms:
                
                if plot:
                    ax = plt.subplot(111)
                    segmentarray = numpy.array(segment.coords)
                    elect0array = numpy.array(elect0.coords)
                    elect1array = numpy.array(elect1.coords)
                    ax.plot(segmentarray[:,0].tolist(), segmentarray[:,1].tolist(), color='k', linewidth=3)
                    ax.plot(elect0array[:,0].tolist(), elect0array[:,1].tolist(), color='b', linewidth=3)
                    ax.plot(elect1array[:,0].tolist(), elect1array[:,1].tolist(), color='r', linewidth=3)
                    ax.axis([min(segmentarray[:,0]),max(segmentarray[:,0]),min(segmentarray[:,1]), max(segmentarray[:,1])])   
                    ax.set_aspect(1)
                    ax.grid(True)
                    plt.show()
                
                if segment.buffer(1e-9).intersects(elect0) and segment.buffer(1e-9).intersects(elect1): #there seems to be a numerical issue with small segments
                    Rline_cleaned[polarity].append(segment)

    return Rline_cleaned


def Enlarge_IDT(IDT_data,phi_0,electrode_type,lambda_approx): 
    plot = False    
    #Rlist = IDT_data['R']
    #Thetalist = IDT_data['Theta']
    Rline_all = IDT_data['Rline']
    Phi_0 = Electrode_pattern(electrode_type,phi_0)
    N = Phi_0[0].size+Phi_0[1].size #number of electrodes over lambda
    IDT = [[],[]]
    for polarity in range(len(Rline_all)):
        for Rline in Rline_all[polarity]:
            #R = Rlist[polarity][i]
            #Theta = Thetalist[polarity][i]
            #xy = pol2cart(Theta,R)
            #Rline = shapely_geom.asLineString(xy)
            IDT[polarity].append(Rline.buffer(lambda_approx/(4*N))) #for a coverage of 50%
    IDT_data['IDT']=IDT
    if plot:
        #fig = plt.figure(1, figsize=SIZE, dpi=90)
        ax = plt.subplot(111)
        for polygon in IDT[0]:
            patch1 = PolygonPatch(polygon, fc=RED, ec=RED, alpha=0.5, zorder=2)
            ax.add_patch(patch1)
        for polygon in IDT[1]:
            patch2 = PolygonPatch(polygon, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2)
            ax.add_patch(patch2)    
        ax.set_title('a) dilation, cap_style=3')
        xrange = [1.1*min(xy[:,0]), 1.1*max(xy[:,0])]
        yrange = [1.1*min(xy[:,1]), 1.1*max(xy[:,1])]
        ax.axis([xrange[0],xrange[1], yrange[0], yrange[1]])
        ax.set_aspect(1)
        plt.show()        

def Draw_IDT(IDT_data):
    IDT = IDT_data['IDT']
    ax = plt.subplot(111)
    branch0 = IDT_data['electrodes']['elect0']
    branch1 = IDT_data['electrodes']['elect1']
    removal0 = IDT_data['removal']['removal0']
    removal1 = IDT_data['removal']['removal1']
    reticule = IDT_data['reticule']
    removal_center = IDT_data['removal']['removal_center']
    for polygon in IDT[0]:
        branch0 = branch0.union(polygon)       
    branch0 = branch0.difference(removal1)
    if 'Hole_elect0' in IDT_data['removal'].keys():
        branch0 = branch0.difference(IDT_data['removal']['Hole_elect0'])   
    branch0 = branch0.difference(removal_center)
    if type(branch0) is shapely_geom.polygon.Polygon:
        patch0 = PolygonPatch(branch0, fc=RED, ec=RED, alpha=0.5, zorder=2)
        ax.add_patch(patch0)
    else:
        tkMessageBox.showwarning(title = 'Warning',message = 'MultiPolygon detected')    
        for polygon in branch0.geoms:
            patch0 = PolygonPatch(polygon, fc=RED, ec=RED, alpha=0.5, zorder=2)
            ax.add_patch(patch0)
    
        
    for polygon in IDT[1]:
        branch1 = branch1.union(polygon)
    branch1 = branch1.difference(removal0)
    if 'Hole_elect1' in IDT_data['removal'].keys():
        branch1 = branch1.difference(IDT_data['removal']['Hole_elect1'])
    branch1 = branch1.difference(removal_center)
    if type(branch1) is shapely_geom.polygon.Polygon:
        patch1 = PolygonPatch(branch1, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2)
        ax.add_patch(patch1)
    else:
        tkMessageBox.showwarning(title = 'Warning',message = 'MultiPolygon detected')    
        for polygon in branch1.geoms:
            patch1 = PolygonPatch(polygon, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2)
            ax.add_patch(patch1)
     
    
    for polygon in reticule:
        patch_ret = PolygonPatch(polygon, fc=GREY, ec=GREY, alpha=0.5, zorder=3)
        ax.add_patch(patch_ret) 
    
    ax.set_title('IDT geometry')
    bounds0 = branch0.bounds
    bounds1 = branch1.bounds
    xrange = [1.1*min([bounds0[0],bounds1[0]]), 1.1*max([bounds0[2],bounds1[2]])]
    yrange = [1.1*min([bounds0[1],bounds1[1]]), 1.1*max([bounds0[3],bounds1[3]])]
    ax.axis([xrange[0],xrange[1], yrange[0], yrange[1]])
    ax.set_aspect(1)
    plt.show()
    #IDT_data['branch0']=branch0
    #IDT_data['branch1']=branch1

def Polygon2gds(IDT_data):
    layout = gds.core.Layout('LIBRARY')
    cell = gds.core.Cell('Main')
    points =  1e6*numpy.array(IDT_data['branch0'].exterior.coords) #best practice: dimensions in um
    boundaries = gds.core.Boundary(points)
    cell.add(boundaries)
    points = 1e6*numpy.array(IDT_data['branch1'].exterior.coords)
    boundaries = gds.core.Boundary(points)
    cell.add(boundaries)
    reticule = IDT_data['reticule']
    for polygon in reticule:
        points = 1e6*numpy.array(polygon.exterior.coords)
        boundaries = gds.core.Boundary(points)
        cell.add(boundaries)
    layout.add(cell)
    layout.save('output.gds')


def Import_reticule(IDT_data,reticule_size,reticule_filename):
    all_points = SVGT.read_path_from_svg(reticule_filename)
    reticule =[]
    for points in all_points:
        reticule.append(shapely_geom.Polygon(points.T))
    reticule = shapely_geom.MultiPolygon(reticule)
    outbox = reticule.bounds
    dx = outbox[2]-outbox[0]
    dy = outbox[3]-outbox[1]
    d = max([dx,dy])
    x0 = outbox[0]+dx/2
    y0 = outbox[1]+dy/2
    factor = 2*reticule_size/d
    reticule = shapely_affinity.translate(reticule, xoff=-x0, yoff=-y0) 
    reticule = shapely_affinity.scale(reticule, xfact = factor, yfact= factor, origin=(0,0,0))
    IDT_data['reticule'] = reticule

def IDT_union(IDT_data):
    IDT = IDT_data['IDT']
    branch0 = IDT_data['electrodes']['elect0']
    branch1 = IDT_data['electrodes']['elect1']
    removal0 = IDT_data['removal']['removal0']
    removal1 = IDT_data['removal']['removal1']
    #reticule = IDT_data['reticule']
    removal_center = IDT_data['removal']['removal_center']
    for polygon in IDT[0]:
        branch0 = branch0.union(polygon)       
    branch0 = branch0.difference(removal1)
    if 'Hole_elect0' in IDT_data['removal'].keys():
        branch0 = branch0.difference(IDT_data['removal']['Hole_elect0'])
    branch0 = branch0.difference(removal_center)
        
    for polygon in IDT[1]:
        branch1 = branch1.union(polygon)
    branch1 = branch1.difference(removal0)
    if 'Hole_elect1' in IDT_data['removal'].keys():
        branch1 = branch1.difference(IDT_data['removal']['Hole_elect1'])
    branch1 = branch1.difference(removal_center)
    
    IDT_data['branch0']=branch0
    IDT_data['branch1']=branch1
    IDT_merged = IDT[0]+IDT[1]
    IDT_bounds = shapely_geom.MultiPolygon(IDT_merged)
    IDT_bounds = IDT_bounds.convex_hull
    IDT_data['convex_hull'] = IDT_bounds
    
    
def IDT2svg(IDT_data,svg_dwg):
    IDT_bounds = IDT_data['convex_hull']
    branch0 = IDT_data['electrodes']['elect0']
    branch1 = IDT_data['electrodes']['elect1']
    Array = numpy.array(IDT_bounds.exterior.coords)
    svg_dwg  = SVGT.Array2svg(svg_dwg,1e5*Array,'simplify',stroke2 = 'green',fill2 = 'green')
    ArrayE0 = numpy.array(branch0.exterior.coords)
    ArrayE1 = numpy.array(branch1.exterior.coords)
    svg_dwg  = SVGT.Array2svg(svg_dwg,1e5*ArrayE0,'original',stroke2 = 'red',fill2 = 'red')
    svg_dwg  = SVGT.Array2svg(svg_dwg,1e5*ArrayE1,'original',stroke2 = 'blue',fill2 = 'blue')
    
def translate_IDT(IDT_data,x0,y0):
    #useful for the SVG drawing
    IDT_data['convex_hull'] = shapely_affinity.translate(IDT_data['convex_hull'],xoff = x0,yoff = y0)
    IDT_data['electrodes']['elect0'] = shapely_affinity.translate(IDT_data['electrodes']['elect0'],xoff = x0,yoff = y0)
    IDT_data['electrodes']['elect1'] = shapely_affinity.translate(IDT_data['electrodes']['elect1'],xoff = x0,yoff = y0)
    #useful for the union before gds drawing
    IDT = IDT_data['IDT'] 
    for ind in range(len(IDT[0])):
        IDT[0][ind] = shapely_affinity.translate(IDT[0][ind],xoff = x0,yoff = y0)
    for ind in range(len(IDT[1])):
        IDT[1][ind] = shapely_affinity.translate(IDT[1][ind],xoff = x0,yoff = y0)
    IDT_data['reticule'] = shapely_affinity.translate(IDT_data['reticule'],xoff = x0,yoff = y0)
    for key in IDT_data['removal'].keys():
        IDT_data['removal'][key] = shapely_affinity.translate(IDT_data['removal'][key],xoff = x0,yoff = y0)
    #IDT_data['removal']['removal1'] = shapely_affinity.translate(IDT_data['removal']['removal1'],xoff = x0,yoff = y0)
    
def Import_SVG_route(IDT_data):
    route_filename = IDT_data['svg_route']
    all_points = SVGT.read_path_from_svg(route_filename)
    route =[]
    for points in all_points:
        route.append(shapely_geom.Polygon(points.T))
    route = shapely_geom.MultiPolygon(route)
    #outbox = route.bounds
    #dx = outbox[2]-outbox[0]
    #dy = outbox[3]-outbox[1]
    #x0 = outbox[0]+dx/2
    #y0 = outbox[1]+dy/2
    x0svg = 4000
    y0svg = 4000
    factor = 1e-5;
    route = shapely_affinity.translate(route, xoff=-x0svg, yoff=-y0svg) 
    route = shapely_affinity.scale(route, xfact = factor, yfact= factor, origin=(0,0,0))
    IDT_data['route'] = route   