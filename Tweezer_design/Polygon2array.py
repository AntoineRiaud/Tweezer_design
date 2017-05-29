# -*- coding: utf-8 -*-
"""
Created on Mon Aug 01 17:31:49 2016

@author: Antoine
"""

from numpy import array
import shapely.geometry as shapely_geom

def Polygon2Array(MyObj):
    if type(MyObj) is dict:
        for key in MyObj.keys():
            try: #here we check if the key is an number, in which case it is not possible to use it to export matlab structures
                float(key) 
                new_key = 'instance_'+key
                MyObj[new_key] = Polygon2Array(MyObj[key])
                MyObj.pop(key,None)
            except ValueError:   #the key is not a number
                MyObj[key] = Polygon2Array(MyObj[key])
    elif type(MyObj) is list:
        for i in range(len(MyObj)):
            MyObj[i] = Polygon2Array(MyObj[i])                            
    elif type(MyObj) is shapely_geom.polygon.Polygon:
        MyObj = array(MyObj.exterior)
    elif type(MyObj) is shapely_geom.multipolygon.MultiPolygon:
        MyObj = Polygon2Array(MyObj.geoms)
    elif type(MyObj) is shapely_geom.base.GeometrySequence:
        MyOldObj = MyObj
        MyObj = []
        for i in range(len(MyOldObj)):
            MyObj.append(Polygon2Array(MyOldObj[i]))
    return MyObj