# -*- coding: utf-8 -*-
"""
Copyright(c) 2017, waylon
All rights reserved.
Distributed under the BSD license.
"""

from shapely.geometry import Polygon,Point  # 比较多边形交叉

def point_in_polygon(point_x, point_y, polygonlist_x,polygonlist_y):
    #如果出问题，注意坐标转换
    point = Point(point_x,point_y)
    i = 0
    polygonlist = []
    while i < len(polygonlist_x):
        t = (polygonlist_x[i], polygonlist_y[i])
        polygonlist.append(t)
        i = i + 1
    poly = Polygon(polygonlist)
    return poly.contains(point)

if __name__=='__main__':
    # 测试点包含在多边形中
    polygonlist_x=[]
    polygonlist_x.append(25.774252)
    polygonlist_x.append(18.466465)
    polygonlist_x.append(32.321384)
    polygonlist_y = []
    polygonlist_y.append(-80.190262)
    polygonlist_y.append(-66.118292)
    polygonlist_y.append(-64.75737)
    result=point_in_polygon(27.254629577800088, -76.728515625,polygonlist_x,polygonlist_y)
    print(result)
    result=point_in_polygon(27.254629577800088, -74.928515625,polygonlist_x,polygonlist_y)
    print(result)