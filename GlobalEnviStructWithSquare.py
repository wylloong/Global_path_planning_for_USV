"""
Copyright(c) 2017, waylon
All rights reserved.
Distributed under the BSD license.
"""

# -*- coding: utf-8 -*-

# 全局环境建模

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import sys
from pylab import figure
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from shapely.geometry import Polygon  # 比较多边形交叉
import collections
import heapq  # heapq是一种子节点和父节点排序的树形数据结构
from math import isnan
from collections import Iterable
import json
import copy

# 类似于struct
Point = collections.namedtuple("Point", ["x", "y"])
OffsetCoord = collections.namedtuple("OffsetCoord", ["row", "col"])
SquareProperty = collections.namedtuple("SquareProperty",
                                        ["offsetCoord", "centerPoint", "cornerPoints", "weight", "squareSize",
                                         "isNavigonal"])
square_directions = [OffsetCoord(1, 1), OffsetCoord(1, -1), OffsetCoord(-1, -1), OffsetCoord(-1, 1)]
square_Neighbordirections = [OffsetCoord(0, -1), OffsetCoord(0, 1), OffsetCoord(-1, 0), OffsetCoord(1, 0),
                             OffsetCoord(1, 1), OffsetCoord(1, -1), OffsetCoord(-1, -1), OffsetCoord(-1, 1)]
square_directions = [OffsetCoord(1, 1), OffsetCoord(1, -1), OffsetCoord(-1, -1), OffsetCoord(-1, 1)]
square_Neighbordirections = [OffsetCoord(0, -1), OffsetCoord(0, 1), OffsetCoord(-1, 0), OffsetCoord(1, 0),
                             OffsetCoord(1, 1), OffsetCoord(1, -1), OffsetCoord(-1, -1), OffsetCoord(-1, 1)]


# 得到某区域正方形网格的所有属性
def SquarePro(ullong, ullati, squaresize, squareColumn, squareRow):
    savesquareprolist = []
    rowCount = 0
    while rowCount < squareRow:
        centerY = ullati - (rowCount + 0.5) * squaresize  # 中心点纬度
        columnCount = 0
        while columnCount < squareColumn:
            centerX = ullong + (columnCount + 0.5) * squaresize  # 中心点经度
            cornerlist = GetSquareCorners(centerX, centerY, squaresize)
            squarepro = SquareProperty(OffsetCoord(rowCount, columnCount), Point(centerX, centerY), cornerlist, 1,
                                       squaresize, True)
            savesquareprolist.append(squarepro)
            columnCount = columnCount + 1
        rowCount = rowCount + 1
    return savesquareprolist

# 得到某方向的一个点
def GetCorner(centerX, centerY, squaresize, cornerCount):
    cornerX = centerX + squaresize / 2 * square_directions[cornerCount].row
    cornerY = centerY + squaresize / 2 * square_directions[cornerCount].col
    return Point(cornerX, cornerY)

def getCenterPoint(ullong, ullati, squaresize, offsetcoord):
    rowCount = offsetcoord.row
    columnCount = offsetcoord.col
    centerY = ullati - (rowCount + 0.5) * squaresize  # 中心点纬度
    centerX = ullong + (columnCount + 0.5) * squaresize  # 中心点经度
    return Point(centerX, centerY)

# 得到一个正方形四个点的坐标
def GetSquareCorners(centerX, centerY, squaresize):
    cornerSum = 4
    cornerCount = 0
    cornerPoints = []
    while cornerCount < cornerSum:
        Points = GetCorner(centerX, centerY, squaresize, cornerCount)
        cornerPoints.append(Points)
        cornerCount = cornerCount + 1
    return cornerPoints

def intersection_recognition(singleCornerslist, lons, lats):
    p1 = Polygon([(singleCornerslist[0].x, singleCornerslist[0].y), (singleCornerslist[1].x, singleCornerslist[1].y),
                  (singleCornerslist[2].x, singleCornerslist[2].y), (singleCornerslist[3].x, singleCornerslist[3].y)])
    i = 0
    polygonlist = []
    while i < len(lats):
        t = (lons[i], lats[i])
        polygonlist.append(t)
        i = i + 1
    p2 = Polygon(polygonlist)
    return p1.intersects(p2)

def getweight(squareCoord, DisNavigonal):
    results = [square_neighbor(squareCoord, 0), square_neighbor(squareCoord, 1), square_neighbor(squareCoord, 2), \
               square_neighbor(squareCoord, 3), square_neighbor(squareCoord, 4), square_neighbor(squareCoord, 5), \
               square_neighbor(squareCoord, 6), square_neighbor(squareCoord, 7)]
    neighbourDisnaCount = 0;
    for x in results:
        if (x in DisNavigonal):
            neighbourDisnaCount += 1
    if (neighbourDisnaCount > 0):
        return 0.5 * (2 ** (neighbourDisnaCount - 1))
    else:
        return 0

def square_neighbor(squareoffset, direction):
    # print(hex_add(hexcube, hex_direction(direction)))
    return square_add(squareoffset, square_Neighbordirections[direction])


def square_add(a, b):
    return OffsetCoord(a.row + b.row, a.col + b.col)

squaresize = 0.002*1.6118548977
# 创建地图
f = figure(figsize=(15, 15))
ax = plt.subplot(111)
# plt.ion()  #interactive mode on
ullong=112.50
ullati=21.98
#下右
drlong=113
drlati=21.5
m = Basemap(projection='merc', llcrnrlat=drlati, urcrnrlat=ullati, \
            llcrnrlon=(ullong + 0.002), urcrnrlon=(drlong - 0.002), lat_ts=1, resolution='i')
m.drawcountries(linewidth=0.1, color='r')
# 正方形行数和列数
# print(drlong-ullong)
squareColumn = (drlong - ullong) // squaresize + 1  # 列数
squareRow = (ullati - drlati) // squaresize + 1  # 行数
print("网格化后矩阵维数 %d * %d" % (squareRow, squareColumn))
# 获得指定大小的正方形网格

squareprolist = SquarePro(ullong, ullati, squaresize, squareColumn, squareRow)

try:
    # 加载XML文件（2种方法,一是加载指定字符串，二是加载指定文件）
    tree = ET.parse("E:\\2017629_9115_ENCResolution.xml")  # 打开xml文档
    # tree= ElementTree.fromstring(text)
    root = tree.getroot()  # 获得root节点
except Exception, e:
    print "Error:cannot parse file"
    sys.exit(1)
for layerInfo in root.findall("LayerInfo"):  # 找到root节点下的所有xx节点
    for layer in layerInfo.findall("Layer"):  # 找到root节点下的所有xx节点
        # print(layer.get('id'))
        # print(layer.find('layerName').text)
        for feature in layer.find("Features"):
            feaid = feature.get('id')  # 子节点下属性name的值
            # print(feaid)
            elevation = feature.find('valdco').text
            # print(elevation)
            for wayPoints in feature.findall("wayPoints"):
                lats = []
                lons = []
                # high=[]
                for waypoint in wayPoints.findall("waypoint"):
                    waypoint_ID = waypoint.find('id').text
                    longitude = waypoint.find('lon').text
                    latitude = waypoint.find('lat').text
                    lons.append(float(longitude))
                    lats.append(float(latitude))
                # 两个多边形交叉判断
                if (len(lats) > 2):
                    compareCount = 0
                    while compareCount < len(squareprolist):
                        if (squareprolist[compareCount].isNavigonal == True):  # 可航时才检查
                            singleCornerslist = squareprolist[compareCount].cornerPoints
                            ans = intersection_recognition(singleCornerslist, lons, lats)  # 两个多边形交叉返回true，否则false
                            if (ans):
                                # 标记不可航行区域
                                temp = squareprolist[int(compareCount)]._replace(isNavigonal=False)
                                squareprolist[int(compareCount)] = temp
                        compareCount = compareCount + 1

# 设置权重
i = 0
Navigonal = []
DisNavigonal = []

while i < len(squareprolist):
    Navigonal.append((squareprolist[i].offsetCoord))
    if (squareprolist[i].isNavigonal == False):
        DisNavigonal.append((squareprolist[i].offsetCoord))
    i = i + 1
# 权重
Naviweight = {}  # 字典
weightIndex = 0
while weightIndex < len(squareprolist):
    weightoffset = 1+getweight(squareprolist[weightIndex].offsetCoord, DisNavigonal)
    temp = squareprolist[int(weightIndex)]._replace(weight=weightoffset)
    squareprolist[int(weightIndex)] = temp
    weightIndex += 1

# 环境建模完成，需要保存到txt中
# 写入到txt中
with open("E:\\squareEnvi5.txt", "w") as f:
    for square in squareprolist:
        f.write(json.dumps(square._asdict()))
        f.write("\n")
