# -*- coding: utf-8 -*-
"""
Copyright(c) 2017, waylon
All rights reserved.
Distributed under the BSD license.
"""

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
import GeometryBase

# 类似于struct
Point = collections.namedtuple("Point", ["x", "y"])
OffsetCoord = collections.namedtuple("OffsetCoord", ["row", "col"])
SquareProperty = collections.namedtuple("SquareProperty",
                                        ["offsetCoord", "centerPoint", "cornerPoints", "weight", "squareSize",
                                         "isNavigonal"])
square_directions = [OffsetCoord(1, 1), OffsetCoord(1, -1), OffsetCoord(-1, -1), OffsetCoord(-1, 1)]
square_Neighbordirections = [OffsetCoord(0, -1), OffsetCoord(0, 1), OffsetCoord(-1, 0), OffsetCoord(1, 0),
                             OffsetCoord(1, 1), OffsetCoord(1, -1), OffsetCoord(-1, -1), OffsetCoord(-1, 1)]

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


# A Star Def
def square_add(a, b):
    return OffsetCoord(a.row + b.row, a.col + b.col)


def square_neighbor(squareoffset, direction):
    return square_add(squareoffset, square_Neighbordirections[direction])


class SquareGrid(object):
    def __init__(self, polygonlists):
        self.polygonlist = polygonlists
        self.walls = []

    def in_bounds(self, squareCoord):
        return squareCoord in self.polygonlist

    def passable(self, squareCoord):
        return squareCoord not in self.walls

    def neighbors(self, squareCoord):
        results = [square_neighbor(squareCoord, 0), square_neighbor(squareCoord, 1), square_neighbor(squareCoord, 2), \
                   square_neighbor(squareCoord, 3), square_neighbor(squareCoord, 4), square_neighbor(squareCoord, 5), \
                   square_neighbor(squareCoord, 6), square_neighbor(squareCoord, 7)]
        results = filter(self.in_bounds, results)  # filter()函数用于过滤序列
        results = filter(self.passable, results)  # 可行区域
        return results


class GridWithWeights(SquareGrid):
    def __init__(self, polygonlist):
        super(GridWithWeights, self).__init__(polygonlist)  # super() 子类调用父类成员
        self.weights = {}

    def cost(self, start, goal, from_node, to_node):
        return heuristic(from_node, to_node)  # 如果不存在，则返回1

class PriorityQueue:
    def __init__(self):
        self.elements = []  # 创建了一个空堆

    def empty(self):
        return len(self.elements) == 0

    def put(self, item, priority):
        heapq.heappush(self.elements, (priority, item))  # 往堆中插入一条新的值

    def get(self):
        return heapq.heappop(self.elements)[1]  # 从堆中弹出最小值

def heuristic(a, b):
    (x1, y1) = a
    (x2, y2) = b
    return np.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))

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

def a_star_search(graph, start, goal):
    frontier = PriorityQueue()  # 新建堆数据结构，类似于OPEN表
    frontier.put(start, 0)  # 初始节点放入open表
    came_from = {}  # 字典数据类型
    cost_so_far = {}  # close表
    came_from[start] = None
    cost_so_far[start] = 0  # open表和closed表综合，代价函数g(x)
    IteraCount = 0
    while not frontier.empty():  # 如果open表为空，则搜索失败，退出
        current = frontier.get()  # 从堆（open表）中弹出最小值，放入closed表
        if current == goal:  # 考察节点n是否为目标节点，若是，则求得了问题的解
            break
        for next in graph.neighbors(current):  # 扩展节点current，在上下左右生成四个可行节点
            # 对节点n的每个子节点i计算估价函数f(i)
            new_cost = cost_so_far[current] + graph.cost(start, goal, current, next)  # 到达next的代价g(x)
            if next not in cost_so_far or new_cost < cost_so_far[next]:
                if next not in cost_so_far:
                    IteraCount += 1
                cost_so_far[next] = new_cost  # 代价函数g(x)
                came_from[next] = current  # next的父节点是current
                # dx1 = next.row - goal.row
                # dy1 = next.col - goal.col
                # dx2 = start.row - goal.row
                # dy2 = start.col - goal.col
                # anglediff = np.abs(dx1 * dy2 - dx2 * dy1) / (
                # np.sqrt(dx1 * dx1 + dy1 * dy1) * np.sqrt(dx2 * dx2 + dy2 * dy2))
                #
                # if (isnan(anglediff)):
                #     print(anglediff)
                #     anglediff = 0;
                # wgh = 9 / (10 - anglediff)  # 效果比较好
                wgh=1
                priority = new_cost + heuristic(goal, next) * wgh  # 估价函数f(x)
                frontier.put(next, priority)  # 放入open表
    return came_from, cost_so_far, IteraCount

# A Star end Def

def dijkstra_search(graph, start, goal):
    frontier = PriorityQueue()
    frontier.put(start, 0)
    came_from = {}
    cost_so_far = {}
    came_from[start] = None
    cost_so_far[start] = 0

    while not frontier.empty():
        current = frontier.get()
        if current == goal:
            break

        for next in graph.neighbors(current):
            new_cost = cost_so_far[current] + graph.cost(start, goal, current, next)
            if next not in cost_so_far or new_cost < cost_so_far[next]:
                cost_so_far[next] = new_cost
                priority = new_cost
                frontier.put(next, priority)
                came_from[next] = current
    return came_from, cost_so_far

# 曲线优化Start
def generallinearequation(firstPnt, thirdPnt):
    # 判断第一点和第二点连线路过的点 Ax+By+c=0
    A = thirdPnt.row - firstPnt.row
    B = firstPnt.col - thirdPnt.col
    C = thirdPnt.col * firstPnt.row - firstPnt.col * thirdPnt.row
    return A, B, C

def curveoptimize(curveIndex, curvelist, DisNavigonal, distance):
    IsCandelete = True;  # 默认可以删除
    while (IsCandelete == True and len(curvelist) > curveIndex + 2):
        firstPnt = curvelist[curveIndex]
        thirdPnt = curvelist[curveIndex + 2]

        A, B, C = generallinearequation(firstPnt, thirdPnt)
        # 求平行于该直线的两条直线C1和C2的值
        lineDistance = np.sqrt(A ** 2 + B ** 2) * distance  # 平行直线的间距
        C1 = C + lineDistance
        C2 = C - lineDistance
        # XY的范围,从小到大
        xsort = sorted((firstPnt.col, thirdPnt.col))
        ysort = sorted((firstPnt.row, thirdPnt.row))
        # print(xsort,ysort)
        # 获得C1和C2夹层中的点坐标，用来后续检查
        if (B != 0):  # 不垂直于x轴
            for x in range(xsort[0], xsort[1] + 2):
                if (IsCandelete == False):
                    break
                y1 = (-C1 - A * x) / B
                y2 = (-C2 - A * x) / B
                ychecksort = sorted((y1, y2))
                # 20170401 x应在xsort范围内
                if (ychecksort[0] < (ysort[0] - 1)):
                    ychecksort[0] = ysort[0] - 1
                if (ychecksort[1] > (ysort[1] + 1)):
                    ychecksort[1] = ysort[1] + 1
                # 得到需要考察的点
                for y in range(int(ychecksort[0] + 0.45), int(ychecksort[1] + 0.55)):
                    # 存在不可航行点，退出
                    if (OffsetCoord(y, x) in DisNavigonal):
                        IsCandelete = False
                        break
        if (A != 0):  # 不垂直于y轴
            for y in range(ysort[0], ysort[1] + 2):
                if (IsCandelete == False):
                    break
                x1 = (-C1 - B * y) / A
                x2 = (-C2 - B * y) / A
                xchecksort = sorted((x1, x2))
                # 20170401 x应在xsort范围内
                if (xchecksort[0] < (xsort[0] - 1)):
                    xchecksort[0] = xsort[0] - 1
                if (xchecksort[1] > (xsort[1] + 1)):
                    xchecksort[1] = xsort[1] + 1
                # 得到需要考察的点
                for x in range(int(xchecksort[0] + 0.45), int(xchecksort[1] + 0.55)):
                    # 存在不可航行点，退出
                    if (OffsetCoord(y, x) in DisNavigonal):
                        IsCandelete = False
                        # print("XChange",firstPnt,thirdPnt,y,x,xchecksort)
                        # print(x,y)
                        break
        if (IsCandelete == True):
            # 可删除这个点
            curvelist.remove(curvelist[curveIndex + 1])
            IsCandelete = True
    return curvelist
# 曲线优化end

def getLineLength(firstwayPnt, secondwayPnt):
    x = [firstwayPnt.y, secondwayPnt.y]
    y = [firstwayPnt.x, secondwayPnt.x]
    # compute map projection coordinates for lat/lon grid.
    xlist, ylist = m(y, x)  # 将经纬度转换为图片所用坐标,米坐标
    return np.sqrt((xlist[0] - xlist[1]) * (xlist[0] - xlist[1]) + (ylist[0] - ylist[1]) * (ylist[0] - ylist[1]))  # 米坐标

# 计算点在多边形的序号
def get_polygonIndexfromPoint(point_x, point_y,squarelist):
    result =False
    polygonlist_x=[]
    polygonlist_y=[]
    for square in squarelist:
        polygonlist_x = []
        polygonlist_y = []
        for corner in square.cornerPoints:
            polygonlist_x.append(corner.x)
            polygonlist_y.append(corner.y)
        result=GeometryBase.point_in_polygon(point_x,point_y,polygonlist_x,polygonlist_y)
        if result==True:
            return square.offsetCoord

# 创建地图
f = figure(figsize=(15, 15))
ax = plt.subplot(111)
# square size
squaresize = 0.002
ullong=112.50
ullati=21.98
#下右
drlong=113
drlati=21.5
m = Basemap(projection='merc', llcrnrlat=drlati, urcrnrlat=ullati, \
            llcrnrlon=(ullong + 0.002), urcrnrlon=(drlong - 0.002), lat_ts=1, resolution='i')
m.drawcountries(linewidth=0.1, color='r')

#绘制从电子海图解析出的陆地等不可航区域
try:
    # 加载XML文件（2种方法,一是加载指定字符串，二是加载指定文件）
    tree=ET.parse("E:\\2017629_9115_ENCResolution.xml") #打开xml文档
    #tree= ElementTree.fromstring(text)
    root=tree.getroot()  #获得root节点
except Exception,e:
    print "Error:cannot parse file"
    sys.exit(1)
for layerInfo in root.findall("LayerInfo"):  #找到root节点下的所有xx节点
    for layer in layerInfo.findall("Layer"):  #找到root节点下的所有xx节点
        #print(layer.get('id'))
        #print(layer.find('layerName').text)
        for feature in layer.find("Features"):
            feaid=feature.get('id') #子节点下属性name的值
            elevation=feature.find('valdco').text
            for wayPoints in feature.findall("wayPoints"):
                lats=[]
                lons=[]
                for waypoint in wayPoints.findall("waypoint"):
                    waypoint_ID=waypoint.find('id').text
                    longitude=waypoint.find('lon').text
                    latitude=waypoint.find('lat').text
                    lons.append(float(longitude))
                    lats.append(float(latitude))
                if(len(lats)>2):
                    x, y = m(lons, lats)
                    shpsegs = []
                    shpsegs.append(zip(x,y))
                    lines = LineCollection(shpsegs,antialiaseds=(1,))
                    #lines.set_facecolors(cm.jet(np.random.rand(1))) #随机颜色
                    lines.set_facecolors(cm.jet(0.02))
                    lines.set_edgecolors('g')
                    lines.set_linewidth(0.3)
                    ax.add_collection(lines)

# 从Txt中读取环境建模保存的信息
squareprolist = []
with open("E:\\squareEnvi4.txt", "r") as r:
    filelist = r.readlines()
    for x in filelist:
        dic = json.loads(x)  # 输出dict类型
        p = SquareProperty(**dic)
        offsetCoord=OffsetCoord(p.offsetCoord[0],p.offsetCoord[1])
        centerPoint=Point(p.centerPoint[0],p.centerPoint[1])
        cornerPoints=[]
        point1=Point(p.cornerPoints[0][0],p.cornerPoints[0][1])
        cornerPoints.append(point1)
        point2=Point(p.cornerPoints[1][0],p.cornerPoints[1][1])
        cornerPoints.append(point2)
        point3=Point(p.cornerPoints[2][0],p.cornerPoints[2][1])
        cornerPoints.append(point3)
        point4=Point(p.cornerPoints[3][0],p.cornerPoints[3][1])
        cornerPoints.append(point4)
        isnavigonal=p.isNavigonal
        squaresize=p.squareSize
        weight=p.weight
        squarepro=SquareProperty(offsetCoord,centerPoint,cornerPoints,weight,squaresize,isnavigonal)
        squareprolist.append(squarepro)

# 输出航行区域与不可航行区域
i = 0
Navigonal = []
DisNavigonal = []
# 权重
Naviweight = {}  # 字典
while i < len(squareprolist):
    Navigonal.append((squareprolist[i].offsetCoord))
    Naviweight[squareprolist[i].offsetCoord] = squareprolist[i].weight
    if (squareprolist[i].isNavigonal == False):
        DisNavigonal.append((squareprolist[i].offsetCoord))
    i = i + 1

#绘制正方形网格以及可航区域与不可航区域
plotsquareCount=0
while plotsquareCount<len(squareprolist):
    squarecorner= squareprolist[plotsquareCount].cornerPoints
    shadowlons=[squarecorner[0].x,squarecorner[1].x,squarecorner[2].x,squarecorner[3].x]
    shadowlats=[squarecorner[0].y,squarecorner[1].y,squarecorner[2].y,squarecorner[3].y]
    x, y = m(shadowlons, shadowlats)
    shpsegs = []
    shpsegs.append(zip(x,y))
    lines = LineCollection(shpsegs,antialiaseds=(1,))
    #不可航行区域
    if(squareprolist[plotsquareCount].isNavigonal==False):
        lines.set_facecolors(cm.jet(0.1))
        lines.set_edgecolors('g')
        lines.set_linewidth(0.6)
        lines.set_alpha(0.6) #设置透明度
        ax.add_collection(lines)  #绘制不可行区域
    else:
        lines.set_facecolors(cm.jet(0.02))
        lines.set_edgecolors('b')
        lines.set_linewidth(1.2)
        #lines.set_alpha(0.1) #设置透明度

        # 设置颜色深度，权重越大深度越大
        weight = Naviweight.get(squareprolist[plotsquareCount].offsetCoord, 1)
        if weight < 1:
            weight = 1
        if weight > 10:
            weight = 10
        lines.set_alpha(0.1 * weight)  # 设置透明度
        ax.add_collection(lines)
    plotsquareCount=plotsquareCount+1
print("环境构建成功")



diagram4 = GridWithWeights(Navigonal)
# 加入不可航性区域
diagram4.walls = DisNavigonal
# 加入每一个点的权重
diagram4.weights = Naviweight
# 起始点

startPoint = get_polygonIndexfromPoint(112.663000,21.830000,squareprolist)
print(startPoint)
# 目标点
targetPoint = get_polygonIndexfromPoint(112.753098,21.664787,squareprolist)
print(targetPoint)
# came_from, cost_so_far = dijkstra_search(diagram4,startPoint,targetPoint)
came_from, cost_so_far, Iteracount = a_star_search(diagram4, startPoint, targetPoint)
print("遍历点数 %d" % (Iteracount))
print("A*算法寻得最短路径")
pathFinderPoints = []
# 回溯父节点
pathFinderPoints.append(targetPoint)
fatherNode = came_from[targetPoint]
while fatherNode != startPoint:
    pathFinderPoints.append(fatherNode)
    fatherNode = came_from[fatherNode]
    # print(fatherNode)
pathFinderPoints.append(startPoint)
print("优化前节点个数： %d" % len(pathFinderPoints))
# 在地图上绘制最短路径
lineCount = 0;
PlanedcenterPnt=[]
while lineCount < (len(pathFinderPoints) - 1):
    last = pathFinderPoints[lineCount]
    new = pathFinderPoints[lineCount + 1]
    # print(last)
    lastCenter = getCenterPoint(ullong, ullati, squaresize, last)  # Point
    PlanedcenterPnt.append(lastCenter)
    newCenter = getCenterPoint(ullong, ullati, squaresize, new)  # Point
    if lineCount == 0:
        y1 = [lastCenter.y]
        x1 = [lastCenter.x]
        x1, y1 = m(x1, y1)
        m.scatter(x1, y1, c='r', marker='o')
    if lineCount == len(pathFinderPoints) - 2:
        y1 = [newCenter.y]
        x1 = [newCenter.x]
        x1, y1 = m(x1, y1)
        m.plot(x1, y1, c='r', marker='o')
    x = [lastCenter.y, newCenter.y]
    y = [lastCenter.x, newCenter.x]
    xlist, ylist = m(y, x)  # 将经纬度转换为图片所用坐标
    #m.plot(xlist, ylist, color='b', linestyle=":")  # 首先画白线
    m.plot(xlist, ylist, color='r')  # 首先画白线
    lineCount = lineCount + 1
last = pathFinderPoints[lineCount]
lastCenter = getCenterPoint(ullong, ullati, squaresize, last)  # Point
PlanedcenterPnt.append(lastCenter)
# 规划路径曲线长度
wayCount = 0
planeddistanceall = 0
while wayCount < (len(PlanedcenterPnt) - 1):
    planeddistanceall = planeddistanceall + getLineLength(PlanedcenterPnt[wayCount], PlanedcenterPnt[wayCount + 1])  # 米坐标
    wayCount += 1
print(planeddistanceall)

# A*曲线优化为近似直线
curvelist = copy.deepcopy(pathFinderPoints)
curvelist.reverse()  # 逆序排列
curveCount = 0
distance = 1.5  # 设置距离，但是这里不太有效吧，仍然有问题
while curveCount < len(curvelist) - 2:
    # print("优化中")
    # 检测是否可以优化
    curvelist = curveoptimize(curveCount, curvelist, DisNavigonal, distance)
    curveCount += 1;
print("路径曲线优化结束")
# 曲线再次优化，算法存在问题，之前不能优化的，现在又可以优化了
curveCount = 0;
distance = 2
while curveCount < len(curvelist) - 2:
    # print("优化中")
    # 检测是否可以优化
    curvelist = curveoptimize(curveCount, curvelist, DisNavigonal, distance)
    curveCount += 1;
print("路径曲线再次优化结束")
# 得到节点中心点和个数
print(curvelist)
centerPnt = []
print("优化后节点个数： %d" % len(curvelist))
for node in curvelist:
    for waypnt in squareprolist:
        if waypnt.offsetCoord == node:
            centerPnt.append(waypnt.centerPoint)
            break
print(centerPnt)
# 规划路径曲线长度
wayCount = 0
distanceall = 0
while wayCount < (len(centerPnt) - 1):
    distanceall = distanceall + getLineLength(centerPnt[wayCount], centerPnt[wayCount + 1])  # 米坐标
    wayCount += 1
print(distanceall)
lineCount = 0
while lineCount < (len(curvelist) - 1):
    last = curvelist[lineCount]
    new = curvelist[lineCount + 1]
    lastCenter = getCenterPoint(ullong, ullati, squaresize, last)  # Point
    newCenter = getCenterPoint(ullong, ullati, squaresize, new)  # Point
    x = [lastCenter.y, newCenter.y]
    y = [lastCenter.x, newCenter.x]
    xlist, ylist = m(y, x)  # 将经纬度转换为图片所用坐标
    #m.plot(xlist, ylist, color='r')  # 首先画白线
    lineCount = lineCount + 1
plt.show()
print("界面显示成功")