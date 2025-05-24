#CENG-487 Assignment-4
# 320201105 HAMMET POLAT
#05/2025

from vec3d import Vec3d
from mat3d import Mat3d
from geometry import Geo

import numpy as np
from math import pi,sin,cos,sqrt,acos


class Object3D:
    def __init__(self, vertices):
        # Original model vertices (list of Vec3d)
        self.original_vertices = vertices

        self.transformed_vertices = vertices

        self.translation = Mat3d()
        self.rotation_x = Mat3d()
        self.rotation_y = Mat3d()
        self.rotation_z = Mat3d()
        self.scaling = Mat3d()
        self.transform_order = "SRT"#RST
        self.a=0
        self.lock=False

        self.face=None
        self.vertex=None

        self.lines=False

    def locked(self):
        self.lock=False
    def unlocked(self):
        self.lock=True
    def get_a(self):
        return self.a
    def inc_a(self):
        if self.lock:
            self.a=self.a+1
    def dec_a(self):
        if self.lock and self.a>0:
            self.a=self.a-1

    def changeLines(self):
        self.lines=not self.lines

    def getLines(self):
        return self.lines

    def setTransformOrder(self, order_str):#default transform order   !!!!!TSR çalışmıyor
        self.transform_order = order_str

    def setTranslation(self, tx, ty, tz):
        self.translation = Mat3d().translationMatrix(tx, ty, tz)

    def setScaling(self, sx, sy, sz):
        self.scaling = Mat3d().scalingMatrix(sx, sy, sz)

    def setRotation(self, angle_x, angle_y, angle_z):
        self.rotation_x = Mat3d().rotationMatrix(angle_x, 'x')
        self.rotation_y = Mat3d().rotationMatrix(angle_y, 'y')
        self.rotation_z = Mat3d().rotationMatrix(angle_z, 'z')


    def applyTransformations(self):
        components = {
            'T': self.translation,
            'R': self.rotation_z.matrixMult(self.rotation_y).matrixMult(self.rotation_x),
            'S': self.scaling
        }
        full_matrix=Mat3d()
        for code in self.transform_order:
            full_matrix = full_matrix.matrixMult(components[code])#components[code].matrixMult(full_matrix)

        self.transformed_vertices = [
            full_matrix.matrixVecMult(v) for v in self.transformed_vertices
        ]

    def getVertices(self):
        return self.transformed_vertices

    def resetFullToIdentity(self):
        # Reset full_matrix to identity matrix
        self.full_Matrix = Mat3d()  # Mat3d() creates an identity matrix
    def resetOriginal(self):
        # Reset full_matrix to identity matrix
        self.transformed_vertices = self.original_vertices  # Mat3d() creates an identity matrix

    def makeCube(self):
        a = self.get_a()
        self.cubeArr=Geo.makeCube(a)
        return self.cubeArr


    def makePlane(self):
        a = self.get_a()
        self.planeArr = Geo.makePlane(a)
        return self.planeArr


    def makePyramid(self):
        a = self.get_a()
        self.pyrArr = Geo.makePyramid(a)
        return self.pyrArr


    def makeSphere(self):
        a=self.get_a()
        if a<0:
            self.inc_a()
            a=a+1
        a=a+1
        self.sphereArr = Geo.makeSphere(a)
        return self.sphereArr


    def makeSphere2(self):
        a = self.get_a()
        self.sphereArr2 = Geo.makeSphere2(a)
        return self.sphereArr2


    def makeTriangle(self):
        a = self.get_a()
        self.triArr = Geo.makeTriangle(a)
        return self.triArr

    def makeParse(self):
        vertex = self.getVertexes()
        faces = self.getFaces()
        subFaces=[]
        for face in faces:
            f0 = vertex[face[0]]
            f1 = vertex[face[1]]
            f2 = vertex[face[2]]
            f3 = vertex[face[3]]
            subFaces.append([f0,f1,f2,f3])

        a = self.get_a()
        # we take the orignal faces and their coordinates as a list
        self.catmullClarkParse = Geo.catmullClark(subFaces,a)

        return self.catmullClarkParse



    def draw(self,arr):
        line=self.getLines()
        Geo.drawObject(arr,line)

    def setVertexes(self,x):
        self.vertex=x

    def getVertexes(self):
        return self.vertex


    def setFaces(self,x):
        self.face=x

    def getFaces(self):
        return  self.face


    @staticmethod
    def resetCam(x,y,z):
        Geo.resetGeoCam(x,y,z)
    @staticmethod
    def rotateAround(x,y):
        Geo.rotateAround(x,y)
    @staticmethod
    def dragZ(z):
        Geo.dragZ(z)
    @staticmethod
    def changeLookInto(x,y,z):
        Geo.changeLookInto(x,y,z)
