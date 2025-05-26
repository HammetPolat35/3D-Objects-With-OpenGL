#HAMMET POLAT
#05/2025

from vec3d import Vec3d
from math import pi,sin,cos,sqrt,acos

import numpy as np


class Mat3d:
    #took it from repo
    @staticmethod
    def create(arg):
        if len(arg) != 16:
            raise Exception("array size must be 16!")

        return Mat3d(arg)

    def __init__(self, rows=None):
        if rows is None:
            self.na = np.identity(4)
        else:
            self.na = np.array(rows)
            self.na = np.reshape(self.na, (4, 4))
            self.na = np.transpose(self.na)

    def fromNumpy(self, a):
        self.na = a

    def __add__(self, other):
        return Mat3d(self.na + other.na)

    def consMult(self, other):
        return Mat3d(self.na*other)

    def transpose(self):
        return Mat3d(np.transpose(self.na))


    def matrixMult(self, other):
        result = np.dot(self.na, other.na)
        return Mat3d(result)

    def matrixVecMult(self, vec):
        # Matrix-vector multiplication: transforms a vector using this matrix
        # Convert Vec3d to a 4D array [x, y, z, w]
        point = np.array([vec.x(), vec.y(), vec.z(), vec.w()])

        # Perform matrix multiplication (4x4 * 4x1)
        result = np.dot(self.na, point)


        return Vec3d(result)  # Only return the 3D part (x, y, z)

    def translationMatrix(self, tx, ty, tz):
        # Returns the translation matrix
        return Mat3d([
            [1, 0, 0, tx],
            [0, 1, 0, ty],
            [0, 0, 1, tz],
            [0, 0, 0, 1]
        ])

    def applyTranslation(self, vec, tx, ty, tz):
        # Apply translation transformation to a vector (Vec3d object)
        translation_matrix = self.translationMatrix(tx, ty, tz)
        return translation_matrix.matrixVecMult(vec)

    def scalingMatrix(self, sx, sy, sz):
        # Returns the scaling matrix
        return Mat3d([
            [sx, 0, 0, 0],
            [0, sy, 0, 0],
            [0, 0, sz, 0],
            [0, 0, 0, 1]
        ])

    def applyScaling(self, vec, sx, sy, sz):
        # Apply scaling transformation to a vector (Vec3d object)
        scaling_matrix = self.scalingMatrix(sx, sy, sz)
        return scaling_matrix.matrixVecMult(vec)

    def rotationMatrix(self, a, axis):
        # Generates a rotation matrix for a given axis and angle in radians
        a = a * pi / 180
        if axis == "x":
            rotation_matrix = Mat3d([
                [1, 0, 0, 0],
                [0, cos(a), -sin(a), 0],
                [0, sin(a), cos(a), 0],
                [0, 0, 0, 1]
            ])
        elif axis == "y":
            rotation_matrix = Mat3d([
                [cos(a), 0, sin(a), 0],
                [0, 1, 0, 0],
                [-sin(a), 0, cos(a), 0],
                [0, 0, 0, 1]
            ])
        elif axis == "z":
            rotation_matrix = Mat3d([
                [cos(a), -sin(a), 0, 0],
                [sin(a), cos(a), 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]
            ])
        else:
            raise ValueError("Invalid axis. Choose 'x', 'y', or 'z'.")
        return rotation_matrix

    def applyRotation(self, vec, angle, axis):
        # Apply rotation to a vector using the specified angle and axis (angle in degrees)
        #angle=angle*pi/180   # Convert to radians
        rotation_matrix = self.rotationMatrix(angle, axis)
        return rotation_matrix.matrixVecMult(vec)


'''def matrixMult(self,other):
        for i in range(len(self.length)):
            for j in range(len(other[i])):
                for k in range(len(other)):
                    self.res[i][j] += self[i][k] * other[k][j]
        return self.res'''