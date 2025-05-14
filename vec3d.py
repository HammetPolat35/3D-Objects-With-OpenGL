#CENG-487 Assignment-3
# 320201105 HAMMET POLAT
#05/2025

import numpy as np
from math import acos

class Vec3d:
    def __init__(self, a):
        a = np.array(a)
        if a.shape[0] == 3:
            a = np.concatenate([a, [1]])
        self.na = a

    def from4d(self, x, y, z, w=1):
        return ([x, y, z, w])

    def x(self):
        return self.na[0]

    def y(self):
        return self.na[1]

    def z(self):
        return self.na[2]

    def w(self):
        return self.na[3]

    def __add__(self, other):
        if self.w() == other.w():
            result = self.na[0:3] + other.na[0:3]
            return Vec3d(np.concatenate([result, [self.w()]]))
        raise ValueError("Error")

    def __sub__(self, other):
        if self.w() == other.w():
            result = self.na[0:3] - other.na[0:3]
            return Vec3d(np.concatenate([result, [self.w()]]))
        raise ValueError("Error")

    def linearDep(self, other):
        if self.w() != other.w():
            raise ValueError("Error")
        ratios = [self.x() / other.x(), self.y() / other.y(), self.z() / other.z()]
        return np.allclose(ratios, ratios[0])

    def scale(self, scalar):
        return Vec3d(np.concatenate([self.na[0:3] * scalar, [self.w()]]))

    def dot(self, other):
        if self.w() == other.w():
            return np.dot(self.na[0:3], other.na[0:3]) + self.w() * other.w()
        raise ValueError("Error")

    def cross(self, other):
        cross = np.cross(self.na[0:3], other.na[0:3])
        return Vec3d(np.concatenate([cross, [0]]))  # w=0 for cross product

    def __length__(self):
        return np.linalg.norm(self.na[0:3])

    def angle(self, other):
        dot_product = self.dot(other)
        length_self = self.__length__()
        length_other = other.__length__()
        return acos(dot_product / (length_self * length_other))

    def projection(self, other):
        a = self.angle(other)
        b = self.__length__() * np.cos(a) / other.__length__()
        return other.scale(b)

    def normalize(self):
        length = self.__length__()
        if length == 0:
            raise ValueError("Cannot normalize a zero vector")
        return Vec3d(np.concatenate([self.na[0:3] / length, [self.w()]]))

    @staticmethod
    def interpolate_bilinear(v0, v1, v2, v3, s, t):
        x_top = v0.x() + (v1.x() - v0.x()) * s
        y_top = v0.y() + (v1.y() - v0.y()) * s
        z_top = v0.z() + (v1.z() - v0.z()) * s

        x_bottom = v3.x() + (v2.x() - v3.x()) * s
        y_bottom = v3.y() + (v2.y() - v3.y()) * s
        z_bottom = v3.z() + (v2.z() - v3.z()) * s

        x = x_top + (x_bottom - x_top) * t
        y = y_top + (y_bottom - y_top) * t
        z = z_top + (z_bottom - z_top) * t

        return Vec3d([x, y, z])













