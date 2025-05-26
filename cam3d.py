#HAMMET POLAT
#05/2025

from mat3d import Mat3d
from vec3d import Vec3d
import numpy as np
from math import pi,sin,cos,sqrt,acos

class Cam3D:
    def __init__(self):
        self.cam=Vec3d([0.0,0.0,0.0,0.0])
        self.target=Vec3d([0.0,0.0,-8.0,0.0])  #where are we looking
        self.mock_up = Vec3d([0.0, 1.0, 0.0, 0.0])

        self.look_dir = self.target.__sub__(self.cam).normalize()
        self.horizon = self.look_dir.cross(self.mock_up).normalize()
        self.up = self.horizon.cross(self.look_dir).normalize()

        self.radius = 8
        self.azimuth = 0
        self.elevation = 0

    def makeViewMatrix(self):
        self.arrView= [
            Vec3d([self.horizon.x(),self.up.x(),-self.look_dir.x(),0]),
            Vec3d([self.horizon.y(),self.up.y(),-self.look_dir.y(),0]),
            Vec3d([self.horizon.z(),self.up.z(),-self.look_dir.z(),0]),
            Vec3d([-self.horizon.dot(self.cam), -self.up.dot(self.cam),self.look_dir.dot(self.cam), 1])
        ]
        return self.arrView

    def returnViewMatrix(self):
        view_rows = self.makeViewMatrix()
        flattened = []
        for vec in view_rows:
            flattened.extend([vec.x(), vec.y(), vec.z(), vec.w()])

        self.viewMatrix = Mat3d(flattened)
        return self.viewMatrix

    def lookInto(self,x,y,z):
        self.target=Vec3d([x,y,z,0.0])
        self.calcAgain()

    def changeCamPos(self,x,y,z):
        self.radius = np.sqrt(x ** 2 + y ** 2 + (z+8) ** 2)
        self.cam=Vec3d([x,y,z,0])
        self.calcAgain()

    def rotateAround(self, dx, dy):#took help with this same logic with sphere
        # Update azimuth and elevation based on mouse movement
        self.azimuth += dx / 700 * (self.radius+4)/4
        self.elevation += dy / 500 * (self.radius+4)/4

        # Limit the vertical movement to avoid flipping over (e.g., between -90 and 90 degrees)
        self.elevation = np.clip(self.elevation, -pi / 2, np.pi / 2)

        # Convert spherical coordinates (radius, azimuth, elevation) to Cartesian coordinates
        x = self.radius * cos(self.elevation) * sin(self.azimuth)
        y = self.radius * sin(self.elevation)
        z = self.radius * cos(self.elevation) * cos(self.azimuth)

        # Update camera position based on the new spherical coordinates
        self.cam = Vec3d([x, y, z, 0])

        # Move the camera to be centered around (0, 0, -8)
        self.cam = self.cam + Vec3d([0, 0, -8, 0])

        # Recalculate view matrix based on the new camera position
        self.calcAgain()

    def dragZ(self, z):
        self.cam = self.cam + self.look_dir.scale(z/100)
        self.radius = sqrt(self.cam.na[0] ** 2 + self.cam.na[1] ** 2 + (self.cam.na[2]+8) ** 2)
        self.calcAgain()


    def calcAgain(self):
        self.look_dir = self.target.__sub__(self.cam).normalize()
        self.horizon = self.look_dir.cross(self.mock_up).normalize()
        self.up = self.horizon.cross(self.look_dir).normalize()

        #calc the matrix
        self.returnViewMatrix()

    def changeLookInto(self,x,y,z):
        #Same logic with rotate around but now we create a coordinate system around the cam not the object.  This is for FPS camera.

        forward = Vec3d([
            cos(self.elevation) * sin(self.azimuth),
            sin(self.elevation),
            cos(self.elevation) * cos(self.azimuth),
            0.0
        ]).normalize()
        right = forward.cross(self.mock_up).normalize()
        up = right.cross(forward).normalize()


        self.target.na[0] = self.target.x() + right.scale(x).x() + up.scale(y).x() + forward.scale(z).x()
        self.target.na[1] = self.target.y() + right.scale(x).y() + up.scale(y).y() + forward.scale(z).y()
        self.target.na[2] = self.target.z() + right.scale(x).z() + up.scale(y).z() + forward.scale(z).z()
        self.cam.na[0] = self.cam.x() + right.scale(x).x() + up.scale(y).x() + forward.scale(z).x()
        self.cam.na[1] = self.cam.y() + right.scale(x).y() + up.scale(y).y() + forward.scale(z).y()
        self.cam.na[2] = self.cam.z() + right.scale(x).z() + up.scale(y).z() + forward.scale(z).z()
        self.radius = np.sqrt(self.cam.x() ** 2 + self.cam.y() ** 2 + (self.cam.z()+8) ** 2)
        self.calcAgain()





