#CENG-487 Assignment-3
# 320201105 HAMMET POLAT
#05/2025

from OpenGL.GLUT import *
from object3d import Object3D
from geometry import Geo

class Renderer:
    def __init__(self, objects, draw_text_func, draw_object_func):
        self.objects = objects
        self.selected_object = 'cube'
        self.draw_text_screen = draw_text_func
        self.drawObject = draw_object_func

    def set_selected_object(self, selected_object):
        self.selected_object = selected_object

    def renderObject(self):
        object_map = {
            'cube': (self.objects[0], self.objects[0].makeCube),
            'plane': (self.objects[1], self.objects[1].makePlane),
            'pyramid': (self.objects[2], self.objects[2].makePyramid),
            'sphere': (self.objects[3], self.objects[3].makeSphere),
            'sphere2': (self.objects[4], self.objects[4].makeSphere2),
            'triangle': (self.objects[5], self.objects[5].makeTriangle),
            'fileObj': (self.objects[6], self.objects[6].makeParse),
        }

        if self.selected_object in object_map:
            obj, make_method = object_map[self.selected_object]
            self.draw_text_screen(10, 460, f"Level of subdiv: {obj.get_a()}")
            self.drawObject(make_method)