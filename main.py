#CENG-487 Assignment-4
# 320201105 HAMMET POLAT
#05/2025

from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
from object3d import Object3D
from mat3d import Mat3d
from renderer3d import Renderer
from tkinter import filedialog
import sys

instructions = "\n0:To add obj from file\n1:Cube\n2:Plane\n3:Pyramid\n4:Sphere 1\n5:Sphere 2\n6:Tetrahedron\n\n+:Increase the number of grids\n-:Decrease the number of grids\nc: Makes lines vis/unvis\nalt+left mouse: Rotate\nalt+right mouse: Zoom in-out\nf: Reset to original camera matrix\nwasd: Moves camera\nqe: Moves alon z axis\nHit ESC key to quit."
ESCAPE = b'\x1b'


def draw_text_screen(x, y, text, font=GLUT_BITMAP_HELVETICA_18, line_height=20):
    glMatrixMode(GL_PROJECTION)
    glPushMatrix()
    glLoadIdentity()
    gluOrtho2D(0, 640, 0, 480)

    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glLoadIdentity()

    glColor3f(1.0, 1.0, 1.0)

    # Draw each line manually
    lines = text.split('\n')
    for i, line in enumerate(lines):
        glRasterPos2i(x, y - i * line_height)
        for ch in line:
            glutBitmapCharacter(font, ord(ch))

    glMatrixMode(GL_MODELVIEW)
    glPopMatrix()
    glMatrixMode(GL_PROJECTION)
    glPopMatrix()

# Number of the glut window.
window = 0

rot_angle = 0
square_rotation_angle = 0


arr=[]
cube = Object3D(Mat3d())
plane = Object3D(Mat3d())
pyr = Object3D(Mat3d())
sphere = Object3D(Mat3d())
sphere2 = Object3D(Mat3d())
tri = Object3D(Mat3d())
fileObj = Object3D(Mat3d())

arr.extend([cube, plane, pyr,sphere ,sphere2,tri,fileObj])



# A general OpenGL initialization function.  Sets all of the initial parameters.
def InitGL(Width, Height):  # We call this right after our OpenGL window is created.
    glClearColor(0.0, 0.0, 0.0, 0.0)  # This Will Clear The Background Color To Black
    glClearDepth(1.0)  # Enables Clearing Of The Depth Buffer
    glDepthFunc(GL_LESS)  # The Type Of Depth Test To Do
    glEnable(GL_DEPTH_TEST)  # Enables Depth Testing
    glShadeModel(GL_SMOOTH)  # Enables Smooth Color Shading

    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()  # Reset The Projection Matrix
    # Calculate The Aspect Ratio Of The Window
    gluPerspective(45.0, float(Width) / float(Height), 0.1, 100.0)

    glMatrixMode(GL_MODELVIEW)


# The function called when our window is resized (which shouldn't happen if you enable fullscreen, below)
def ReSizeGLScene(Width, Height):
    if Height == 0:  # Prevent A Divide By Zero If The Window Is Too Small
        Height = 1

    glViewport(0, 0, Width, Height)  # Reset The Current Viewport And Perspective Transformation
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(45.0, float(Width) / float(Height), 0.1, 100.0)
    glMatrixMode(GL_MODELVIEW)


last_mouse_x = 0
last_mouse_y = 0
dragging = False
pressed = None

def mouse(button, state, x, y):
    global last_mouse_x, last_mouse_y, dragging,pressed
    modifiers = glutGetModifiers()
    if modifiers & GLUT_ACTIVE_ALT:
        if button == GLUT_LEFT_BUTTON and state == GLUT_DOWN:
            last_mouse_x = x
            last_mouse_y = y
            dragging=True
            pressed = GLUT_LEFT_BUTTON
        elif button == GLUT_RIGHT_BUTTON and state == GLUT_DOWN:
            last_mouse_x = x
            last_mouse_y = y
            dragging = True
            pressed = GLUT_RIGHT_BUTTON
        else:
            dragging = False
    else:
        dragging = False

def motion(x, y):
    global last_mouse_x, last_mouse_y, dragging,pressed

    dx = x - last_mouse_x
    dy = y - last_mouse_y
    last_mouse_x = x
    last_mouse_y = y

    if dragging==True:
        if pressed == GLUT_LEFT_BUTTON:
            Object3D.rotateAround(-dx,dy)
        if pressed == GLUT_RIGHT_BUTTON:
            Object3D.dragZ(dx)


arr[0].unlocked()
current_index = 0

vertex = []
faces = []

def keyPressed(*args):
    global  current_index, vertex, faces
    key = args[0]

    if key == b'\x1b':
        glutLeaveMainLoop()

    if key == b'0':
        filePath = filedialog.askopenfilename()
        if not filePath:
            return  # user cancelled

        vertex = []
        faces = []


        with open(filePath, 'r') as file:
            for line in file:
                if line.startswith('v '):
                    parts = line.strip().split()
                    vertex.append([float(p) for p in parts[1:]])
                elif line.startswith('f '):
                    parts = line.strip().split()
                    face_indices = []
                    for p in parts[1:]:
                        index = p.split('/')[0]
                        face_indices.append(int(index) - 1)
                    faces.append(face_indices)
        file.close()

        fileObj.setVertexes(vertex)
        fileObj.setFaces(faces)
        # Lock the previous object
        arr[current_index].locked()

        # Add the new object to the scene

        current_index = 6
        render.set_selected_object('fileObj')
        arr[current_index].unlocked()

    object_map = {
        b'1': ('cube', 0),
        b'2': ('plane', 1),
        b'3': ('pyramid', 2),
        b'4': ('sphere', 3),
        b'5': ('sphere2', 4),
        b'6': ('triangle', 5),
    }


    if key in object_map:
        new_object, new_index = object_map[key]
        arr[current_index].locked()
        arr[new_index].unlocked()
        render.set_selected_object(new_object)
        current_index = new_index

    elif key == b'+':
        for obj in arr:
            obj.inc_a()

    elif key == b'-':
        for obj in arr:
            obj.dec_a()
    elif key== b'c':
        for obj in arr:
            obj.changeLines()
    elif key == b'f':
        Object3D.resetCam(0, 0, 0.0001)
    elif key == b'w':
        Object3D.changeLookInto(0,0.5,0)
    elif key == b's':
        Object3D.changeLookInto(0,-0.5,0)
    elif key == b'a':
        Object3D.changeLookInto(-0.5,0,0)
    elif key == b'd':
        Object3D.changeLookInto(0.5,0,0)
    elif key == b'q':
        Object3D.changeLookInto(0.0,0,-0.5)
    elif key == b'e':
        Object3D.changeLookInto(0.0,0,0.5)
    else:
        print(f"{key}")





def drawObject( makeMethod):
    firstObject = makeMethod()
    for x in firstObject:

        x.setRotation(0, 0, 0)
        x.setScaling(1, 1, 1)
        x.setTranslation(0.0, 0.0, -8)
        x.applyTransformations()
    arr[0].draw(firstObject)




# The main drawing function.
#global rot_angle, vertex, faces

#rot_angle += 1/4

render = Renderer(arr,draw_text_screen,drawObject)

def DrawGLScene():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)  # Clear the screen


    draw_text_screen(10, 460, instructions)


    render.renderObject()

    glutMouseFunc(mouse)
    glutMotionFunc(motion)

    # Swap buffers to display the drawn object
    glutSwapBuffers()

def main():
    global window
    # For now we just pass glutInit one empty argument. I wasn't sure what should or could be passed in (tuple, list, ...)
    # Once I find out the right stuff based on reading the PyOpenGL source, I'll address this.
    glutInit(sys.argv)

    # Select type of Display mode:
    #  Double buffer
    #  RGBA color
    # Alpha components supported
    # Depth buffer
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)

    # get a 640 x 480 window
    glutInitWindowSize(780, 560)

    # the window starts at the upper left corner of the screen
    glutInitWindowPosition(0, 0)

    # Create the window and store the window id
    window = glutCreateWindow(b"CENG487 Hello Immediate Mode")


    # Register the drawing function with GLUT
    glutDisplayFunc(DrawGLScene)

    # Uncomment this line to get full screen.
    glutFullScreen()

    # When we are doing nothing, redraw the scene.
    glutIdleFunc(DrawGLScene)

    # Register the function called when our window is resized.
    glutReshapeFunc(ReSizeGLScene)

    # Register the function called when the keyboard is pressed.
    glutKeyboardFunc(keyPressed)

    # Initialize OpenGL
    InitGL(640, 480)

    # Start the event processing engine
    glutMainLoop()


# Print message to console, and kick off the main to get it rolling.
main()

