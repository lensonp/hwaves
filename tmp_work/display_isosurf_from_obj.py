import os
import ctypes

import pywavefront
from pywavefront import visualization
import pyglet
from pyglet.gl import *
from hwaves import hwf_plot

win = pyglet.window.Window()
lightfv = ctypes.c_float * 4
obj_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'hwaves_obj','1_0_0.obj')
view_distance = 1.
print(obj_path)


scene = pywavefront.Wavefront(obj_path,strict=True,create_materials=True)
rotation = 0

@win.event
def on_draw():
    win.clear()
    glLoadIdentity()
    glLightfv(GL_LIGHT0, GL_POSITION, lightfv(-1.0, 1.0, 1.0, 0.0))
    glEnable(GL_LIGHT0)
    glTranslated(0.0, 0.0, -1*view_distance)
    glRotatef(rotation, 0.0, 1.0, 0.0)
    glRotatef(-25.0, 1.0, 0.0, 0.0)
    glRotatef(45.0, 0.0, 0.0, 1.0)
    glEnable(GL_LIGHTING)
    visualization.draw(scene)

@win.event
def on_resize(width, height):
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(60., float(width)/height, 1., 100.)
    glMatrixMode(GL_MODELVIEW)
    return True

@win.event
def update(dt):
    global rotation
    rotation += 90.0 * dt
    if rotation > 720.0:
        rotation = 0.0

pyglet.clock.schedule(update)
pyglet.app.run()


