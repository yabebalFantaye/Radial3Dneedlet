# Create a new mayavi scene.
mayavi.new_scene()

# Get the current active scene.
s = mayavi.engine.current_scene

# Read a data file.
d = mayavi.open('fire_ug.vtu')

# Import a few modules.
from mayavi.modules.api import Outline, IsoSurface, Streamline

# Show an outline.
o = Outline()
mayavi.add_module(o)
o.actor.property.color = 1, 0, 0 # red color.

# Make a few contours.
iso = IsoSurface()
mayavi.add_module(iso)
iso.contour.contours = [450, 570]
# Make them translucent.
iso.actor.property.opacity = 0.4
# Show the scalar bar (legend).
iso.module_manager.scalar_lut_manager.show_scalar_bar = True

# A streamline.
st = Streamline()
mayavi.add_module(st)
# Position the seed center.
st.seed.widget.center = 3.5, 0.625, 1.25
st.streamline_type = 'tube'

# Save the resulting image to a PNG file.
s.scene.save('test.png')

# Make an animation:
for i in range(36):
    # Rotate the camera by 10 degrees.
    s.scene.camera.azimuth(10)

    # Resets the camera clipping plane so everything fits and then
    # renders.
    s.scene.reset_zoom()

    # Save the scene.
    s.scene.save_png('anim%d.png'%i)
