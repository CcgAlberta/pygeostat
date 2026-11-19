(view3d)=
# 3D Viewer

PyGeoStat can show drill holes, block models and sub-blocks in 3D. It is built on
[PyVista](https://docs.pyvista.org/), so you do not need to know PyVista to use it. The 3D functions work like
the 2D plots: `gs.view3d()` takes a `DataFile` or `DataFrame` and the name of a variable, and
several calls can share one scene the way 2D plots share `ax`.

:::{seealso}
For every argument, see the {ref}`3D Viewer API Reference <plotting3d>`.
:::

## Installation

PyVista is an optional dependency, so the normal PyGeoStat install does not include it:

```bash
pip install pygeostat[view3d]
```

## Drill holes

`gs.view3d()` returns a `Viewer`. Call `show()` to open it in an interactive window, or `show_static()`
to draw it as a picture on a matplotlib axis. This example uses the oil sands data that ships with
PyGeoStat. The samples of each drill hole are joined with `trace=True`:

```python
import pygeostat as gs

data = gs.ExampleData('oilsands')
viewer = gs.view3d(data, 'Bitumen', trace=True, radius=15)
viewer.add_ruler()
viewer.set_vertical_exaggeration(10)
viewer.show()                    # interactive window
fig, ax = viewer.show_static()   # or a matplotlib figure
```

![Drill holes coloured by bitumen](../figures/view3d/drillholes.png)

How the drill holes are drawn depends on the columns of your data:

- **Intervals** with a start and an end: `gs.view3d(data, 'Au', start=('xs', 'ys', 'zs'), end=('xe', 'ye', 'ze'))`
  draws one cylinder per interval.
- **Samples** along a hole: `trace=True` joins the samples of each drill hole (from `DataFile.dh`) with
  cylinders, each with the value of its upper sample. Samples further apart than `max_gap` (by default
  three times the typical spacing) are not joined.
- **Points**: with neither, each row is a point at the coordinates of the `DataFile`.

The coordinates, drill hole id and variable are read from the `DataFile` when it knows them. For a
`DataFrame`, pass `x`, `y`, `z` and `dh`.

## Block models

A block model on a regular grid needs its `GridDef` (it is read from `DataFile.griddef` when set).
Sub-blocks of any size need one row per block, with the block centre and its three sizes:

```python
blocks = gs.ExampleData('3d_grid')
griddef = gs.GridDef([40, 1, 2, 40, 1, 2, 40, 0.5, 1])
gs.view3d(blocks, 'value', griddef=griddef)

# sub-blocks: columns x, y, z (centre) and xsize, ysize, zsize
gs.view3d(subblocks, 'Au', x='x', y='y', z='z', size=('xsize', 'ysize', 'zsize'))
```

Every row is drawn where its coordinates say. A block with a wrong coordinate shows up away from the
others; it is not an error.

## Layers and slices

Pass `viewer=` to add to a scene. `orient`, `slice_number` and `slice_thickness` work as in
{func}`gs.slice_plot() <pygeostat.plotting.slice_plot.slice_plot>`: `orient` is `'xy'`, `'xz'` or `'yz'`.
This puts a slice of the model inside a see-through copy of it, next to the 2D plot of the same slice:

```python
import matplotlib.pyplot as plt

viewer = gs.view3d(blocks, 'value', griddef=griddef, opacity=0.15)
gs.view3d(blocks, 'value', griddef=griddef, orient='xz', slice_number=20, viewer=viewer)
viewer.add_ruler()

fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
gs.slice_plot(blocks, var='value', griddef=griddef, orient='xz', slice_number=20, ax=left)
viewer.show_static(ax=right, title='3D')
```

![A slice of a block model, in 2D and in 3D](../figures/view3d/layers_and_slice.png)

- **Data without a grid** (drill holes, sub-blocks) are sliced with `slice_coordinate=` (the position along
  the axis not in `orient`) and `pointtol=` (how far either side of it to show). Drill holes are cut at
  that distance.
- **`threshold=(min, max)`** shows only values in that range.
- A slice or a threshold keeps the colours of the full range of the variable (unless you pass `clim`), and
  the ruler keeps the whole model, so moving the slice does not change what the colours mean.

## Categories

Pass `catdict` to colour a categorical variable (a rock type, for example). The colours match
{func}`gs.location_plot() <pygeostat.plotting.location_plot.location_plot>`. If the `DataFile` has a `cat`
column, `gs.view3d(data, data.cat)` uses its `catdict` without anything else:

```python
data = gs.ExampleData('point3d_ind_mv')
catdict = {1: 'Sand', 2: 'Breccia', 3: 'Silt', 4: 'Shale', 5: 'Coal'}
viewer = gs.view3d(data, 'Lithofacies', trace=True, radius=4, catdict=catdict)
viewer.add_ruler()
viewer.set_vertical_exaggeration(10)
```

![Drill holes coloured by rock type](../figures/view3d/categories.png)

## Axes, exaggeration and views

- `viewer.add_ruler()` draws axes titled from `gs.Parameters` (`plotting.xname`, `plotting.unit`, ...), with
  labels on round numbers.
- `viewer.set_vertical_exaggeration(10)` stretches the scene vertically. The ruler keeps the real elevations.
- `viewer.get_camera()` returns the current view as a dictionary of plain numbers. Save it with `json`,
  and give it back to `viewer.set_camera()` to get the same view again, in `show()` and in `show_static()`.

## Static figures

`viewer.show_static()` returns a matplotlib `(fig, ax)`, so a 3D view can sit in a figure with 2D plots,
have a legend and be saved like any matplotlib figure. It draws the scene off screen, with the ruler, the
vertical exaggeration and the camera of the viewer. It works without a display, which is also how the
tests run.

## Interactive panel

`viewer.panel()` shows the scene with a panel of controls beside it. It needs a few more packages
(trame), which are optional too:

```bash
pip install "pygeostat[view3d-panel]"
```

```python
import pygeostat as gs

blocks = gs.ExampleData('3d_grid')
griddef = gs.GridDef([40, 1, 2, 40, 1, 2, 40, 0.5, 1])
viewer = gs.view3d(blocks, 'value', griddef=griddef)
gs.view3d(gs.ExampleData('3d_estimate'), 'Estimate', viewer=viewer)
viewer.add_ruler()
viewer.panel()    # last line of a notebook cell
```

For each layer the panel has:

- a checkbox to show or hide it, and an opacity slider;
- the colour limits (for a layer coloured by a continuous variable), with the range of the data in
  each box label;
- a **Threshold** checkbox and a typed minimum: with the box ticked, only the values from the minimum
  up are shown. The box starts at the lowest value of the data.
- a **Colormap** list with a small preview of each colormap (for a layer coloured by categories, a
  **Palette** list of the palettes that have enough colours);
- **Colour by**, when the layer has other variables: it is coloured by the one you choose, with its own
  colour bar. The threshold keeps using the variable the layer was made with. By default a layer made
  from a `DataFile` keeps its numeric variables and its category column; use
  `gs.view3d(..., other_vars=[...])` to choose them;
- for a layer coloured by categories, a button for each category in its colour: click one to show only
  that category, click more to add them, and **Show all** shows every category again;
- **Show cell edges**, for a block layer;
- a **Tube radius** slider, for a layer of drill holes. The tubes are drawn again when you let go.

Values that depend on the data are typed, and applied on Enter or when you leave the box. A value that
does not work, such as a minimum above the maximum or a threshold that leaves nothing, is explained in
the panel and changes nothing.

**Slicer.** One plane for the whole scene. Tick **Enable slicer** and move the plane with the Origin X,
Y and Z sliders, tilt it with the Normal X, Y and Z sliders, or press **XY**, **XZ** or **YZ** to face
it along an axis. Block layers then show only the cut through the plane, and the drill holes are
hidden. **Show drill holes on slice** shows the holes that cross the plane: only the holes whose tube
reaches it, and of those the part within the **Tolerance** (the thickness, in the units of your data,
of the band around the plane). The cut follows the sliders as you drag them.

With a ruler (`viewer.add_ruler()`) the panel also has boxes for East, North and Elev, and **All axes
from 0**. The ruler then counts from that point instead of from the coordinates, which helps to
compare two models whose coordinates are far apart. The axis titles say where the labels start.

The name of a layer is its `label` (or its variable), so the panel shows the names of your own data,
for example `gs.view3d(blocks, 'value', griddef=griddef, label='Block model value')`.

**Scene.** The **Vertical exaggeration** slider stretches the scene, and **Background** is white or
black (the text of the ruler and the colour bars follows). **Outline** draws a box around the scene,
**Axis** the ruler, **X/Y/Z arrows in the scene** large arrows at the lowest corner and **Coordinate
legend (corner)** the small marker that shows the directions. **Axis ticks** sets the spacing of the
ruler's ticks along East/North and along Elevation: the same spacing in two pictures gives grid lines the
same distance apart. A spacing that does not fit is replaced by a multiple of it, and the panel says so.

**Camera.** **Top 2D**, **East-Elev 2D** and **North-Elev 2D** are flat views (a parallel projection)
that leave out the axis they look along. **Iso 3D**, **X 3D**, **Y 3D** and **Z 3D** are perspective
views. **Parallel projection** switches the projection of the current view, and **Fit view** frames
the whole scene again.

**Picking.** Click a drill hole (a click, not a drag) to see its values under **Picking**: the hole, the
elevation of the interval, the value of each variable and the name of its category. The interval is
marked with a wider tube. Blocks cannot be picked, so a click reaches a hole inside or behind a block
model.

**Export.** **Screenshot (PNG)** downloads the picture, and **Interactive HTML** downloads one web page
that opens in any browser and can be turned and zoomed. The page has the colour bars and the ruler with
the real coordinates, but not the slicer, the panel or picking.

**Saved views.** **Save** writes the whole view (the camera, the layers' settings and the ruler) to a
json file in the folder `views/`, or in `viewer.panel(views_folder='...')`, and **Load** shows a saved
view again. Without a name, a name is made from the time.

### The same from a script

Each control is a method of the viewer, so nothing needs to be clicked to get a view:

```python
viewer.set_layer('value', opacity=0.3, clim=(-2, 2))
viewer.set_layer('value', threshold=(0.5, 1.0))                   # False shows all values again
viewer.set_layer('Estimate', visible=False)
viewer.set_slicer(enabled=True, origin=(40, 40, 10), normal=(1, 0, 1))   # the plane of the slicer
viewer.set_slicer(tolerance=20, show_holes=True)
viewer.set_slicer(enabled=False)
viewer.set_axis_origin((0, 0, 0), from_zero=True)
viewer.set_axis_ticks(xy=20, z=5)                                 # no arguments: round spacings that suit the scene
viewer.set_vertical_exaggeration(3)
viewer.set_background('black')
viewer.set_scene(outline=True, axis_arrows=True, corner_axes=True, ruler=True)
viewer.set_camera_view('top')                                     # 'iso', 'x', 'y', 'z', 'top', 'east' or 'north'
viewer.set_parallel_projection(True)
viewer.fit_view()

viewer.set_layer('value', cmap='plasma', edges=True)              # colormap, cell edges of a block layer
# for a layer of drill holes called 'holes', with the variables 'Cu' and 'Rock' (categories)
viewer.set_layer('holes', radius=4, color_by='Cu')                # tube radius, another variable
viewer.set_layer('holes', color_by='Rock', cmap='cat_pastel')     # a variable with categories and its palette
viewer.set_layer('holes', categories=['Granite'])                 # show only these; [] shows all again

viewer.save_screenshot('scene.png')
viewer.export_html('scene.html')

import json
view = viewer.get_view()                  # plain numbers and strings, so json works
json.dump(view, open('view.json', 'w'))
viewer.set_view(json.load(open('view.json')))
```

`set_view` matches layers by name and refuses a view that names a layer the scene does not have.

A single layer can also be cut without the slicer, with the arguments of `gs.view3d()`:
`viewer.set_layer('value', orient='xz', slice_number=20, slice_thickness=2)`, and
`orient=False` shows the whole layer again.

## What is not included yet

Surfaces from DXF files, picking a block to see its values, and display options (drawing styles,
choosing categories from the legend, the background colour, a flat 2D view of an axis) will follow
in separate issues.
