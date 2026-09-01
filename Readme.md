# FTPlot

`FTPlot` is a Matplotlib-based utility for creating aircraft flight-test and simulation time-history plots with multiple independent Y axes.

It provides:

- A common X/time axis.
- Multiple independently scaled Y axes.
- Left- and right-side Y axes.
- Automatic Y-axis scaling.
- Custom axis vertical positions and heights.
- Optional curve labels.
- An interactive vertical cursor.
- An optional slider for moving the cursor.
- Value boxes showing interpolated curve values at the cursor position.

## Requirements

```text
numpy
matplotlib
```

Install with:

```bash
pip install numpy matplotlib
```

## Basic Usage

```python
import numpy as np
import matplotlib.pyplot as plt

from ftplot import FTPlot

fig, ax = plt.subplots()

plot = FTPlot(
    fig,
    ax,
    XLim=(0, 10),
    slider=True,
    Xunit="s",
)

plot.AddAxis(
    Name="Airspeed",
    GridHeight=5,
    GridPos=0,
    Unit="m/s",
)

time = np.linspace(0, 10, 500)
airspeed = 20 + 3 * np.sin(time)

plot.AddCurve(
    Name="V",
    Axis="Airspeed",
    Xdata=time,
    Ydata=airspeed,
)

plt.show()
```

---

# Class Reference

## `FTPlot`

```python
FTPlot(
    fig,
    ax,
    Ngridx=10,
    Ngridy=13,
    XLim=(0, 1),
    slider=False,
    Xunit="s",
)
```

Creates and configures an `FTPlot` object.

### Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fig` | `matplotlib.figure.Figure` | — | Matplotlib figure containing the plot. |
| `ax` | `matplotlib.axes.Axes` | — | Main Matplotlib axis used for the common X axis and grid. |
| `Ngridx` | `int` | `10` | Number of major intervals on the X axis. |
| `Ngridy` | `int` | `13` | Number of vertical grid units used to position additional axes. |
| `XLim` | tuple | `(0, 1)` | Limits of the common X axis. |
| `slider` | `bool` | `False` | Enables the interactive vertical cursor and slider when `True`. |
| `Xunit` | `str` | `"s"` | Unit displayed on the X axis and slider. |

### Main Attributes

| Attribute | Description |
|---|---|
| `fig` | Reference to the Matplotlib figure. |
| `ax` | Reference to the main axis. |
| `XLim` | Common X-axis limits. |
| `Ngridx` | Number of X grid intervals. |
| `Ngridy` | Number of Y grid intervals. |
| `Xunit` | X-axis unit. |
| `extent` | Position of the original axis within the figure. |
| `Axis` | Dictionary containing all additional Y axes. |
| `Curve` | Dictionary containing all curves. |
| `vline` | Interactive vertical cursor, when enabled. |
| `slider_ax` | Axis containing the slider, when enabled. |
| `x_slider` | Matplotlib `Slider` object, when enabled. |

The constructor also configures major/minor grids and removes Y tick labels from the main axis.

---

# Adding Y Axes

## `AddAxis`

```python
AddAxis(
    Name,
    GridHeight=4,
    GridPos=0,
    Unit="",
    Position="Left",
    YLims="Auto",
    offset=0.02,
)
```

Adds an independent Y axis to the figure.

The additional axes are implemented as separate Matplotlib `Axes` objects overlaid on the original plotting area.

### Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `Name` | `str` | — | Unique identifier for the axis. |
| `GridHeight` | `int` | `4` | Height of the axis in units of the internal Y grid. |
| `GridPos` | `float` | `0` | Bottom position of the axis in units of the internal Y grid. |
| `Unit` | `str` | `""` | Unit displayed below the axis name. |
| `Position` | `str` | `"Left"` | `"Left"` places the axis on the left; any other value places it on the right. |
| `YLims` | tuple/list or `"Auto"` | `"Auto"` | Initial Y limits. `"Auto"` enables automatic scaling. |
| `offset` | `float` | `0.02` | Offset of the Y-axis spine outside the main plotting area. |

### Grid-based positioning

`Ngridy` defines a virtual vertical grid.

For example:

```python
plot = FTPlot(fig, ax, Ngridy=13)

plot.AddAxis(
    Name="Velocity",
    GridHeight=5,
    GridPos=0,
    Unit="m/s",
)
```

places the axis over approximately the bottom `5/13` of the plotting area.

A second axis can be placed above it:

```python
plot.AddAxis(
    Name="Altitude",
    GridHeight=5,
    GridPos=5,
    Unit="m",
)
```

This makes it possible to construct vertically separated flight-test plot bands while retaining a common X axis.

### Automatic scaling

If:

```python
YLims="Auto"
```

the axis is automatically scaled based on the curves added to it.

The implementation uses Matplotlib's:

```python
relim()
autoscale_view()
```

followed by `MaxNLocator` to control the number of major tick intervals.

If explicit limits are provided:

```python
YLims=(-10, 30)
```

automatic scaling is disabled.

### Duplicate axis names

Axis names must be unique. If `Name` already exists, the method prints:

```text
Name already exists
```

and does not create a new axis.

---

# Adding Curves

## `AddCurve`

```python
AddCurve(
    Name,
    Axis,
    Xdata,
    Ydata,
    Label=None,
    **kwargs,
)
```

Adds a curve to an existing Y axis.

### Parameters

| Parameter | Type | Description |
|---|---|---|
| `Name` | `str` | Unique identifier for the curve. |
| `Axis` | `str` | Name of the axis on which the curve is plotted. |
| `Xdata` | array-like | X coordinates. |
| `Ydata` | array-like | Y coordinates. |
| `Label` | optional | If not `None`, a text label is created on the curve. |
| `**kwargs` | dict | Additional arguments passed to Matplotlib's `plot()`. |

### Example

```python
plot.AddCurve(
    Name="V",
    Axis="Airspeed",
    Xdata=time,
    Ydata=airspeed,
    linewidth=2,
)
```

Matplotlib styling options can be passed directly:

```python
plot.AddCurve(
    Name="V",
    Axis="Airspeed",
    Xdata=time,
    Ydata=airspeed,
    color="r",
    linestyle="--",
    linewidth=2,
)
```

### Curve labels

When `Label` is not `None`, a text object is placed at a randomly selected location between approximately 25% and 75% of the data.

**Current implementation detail:** the supplied `Label` value is only used as a trigger. The displayed text is actually `Name`.

For example:

```python
plot.AddCurve(
    Name="Nz",
    Axis="Load Factor",
    Xdata=time,
    Ydata=nz,
    Label="$N_z$",
)
```

currently displays `Nz`, not `$N_z$`.

If the intention is to display the supplied `Label`, the implementation should replace:

```python
Name
```

with:

```python
Label
```

when creating the text object.

### Duplicate curve names

Curve names must be unique. A duplicate produces:

```text
Name already exists
```

and the second curve is not added.

---

# Interactive Cursor

The interactive functionality is enabled with:

```python
plot = FTPlot(fig, ax, slider=True)
```

When enabled, the class creates:

1. A vertical dashed cursor.
2. A horizontal slider below the plot.
3. A callback connecting the slider to `update()`.

The initial cursor position is:

```python
XLim[0]
```

The slider displays the current X value using:

```text
%.2f <Xunit>
```

For example:

```text
5.27 s
```

---

# `update`

```python
update(val)
```

Callback used by the slider.

It:

1. Moves the vertical cursor to `val`.
2. Calls `updateDataBoxes(val)`.
3. Requests a redraw with `draw_idle()`.

The method is normally called automatically by Matplotlib and does not need to be called manually.

---

# `updateDataBoxes`

```python
updateDataBoxes(x_value)
```

Updates the value box associated with every curve.

For each curve, the method obtains its plotted data and computes:

```python
y_value = np.interp(x_value, x, y)
```

The value box is moved to:

```text
(x_value, y_value)
```

and displays the value with one decimal place:

```text
12.3
```

If the cursor is outside the current X limits, the value box is hidden.

### Important requirement

Because `np.interp()` is used, curve X data should normally be monotonically increasing:

```python
time[0] < time[1] < ... < time[-1]
```

This is naturally satisfied by most flight-test time histories.

---

# `UpdateAxis`

```python
UpdateAxis(Axis)
```

Recalculates the Y-axis scaling for the specified axis.

This method is automatically called by `AddCurve()`.

For automatically scaled axes it:

1. Recalculates the data limits.
2. Determines suitable Y-axis tick positions.
3. Updates the Y limits.
4. Creates three main Y ticks.
5. Adjusts the visible axis spine to match the central portion of the scale.

The number of requested major tick intervals is:

```python
GridHeight + 1
```

using Matplotlib's `MaxNLocator`.

---

# Internal Data Structures

## `Axis`

`Axis` is a dictionary indexed by axis name.

For example:

```python
plot.Axis["Airspeed"]
```

contains:

```python
{
    "ax": <matplotlib.axes.Axes>,
    "Name": "Airspeed",
    "Unit": "m/s",
    "AutoScale": True,
    "Position": "Left",
    "GridHeight": 5,
    "GridPos": 0,
}
```

### Keys

| Key | Description |
|---|---|
| `"ax"` | Matplotlib `Axes` object. |
| `"Name"` | Axis identifier/name. |
| `"Unit"` | Axis unit. |
| `"AutoScale"` | Whether automatic Y scaling is enabled. |
| `"Position"` | Left/right placement. |
| `"GridHeight"` | Height in virtual grid units. |
| `"GridPos"` | Vertical position in virtual grid units. |

---

## `Curve`

`Curve` is a dictionary indexed by curve name.

For example:

```python
plot.Curve["V"]
```

contains:

```python
{
    "Curve": <matplotlib.lines.Line2D>,
    "Axis": "Airspeed",
    "Label": <matplotlib.text.Text>,
    "ValueBox": <matplotlib.text.Text>,
}
```

### Keys

| Key | Description |
|---|---|
| `"Curve"` | Matplotlib `Line2D` object. |
| `"Axis"` | Axis containing the curve. |
| `"Label"` | Curve label text object, or `None`. |
| `"ValueBox"` | Text object used for the interactive value display. |

---

# Complete Example

```python
import numpy as np
import matplotlib.pyplot as plt

from ftplot import FTPlot

plt.close("all")

fig, ax = plt.subplots()

plot = FTPlot(
    fig,
    ax,
    Ngridx=10,
    Ngridy=13,
    XLim=(0, 1),
    slider=True,
    Xunit="s",
)

plot.AddAxis(
    Name="Velocity",
    GridHeight=6,
    GridPos=0.5,
    Unit="m/s",
)

plot.AddAxis(
    Name="Altitude",
    GridHeight=4,
    GridPos=3,
    Unit="m",
    Position="Right",
)

plot.AddAxis(
    Name="Angle",
    GridHeight=4,
    GridPos=5,
    Unit="deg",
    offset=0.1,
)

Xdata = np.linspace(0, 1, 500)

velocity = 3.7 * np.cos(13 * Xdata)
altitude = 100 + 20 * np.sin(5 * Xdata)
angle = 5 * np.sin(10 * Xdata)

plot.AddCurve(
    Name="V",
    Axis="Velocity",
    Xdata=Xdata,
    Ydata=velocity,
    Label="V",
)

plot.AddCurve(
    Name="H",
    Axis="Altitude",
    Xdata=Xdata,
    Ydata=altitude,
)

plot.AddCurve(
    Name="Theta",
    Axis="Angle",
    Xdata=Xdata,
    Ydata=angle,
)

plt.show()
```

---

# Typical Flight-Test Application

A typical flight-test plot can use one axis for each group of quantities with similar units and scale:

```python
fig, ax = plt.subplots()

plot = FTPlot(
    fig,
    ax,
    XLim=(time[0], time[-1]),
    Ngridx=10,
    Ngridy=13,
    slider=True,
    Xunit="s",
)

plot.AddAxis(
    Name="Airspeed",
    GridHeight=4,
    GridPos=0,
    Unit="m/s",
)

plot.AddAxis(
    Name="Altitude",
    GridHeight=4,
    GridPos=4,
    Unit="m",
    Position="Right",
)

plot.AddAxis(
    Name="Pitch",
    GridHeight=4,
    GridPos=8,
    Unit="deg",
)

plot.AddCurve("V", "Airspeed", time, airspeed)
plot.AddCurve("H", "Altitude", time, altitude)
plot.AddCurve("Theta", "Pitch", time, pitch_angle)

plt.show()
```

This architecture is useful for displaying aircraft parameters such as:

- Airspeed
- Altitude
- Angle of attack
- Pitch attitude
- Roll attitude
- Load factor
- Angular rates
- Control-surface positions
- Engine parameters

against a common time reference.

---

# Method Summary

| Method | Purpose |
|---|---|
| `FTPlot(...)` | Create and configure the plotting object. |
| `AddAxis(...)` | Add an independent Y axis. |
| `AddCurve(...)` | Plot a curve on an existing Y axis. |
| `UpdateAxis(...)` | Update automatic scaling of a Y axis. |
| `updateDataBoxes(...)` | Update curve values at a selected X position. |
| `update(...)` | Slider callback that moves the cursor and updates values. |

---

# Dependencies

The implementation uses:

```python
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.ticker import AutoMinorLocator, MaxNLocator
```

### NumPy

Used for:

- Array operations.
- `linspace()`.
- `interp()`.
- Mean/range calculations.
- Tick generation.

### Matplotlib

Used for:

- Figure and axis management.
- Curve plotting.
- Axis creation.
- Text labels.
- Grid lines.
- Axis scaling.
- Interactive slider functionality.

---

# Design Notes

The class does **not** use Matplotlib's conventional `twinx()` mechanism. Instead, each Y axis is a separate `Axes` object positioned explicitly over the original plotting region.

This provides greater control over:

- Axis vertical position.
- Axis height.
- Left/right placement.
- Axis spine offsets.
- Independent scaling.

The resulting structure is well suited to flight-test plots where many parameters have different units and numerical ranges but share a common time axis.

## Current Limitations

The current implementation has a few assumptions worth keeping in mind:

1. `AddCurve()` assumes the specified axis already exists.
2. `np.interp()` assumes appropriately ordered X data.
3. The curve-label position is selected randomly, so the figure can change slightly between executions.
4. `Label` currently acts as an enable/disable flag rather than the actual displayed label text.
5. The right-side axis is selected by any `Position` value other than `"Left"`; there is no explicit validation.
6. The class assumes that the input X/Y arrays are compatible with Matplotlib's plotting interface.
7. Duplicate axis and curve names are reported with `print()` rather than raising an exception.

---

# License

No license information is included in the supplied source code.

If `FTPlot` is distributed as part of a software package, a project license should be added separately.
