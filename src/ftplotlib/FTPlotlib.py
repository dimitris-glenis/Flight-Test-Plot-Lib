# -*- coding: utf-8 -*-
"""
Created on Tue May 20 11:51:17 2025

@author: aantonak
"""
import matplotlib
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.ticker import AutoMinorLocator, MaxNLocator



class FTPlot:
    """
    Aircraft Simulator for Spin
    """

    def __init__(self, fig, ax,
                 Ngridx=10, Ngridy=13,
                 XLim=(0, 1), slider=False, Xunit="s"):
        # Core references
        self.fig = fig
        self.ax = ax
        self.XLim = XLim
        self.Ngridx = Ngridx
        self.Ngridy = Ngridy
        self.Xunit = Xunit
        self.slider = slider

        # Reserve space for grids and optional slider
        plt.subplots_adjust(left=0.15, right=0.85, top=0.95)
        if slider:
            ext = ax.get_position()
            plt.subplots_adjust(
                bottom=0.15,
                top=ext.y1 + 0.15 - ext.y0
            )

        # Main axes grid
        ax.set_xlim(self.XLim)
        ax.set_ylim(0, 1)
        ax.set_xticks(np.linspace(self.XLim[0], self.XLim[1], self.Ngridx + 1))
        ax.set_yticks(np.linspace(0, 1, self.Ngridy + 1))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.set_yticklabels([])
        ax.set_xlabel(f"Time ({self.Xunit})")
        ax.grid(True, which='major', ls='-', color='k', alpha=0.3)
        ax.grid(True, which='minor', ls=':', color='k', alpha=0.3)

        # Slider line
        if slider:
            self.initial_x = self.XLim[0]
            self.vline = ax.axvline(self.initial_x, color='red', linestyle='--')
            ext = ax.get_position()
            self.slider_ax = fig.add_axes([ext.x0, 0.02, ext.width, 0.03])
            self.x_slider = Slider(
                self.slider_ax, "",
                self.XLim[0], self.XLim[1],
                valinit=self.initial_x,
                valfmt=f"%.2f {self.Xunit}"
            )
            self.x_slider.on_changed(self.update)

        # Containers
        self.Axis = {}
        self.Curve = {}
        self.fig.canvas.mpl_connect('resize_event', self._on_resize)

    def updateDataBoxes(self, t):
        for name, data in self.Curve.items():
            # Skip vertical lines - their value boxes should remain fixed
            if data.get('IsVerticalLine', False):
                continue

            # For regular curves, update as before
            x, y = data['Curve'].get_xdata(), data['Curve'].get_ydata()
            yv = np.interp(t, x, y)
            vb = data['ValueBox']
            vb.set_position((t, yv))
            vb.set_text(f"{yv:.3f}" if self.XLim[0] <= t <= self.XLim[1] else "")

    def update(self, val):
        self.vline.set_xdata([val, val])
        self.updateDataBoxes(val)
        self.fig.canvas.draw_idle()

    def AddAxis(self, Name,
                GridHeight=1, GridPos=0,
                Unit="", Position="Left",
                YLims="Auto", offset=0.02):
        """
        Add a sub-axis spanning GridHeight rectangles,
        centered on the GridPos-th rectangle.
        """
        if Name in self.Axis:
            print(f"Axis '{Name}' already exists")
            return

        self.fig.canvas.draw()
        ext = self.ax.get_position()
        cell_h = ext.height / self.Ngridy
        center = ext.y0 + (GridPos + 0.5) * cell_h
        height = GridHeight * cell_h
        bottom = center - height / 2
        ax2 = self.fig.add_axes([ext.x0, bottom, ext.width, height])

        # Initial limits
        lims = [-1, 1] if YLims == 'Auto' else YLims
        ax2.set_ylim(lims)
        ax2.set_autoscaley_on(True)
        lo, hi = lims
        ax2.set_yticks([lo, (lo + hi) / 2, hi])

        ax2.set_xlim(self.XLim)
        ax2.set_xticks([]); ax2.set_xticklabels([])
        ax2.set_ylabel(f"{Name}\n{Unit}", labelpad=0)
        for s in ('top', 'bottom'):
            ax2.spines[s].set_visible(False)
        ax2.patch.set_visible(False)
        if Position.lower() == 'left':
            ax2.spines['right'].set_visible(False)
            ax2.spines['left'].set_position(('axes', -offset))
        else:
            ax2.yaxis.tick_right()
            ax2.yaxis.set_label_position('right')
            ax2.spines['left'].set_visible(False)
            ax2.spines['right'].set_position(('axes', 1 + offset))

        self.Axis[Name] = {
            'ax': ax2,
            'AutoScale': True,
            'GridHeight': GridHeight,
            'GridPos': GridPos,
            'Position': Position
        }

    def AddCurve(self, Name, Axis, Xdata, Ydata, **kwargs):
        """
        Plot a curve on a sub-axis,
        then autoscale to integer bounds plus one-unit padding both sides.
        """
        ax2 = self.Axis[Axis]['ax']
        (line,) = ax2.plot(Xdata, Ydata, **kwargs)
        # Label
        idx = len(Xdata) // 2
        ax2.text(
            Xdata[idx], Ydata[idx], Name,
            ha='center', va='center', color=line.get_color(),
            bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1')
        )
        vb = ax2.text(
            Xdata[0], Ydata[0], '',
            ha='center', va='center', color=line.get_color(),
            bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1', alpha=1.0)
        )
        self.Curve[Name] = {'Curve': line, 'ValueBox': vb}

        # Auto-rescale with padding
        info = self.Axis[Axis]
        if info['AutoScale']:
            dmin = np.min(Ydata)
            dmax = np.max(Ydata)
            # Integer floor/ceil
            flo = np.floor(dmin)
            cei = np.ceil(dmax)
            # Pad one unit both sides
            lo = flo - 1
            hi = cei + 1
            ax2.set_ylim(lo, hi)
            ax2.set_yticks([lo, (lo + hi) / 2, hi])

    def AddVerticalLine(self, Name, Xpos, YLims=None, linestyle='--', marker=None, **kwargs):
        """
        Add a vertical line at Xpos on the main axis.
        Optionally specify YLims to set the vertical limits of the line.
        """
        if Name in self.Curve:
            print(f"Vertical line '{Name}' already exists")
            return

        if YLims is None:
            YLims = self.ax.get_ylim()

        # Apply default linestyle if not specified in kwargs
        if 'linestyle' not in kwargs and 'ls' not in kwargs:
            kwargs['linestyle'] = linestyle

        line = self.ax.axvline(Xpos, ymin=YLims[0], ymax=YLims[1], **kwargs)

        y_pos = YLims[0] - 0.05
        vb = self.ax.text(
            Xpos, y_pos, f"{Xpos:.3f}",
            ha='center', va='top', color=line.get_color(),
            bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1', alpha=1.0)
        )

        # Initialize intersections dictionary
        intersections = {}

        # Create dictionary to store intersection data
        self.Curve[Name] = {
            'Curve': line,
            'ValueBox': vb,
            'IsVerticalLine': True,
            'XPosition': Xpos,
            'Intersections': intersections
        }

        # Add markers at intersections with existing curves
        for curve_name, data in self.Curve.items():
            # Skip the vertical line itself and other vertical lines
            if curve_name == Name or data.get('IsVerticalLine', False):
                continue

            try:
                # Get curve data
                curve = data['Curve']
                x_data = curve.get_xdata()
                y_data = curve.get_ydata()

                # Check if Xpos is within curve x range
                if Xpos >= min(x_data) and Xpos <= max(x_data):
                    # Interpolate to find y value at the vertical line position
                    y_val = np.interp(Xpos, x_data, y_data)

                    # Add marker at intersection point if marker is specified
                    mark = None
                    if marker is not None:
                        mark = curve.axes.plot([Xpos], [y_val], marker, color=curve.get_color(),
                                               markersize=6, zorder=10)[0]

                    # Create a text annotation for the intersection
                    intersection_box = curve.axes.text(
                        Xpos, y_val, f"{y_val:.3f}",
                        ha='left', va='bottom', color=curve.get_color(),
                        bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1', alpha=1.0)
                    )

                    # Store references to the marker and text
                    intersections[curve_name] = {
                        'marker': mark,
                        'text': intersection_box
                    }
            except:
                pass

        self.fig.canvas.draw_idle()

    def RemoveCurve(self, Name):
        if Name not in self.Curve:
            print(f"No curve named '{Name}' exists")
            return

        entry = self.Curve.pop(Name)

        # Remove the curve itself
        if 'Curve' in entry and entry['Curve'] is not None:
            entry['Curve'].remove()

        # Remove the value box text
        if 'ValueBox' in entry and entry['ValueBox'] is not None:
            entry['ValueBox'].remove()

        # Remove intersection markers and value boxes if this is a vertical line
        if 'Intersections' in entry:
            for intersection in entry['Intersections'].values():
                if 'marker' in intersection and intersection['marker'] is not None:
                    intersection['marker'].remove()
                if 'text' in intersection and intersection['text'] is not None:
                    intersection['text'].remove()

        # Remove intersection references if this is a regular curve
        for other_name, other_data in self.Curve.items():
            if other_data.get('IsVerticalLine', False) and 'Intersections' in other_data:
                if Name in other_data['Intersections']:
                    other_data['Intersections'][Name]['marker'].remove()
                    other_data['Intersections'][Name]['text'].remove()
                    del other_data['Intersections'][Name]

        # Find and remove the curve label
        for ax in self.fig.axes:
            for text in ax.texts:
                if text.get_text() == Name:
                    text.remove()
                    break

        self.fig.canvas.draw_idle()

    def enable_autoscale(self, Name):
        """Re-enable autoscale for sub-axis and apply padding."""
        if Name not in self.Axis:
            print(f"No such axis '{Name}'")
            return
        info = self.Axis[Name]
        info['AutoScale'] = True
        ax2 = info['ax']
        lines = ax2.get_lines()
        if not lines:
            return
        ys = np.hstack([ln.get_ydata() for ln in lines])
        flo = np.floor(ys.min())
        cei = np.ceil(ys.max())
        lo = flo - 1
        hi = cei + 1
        ax2.set_ylim(lo, hi)
        ax2.set_yticks([lo, (lo + hi) / 2, hi])

    def _on_resize(self, event):
        for info in self.Axis.values():
            if not info['AutoScale']:
                continue
            ax2 = info['ax']
            lines = ax2.get_lines()
            if not lines:
                continue
            ys = np.hstack([ln.get_ydata() for ln in lines])
            flo = np.floor(ys.min())
            cei = np.ceil(ys.max())
            lo = flo - 1
            hi = cei + 1
            ax2.set_ylim(lo, hi)
            ax2.set_yticks([lo, (lo + hi) / 2, hi])




if __name__ == "__main__":

    matplotlib.use('Qt5Agg')
    plt.ion()
    plt.close('all')
    fig, ax = plt.subplots()


    Fdr = FTPlot(fig, ax,slider=True)
    Xdata = np.linspace(0, 1)
    Ydata = 5 * np.cos(13 * Xdata)
    Fdr.AddAxis(Name='Axis 1',GridHeight=1,GridPos=1,Unit='m/s')
    Fdr.AddCurve('C1', 'Axis 1', Xdata, Ydata)
    Fdr.AddAxis(Name='Axis 2',GridPos=3,Unit='m/s',Position='Right')
    Fdr.AddCurve('C2', 'Axis 2', Xdata, Ydata, color='r')
    Fdr.AddAxis(Name='Axis 3',GridPos=5,Unit='deg',offset=.1)
    Fdr.AddCurve('C3', 'Axis 3', Xdata, Ydata, color='m')

    Fdr.RemoveCurve('C1')
    Fdr.AddVerticalLine(Name='0.5', Xpos=0.5, linestyle='--')
    Fdr.AddVerticalLine(Name='0.7', Xpos=0.7, linestyle=':', color='green')
    Fdr.AddVerticalLine(Name='0.3', Xpos=0.3, linestyle='-.', marker='s', color='blue')

