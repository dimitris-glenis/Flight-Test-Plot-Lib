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
        if self.slider:
            ext = ax.get_position()
            plt.subplots_adjust(
                bottom=0.15,
                top=ext.y1 + 0.15 - ext.y0
            )

        # Main axes grid setup
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

        # Slider control
        if self.slider:
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

        # Containers for sub-axes and curves
        self.Axis = {}
        self.Curve = {}

        # Connect resize event for auto-rescale
        self.fig.canvas.mpl_connect('resize_event', self._on_resize)

    def updateDataBoxes(self, t):
        for data in self.Curve.values():
            x = data['Curve'].get_xdata()
            y = data['Curve'].get_ydata()
            yv = np.interp(t, x, y)
            vb = data['ValueBox']
            vb.set_position((t, yv))
            vb.set_text(f"{yv:.3f}" if self.XLim[0] <= t <= self.XLim[1] else "")

    def update(self, val):
        # Move slider line and update value boxes
        self.vline.set_xdata([val, val])
        self.updateDataBoxes(val)
        self.fig.canvas.draw_idle()

    def AddAxis(self, Name,
                GridHeight=1, GridPos=0,
                Unit="", Position="Left",
                YLims="Auto", offset=0.02):
        """
        Add a sub-axis GridHeight cells tall, centered on the GridPos-th cell.
        """
        if Name in self.Axis:
            print(f"Axis '{Name}' already exists")
            return

        # Force draw and measure main axes
        self.fig.canvas.draw()
        ext = self.ax.get_position()
        cell_h = ext.height / self.Ngridy
        center = ext.y0 + (GridPos + 0.5) * cell_h
        height = GridHeight * cell_h
        bottom = center - height / 2

        ax2 = self.fig.add_axes([ext.x0, bottom, ext.width, height])

        # Set initial y-limits and enable autoscale
        lims = [-1, 1] if YLims == 'Auto' else YLims
        ax2.set_ylim(lims)
        ax2.set_autoscaley_on(True)

        # Major ticks at bottom, middle, top
        lo, hi = lims
        mid = (lo + hi) / 2
        ax2.set_yticks([lo, mid, hi])

        # Share x-axis
        ax2.set_xlim(self.XLim)
        ax2.set_xticks([])
        ax2.set_xticklabels([])

        # Label and spine styling
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

        # Save metadata
        self.Axis[Name] = {
            'ax': ax2,
            'AutoScale': True,
            'GridHeight': GridHeight,
            'GridPos': GridPos,
            'Position': Position
        }

    def AddCurve(self, Name, Axis, Xdata, Ydata, **kwargs):
        """
        Plot a curve on the specified sub-axis and auto-rescale to integer bounds
        with a one-unit padding beyond any exact-integer data endpoints.
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
        # Value box
        vb = ax2.text(
            Xdata[0], Ydata[0], '',
            ha='center', va='center', color=line.get_color(),
            bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1', alpha=1.0)
        )
        self.Curve[Name] = {'Curve': line, 'ValueBox': vb}

        # Auto-rescale if flagged
        info = self.Axis[Axis]
        if info['AutoScale']:
            dmin = Ydata.min()
            dmax = Ydata.max()
            # Compute integer bounds
            raw_lo = np.floor(dmin)
            raw_hi = np.ceil(dmax)
            # pad by 1 unit if data touches integer
            lo = raw_lo - 1 if np.isclose(dmin, raw_lo) else raw_lo
            hi = raw_hi + 1 if np.isclose(dmax, raw_hi) else raw_hi
            ax2.set_ylim(lo, hi)
            ax2.set_yticks([lo, (lo + hi) / 2, hi])

    def enable_autoscale(self, Name):
        """Re-enable autoscaling for sub-axis 'Name'."""
        if Name not in self.Axis:
            print(f"No such axis '{Name}'")
            return
        info = self.Axis[Name]
        info['AutoScale'] = True
        ax2 = info['ax']
        lines = ax2.get_lines()
        if not lines:
            return
        all_y = np.hstack([ln.get_ydata() for ln in lines])
        dmin = all_y.min()
        dmax = all_y.max()
        raw_lo = np.floor(dmin)
        raw_hi = np.ceil(dmax)
        lo = raw_lo - 1 if np.isclose(dmin, raw_lo) else raw_lo
        hi = raw_hi + 1 if np.isclose(dmax, raw_hi) else raw_hi
        ax2.set_ylim(lo, hi)
        ax2.set_yticks([lo, (lo + hi) / 2, hi])

    def _on_resize(self, event):
        # Reapply padded integer bounds on resize
        for info in self.Axis.values():
            if not info['AutoScale']:
                continue
            ax2 = info['ax']
            lines = ax2.get_lines()
            if not lines:
                continue
            all_y = np.hstack([ln.get_ydata() for ln in lines])
            dmin = all_y.min()
            dmax = all_y.max()
            raw_lo = np.floor(dmin)
            raw_hi = np.ceil(dmax)
            lo = raw_lo - 1 if np.isclose(dmin, raw_lo) else raw_lo
            hi = raw_hi + 1 if np.isclose(dmax, raw_hi) else raw_hi
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
    Fdr.AddAxis(Name='Axis 3',GridPos=5,Unit='deg',offset=.1)




    Fdr.AddCurve('C2','Axis 2',Xdata,Ydata,color='r')
    Fdr.AddCurve('C3','Axis 3',Xdata,Ydata,color='m')

