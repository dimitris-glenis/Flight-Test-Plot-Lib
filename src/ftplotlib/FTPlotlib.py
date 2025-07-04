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
                 Ngridx=10, Ngridy=10,
                 XLim=(0, 1), slider=False, Xunit="s"):
        # Initialize core attributes
        self.fig = fig
        self.ax = ax
        self.XLim = XLim
        self.Ngridx = Ngridx
        self.Ngridy = Ngridy
        self.Xunit = Xunit

        # Adjust margins for grid and slider space
        plt.subplots_adjust(left=0.15, right=0.85, top=0.95)
        if slider:
            ext = ax.get_position()
            plt.subplots_adjust(
                bottom=0.15,
                top=ext.y1 + 0.15 - ext.y0
            )

        # Configure main axes grid
        self.ax.set_xlim(self.XLim)
        self.ax.set_ylim([0, 1])
        self.ax.set_xticks(np.linspace(self.XLim[0], self.XLim[1], Ngridx + 1))
        self.ax.set_yticks(np.linspace(0, 1, Ngridy + 1))
        self.ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.set_yticklabels([])
        self.ax.set_xlabel(f"Time ({self.Xunit})")
        self.ax.grid(True, which="major", ls="-", color="k", alpha=0.3)
        self.ax.grid(True, which="minor", ls=":", color="k", alpha=0.3)

        # Slider setup
        self.slider = slider
        if slider:
            self.initial_x = XLim[0]
            self.vline = self.ax.axvline(self.initial_x, color="red", linestyle="--")
            pos = self.ax.get_position()
            self.slider_ax = self.fig.add_axes([
                pos.x0,
                0.02,
                pos.width,
                0.03
            ])
            self.x_slider = Slider(
                self.slider_ax, "",
                self.XLim[0], self.XLim[1],
                valinit=self.initial_x,
                valfmt="%.2f " + self.Xunit
            )
            self.x_slider.on_changed(self.update)

        # Containers for extra axes and curves
        self.Axis = {}
        self.Curve = {}

    def updateDataBoxes(self, x_value):
        for data in self.Curve.values():
            x = data['Curve'].get_xdata()
            y = data['Curve'].get_ydata()
            yv = np.interp(x_value, x, y)
            vb = data['ValueBox']
            vb.set_position((x_value, yv))
            if self.XLim[0] <= x_value <= self.XLim[1]:
                vb.set_text(f"{yv:.3f}")
            else:
                vb.set_text("")

    def update(self, val):
        # Move vertical line and update value boxes
        self.vline.set_xdata([val, val])
        self.updateDataBoxes(val)
        self.fig.canvas.draw_idle()

    def AddAxis(self, Name,
                GridHeight=1, GridPos=0,
                Unit="", Position="Left",
                YLims="Auto", offset=0.02):
        """
        Adds a sub-axis GridHeight data-grid cells tall,
        centered on the GridPos'th cell (0-based).
        GridHeight=1 -> one full rectangle high.
        GridPos=k -> axis center at middle of k-th rectangle.
        """
        if Name in self.Axis:
            print(f"Axis '{Name}' already exists")
            return

        # Ensure canvas is up-to-date before measuring
        self.fig.canvas.draw()

        # Get current main axes position (in figure coords)
        ext = self.ax.get_position()  # Bbox([x0, y0, x1, y1])
        total_h = ext.y1 - ext.y0
        cell_h = total_h / self.Ngridy

        # Compute bottom so that the axis center lies at (GridPos + 0.5) * cell_h
        bottom = ext.y0 + (GridPos + 0.5 - GridHeight/2) * cell_h
        height = GridHeight * cell_h

        # Create the new axes at the computed position
        ax2 = self.fig.add_axes([
            ext.x0,
            bottom,
            ext.width,
            height
        ])

        # Label & styling
        ax2.set_ylabel(f"{Name}\n{Unit}", labelpad=0)
        auto = (YLims == "Auto")
        lims = [-1, 1] if auto else YLims
        ax2.set_ylim(lims)

        # Three ticks: bottom, middle, top
        yr = (lims[1] - lims[0]) / 2
        ym = np.mean(lims)
        ax2.set_yticks([ym - yr, ym, ym + yr])

        # Share X-axis limits and hide labels
        ax2.set_xlim(self.XLim)
        ax2.set_xticks([])
        ax2.set_xticklabels([])

        # Remove top/bottom spines and background
        for s in ('top', 'bottom'):
            ax2.spines[s].set_visible(False)
        ax2.patch.set_visible(False)

        # Position left or right spine with offset
        if Position.lower() == 'left':
            ax2.spines['right'].set_visible(False)
            ax2.spines['left'].set_bounds(ym - yr, ym + yr)
            ax2.spines['left'].set_position(('axes', -offset))
        else:
            ax2.yaxis.tick_right()
            ax2.yaxis.set_label_position('right')
            ax2.spines['left'].set_visible(False)
            ax2.spines['right'].set_bounds(ym - 0.5*yr, ym + 0.5*yr)
            ax2.spines['right'].set_position(('axes', 1 + offset))

        # Store reference for future curves & autoscale
        self.Axis[Name] = {
            'ax': ax2,
            'Name': Name,
            'Unit': Unit,
            'AutoScale': auto,
            'Position': Position,
            'GridHeight': GridHeight,
            'GridPos': GridPos
        }

    def AddCurve(self, Name, Axis, Xdata, Ydata, **kwargs):
        if Name in self.Curve:
            print(f"Curve '{Name}' already exists")
            return
        ax2 = self.Axis[Axis]['ax']
        line, = ax2.plot(Xdata, Ydata, **kwargs)

        # Static label in the middle
        idx = np.random.randint(len(Xdata)//4, int(len(Xdata)*0.75))
        ax2.text(
            Xdata[idx], Ydata[idx], Name,
            ha='center', va='center', color=line.get_color(),
            bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1')
        )

        # Dynamic value box at start
        vb = ax2.text(
            Xdata[0], Ydata[0], '',
            ha='center', va='center', color=line.get_color(),
            bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1', alpha=1.0)
        )
        self.Curve[Name] = {'Curve': line, 'Axis': Axis, 'ValueBox': vb}

        # Auto-scale if needed
        if self.Axis[Axis]['AutoScale']:
            self.UpdateAxis(Axis)

    def UpdateAxis(self, Axis):
        info = self.Axis[Axis]
        ax2 = info['ax']
        if not info['AutoScale']:
            return
        ax2.set_autoscaley_on(True)
        ax2.relim()
        ax2.autoscale_view()

        nb = info['GridHeight'] + 1
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=nb))
        yt = ax2.get_yticks()
        ax2.set_ylim(yt[0], yt[-1])

        # Recompute three-tick scheme
        YL = ax2.get_ylim()
        rng = (YL[1] - YL[0]) / 2
        mid = np.mean(YL)
        ax2.set_yticks([mid - 0.5*rng, mid, mid + 0.5*rng])

        # Adjust spine length
        spine = 'left' if info['Position'].lower()=='left' else 'right'
        ax2.spines[spine].set_bounds(mid - 0.5*rng, mid + 0.5*rng)


if __name__ == "__main__":

    matplotlib.use('Qt5Agg')
    plt.ion()
    plt.close('all')
    fig, ax = plt.subplots()


    Fdr = FTPlot(fig, ax,slider=True)

    Fdr.AddAxis(Name='Axis 1',GridHeight=1,GridPos=1,Unit='m/s')
    Fdr.AddAxis(Name='Axis 2',GridPos=3,Unit='m/s',Position='Right')
    Fdr.AddAxis(Name='Axis 3',GridPos=5,Unit='deg',offset=.1)

    Xdata = np.linspace(0,1)
    Ydata = 3.7 * np.cos(13 * Xdata)

    Fdr.AddCurve('C1','Axis 1',Xdata,Ydata)
    Fdr.AddCurve('C2','Axis 2',Xdata,Ydata,color='r')
    Fdr.AddCurve('C3','Axis 3',Xdata,Ydata,color='m')

