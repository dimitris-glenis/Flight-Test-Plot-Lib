# -*- coding: utf-8 -*-
"""
Created on Tue May 20 11:51:17 2025

@author: aantonak
"""

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.ticker import AutoMinorLocator, MaxNLocator


class FTPlot:
    """
    Aircraft Simulator for Spin
    """

    def __init__(self, fig, ax, Ngridx=10, Ngridy=13, XLim=(0, 1), slider=False, Xunit="s"):

        self.fig = fig
        self.ax = ax
        self.XLim = XLim  # time axis limits
        self.Ngridx = Ngridx
        self.Ngridy = Ngridy
        self.Xunit = Xunit

        plt.subplots_adjust(left=0.15, right=0.85, top=0.95)

        if slider:

            extent = self.ax.get_position()

            plt.subplots_adjust(
                bottom=0.15, top=extent.y1 + 0.15 - extent.y0
            )  # Make space for slider

        self.extent = ax.get_position()  # x0 = left, y0 = bottom, x1 = right, y1 = top

        if slider:

            # Add initial vertical line
            self.initial_x = XLim[0]  # .5*XLim[0]+.5*XLim[1]
            self.vline = ax.axvline(self.initial_x, color="red", linestyle="--")

            self.slider_ax = self.fig.add_axes(
                [self.extent.x0, 0.02, self.extent.x1 - self.extent.x0, 0.03]
            )
            self.x_slider = Slider(
                self.slider_ax,
                "",
                self.XLim[0],
                self.XLim[1],
                valinit=self.initial_x,
                valfmt="%.2f " + str(self.Xunit),
            )

            self.x_slider.on_changed(self.update)

        self.ax.set_xlim(self.XLim)
        self.ax.set_ylim([0, 1])

        self.ax.set_xticks(np.linspace(self.XLim[0], self.XLim[1], Ngridx + 1))
        self.ax.set_yticks(np.linspace(0, 1, Ngridy + 1))
        self.ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
        self.ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))

        self.ax.set_yticklabels([])  # no y tick labels

        self.ax.set_xlabel("Time ({})".format(self.Xunit))

        self.ax.grid("on", which="major", ls="-", color="k", alpha=0.3)
        self.ax.grid("on", which="minor", ls=":", color="k", alpha=0.3)

        self.Axis = {}
        self.Curve = {}

    def updateDataBoxes(self, x_value):

        for curve in self.Curve.keys():

            # get curve x, y data
            x = self.Curve[curve]["Curve"].get_xdata()
            y = self.Curve[curve]["Curve"].get_ydata()

            y_value = np.interp(x_value, x, y)

            self.Curve[curve]["ValueBox"].set_position((x_value, y_value))

            if x_value > self.ax.get_xlim()[0] and x_value < self.ax.get_xlim()[1]:

                self.Curve[curve]["ValueBox"].set_text("{:.3f}".format(y_value))

            else:

                self.Curve[curve]["ValueBox"].set_text("")

    def update(self, val):

        self.vline.set_xdata([val, val])  # Update x-position of the vertical line
        self.updateDataBoxes(val)
        self.fig.canvas.draw_idle()  # Redraw the plot

    def AddAxis(
        self, Name, GridHeight=4, GridPos=0, Unit="", Position="Left", YLims="Auto", offset=0.02
    ):

        if Name in self.Axis.keys():

            print("Name already exists")

        else:

            self.Axis[Name] = {
                "ax": self.fig.add_axes(
                    [
                        self.extent.x0,
                        self.extent.y0
                        + GridPos / (self.Ngridy - 0) * (self.extent.y1 - self.extent.y0),
                        self.extent.x1 - self.extent.x0,
                        GridHeight / (self.Ngridy - 0) * (self.extent.y1 - self.extent.y0),
                    ]
                ),
                "Name": Name,
                "Unit": Unit,
                "AutoScale": YLims == "Auto",
                "Position": Position,
                "GridHeight": GridHeight,
                "GridPos": GridPos,
            }

            self.Axis[Name]["ax"].set_ylabel(
                self.Axis[Name]["Name"] + "\n" + self.Axis[Name]["Unit"], labelpad=0
            )

            # Adjust y scale

            if YLims == "Auto":

                YLims = [-1, 1]

            self.Axis[Name]["ax"].set_ylim(YLims)

            YRange = (YLims[1] - YLims[0]) / 2
            Ymean = np.mean(YLims)

            self.Axis[Name]["ax"].set_yticks([-0.5 * YRange + Ymean, Ymean, 0.5 * YRange + Ymean])

            self.Axis[Name]["ax"].set_xlim(self.XLim)

            self.Axis[Name]["ax"].set_xticks([])
            self.Axis[Name]["ax"].set_xticklabels([])

            self.Axis[Name]["ax"].spines["bottom"].set_visible(False)
            self.Axis[Name]["ax"].spines["top"].set_visible(False)
            self.Axis[Name]["ax"].patch.set_visible(False)

            if Position == "Left":

                self.Axis[Name]["ax"].spines["left"].set_bounds(
                    -0.5 * YRange + Ymean, 0.5 * YRange + Ymean
                )
                self.Axis[Name]["ax"].spines["right"].set_visible(False)

                self.Axis[Name]["ax"].spines["left"].set_position(("axes", -offset))

            else:

                self.Axis[Name]["ax"].yaxis.tick_right()
                self.Axis[Name]["ax"].yaxis.set_label_position("right")
                self.Axis[Name]["ax"].spines["right"].set_bounds(
                    -0.5 * YRange + Ymean, 0.5 * YRange + Ymean
                )
                self.Axis[Name]["ax"].spines["left"].set_visible(False)

                self.Axis[Name]["ax"].spines["right"].set_position(("axes", 1 + offset))

    def AddCurve(self, Name, Axis, Xdata, Ydata, **kwargs):

        if Name in self.Curve.keys():

            print("Name already exists")

        else:

            (C,) = self.Axis[Axis]["ax"].plot(Xdata, Ydata, **kwargs)

            idx = np.random.randint(len(Xdata) // 4, int(len(Xdata) * 0.75))

            Label = self.Axis[Axis]["ax"].text(
                Xdata[idx],
                Ydata[idx],
                Name,
                ha="center",
                va="center",
                color=C.get_color(),
                bbox=dict(facecolor="w", edgecolor="w", boxstyle="round,pad=0.1"),
            )

            ValueBox = self.Axis[Axis]["ax"].text(
                Xdata[0],
                Ydata[0],
                "",
                ha="center",
                va="center",
                color=C.get_color(),
                bbox=dict(facecolor="w", edgecolor="w", boxstyle="round,pad=0.1", alpha=1.0),
            )

            self.Curve[Name] = {"Curve": C, "Axis": Axis, "Label": Label, "ValueBox": ValueBox}

            self.UpdateAxis(Axis)

    def UpdateAxis(self, Axis):

        if self.Axis[Axis]["AutoScale"]:

            self.Axis[Axis]["ax"].set_autoscaley_on(True)
            self.Axis[Axis]["ax"].relim()

            Nbins = self.Axis[Axis]["GridHeight"] + 1

            self.Axis[Axis]["ax"].yaxis.set_major_locator(MaxNLocator(nbins=Nbins))
            self.Axis[Axis]["ax"].autoscale_view()

            yticks = self.Axis[Axis]["ax"].get_yticks()

            self.Axis[Axis]["ax"].set_ylim([yticks[0], yticks[-1]])

            YLims = self.Axis[Axis]["ax"].get_ylim()
            YRange = (YLims[1] - YLims[0]) / 2
            Ymean = np.mean(YLims)

            self.Axis[Axis]["ax"].set_yticks([-0.5 * YRange + Ymean, Ymean, 0.5 * YRange + Ymean])

            if self.Axis[Axis]["Position"] == "Left":

                self.Axis[Axis]["ax"].spines["left"].set_bounds(
                    -0.5 * YRange + Ymean, 0.5 * YRange + Ymean
                )

            else:

                self.Axis[Axis]["ax"].spines["right"].set_bounds(
                    -0.5 * YRange + Ymean, 0.5 * YRange + Ymean
                )


"""
plt.close('all')


fig,ax = plt.subplots()


Fdr = FDRplot(ax,slider=True)

Fdr.AddAxis(Name='Axis 1',GridHeight=6,GridPos=.5,Unit='m/s')
Fdr.AddAxis(Name='Axis 2',GridPos=3,Unit='m/s',Position='Right')
Fdr.AddAxis(Name='Axis 3',GridPos=5,Unit='deg',offset=.1)

Xdata = np.linspace(0,1)
Ydata = 3.7 * np.cos(13 * Xdata)

Fdr.AddCurve('C1','Axis 1',Xdata,Ydata)
Fdr.AddCurve('C2','Axis 2',Xdata,Ydata,color='r')
Fdr.AddCurve('C3','Axis 3',Xdata,Ydata,color='m')
"""
