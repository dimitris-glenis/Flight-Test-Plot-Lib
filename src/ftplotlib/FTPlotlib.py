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

        # Reserve space
        plt.subplots_adjust(left=0.15, right=0.85, top=0.95)
        if slider:
            ext = ax.get_position()
            plt.subplots_adjust(bottom=0.15, top=ext.y1 + 0.15 - ext.y0)

        # Main grid
        ax.set_xlim(XLim)
        ax.set_ylim(0, 1)
        ax.set_xticks(np.linspace(XLim[0], XLim[1], Ngridx + 1))
        ax.set_yticks(np.linspace(0, 1, Ngridy + 1))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.set_yticklabels([])
        ax.set_xlabel(f"Time ({Xunit})")
        ax.grid(True, which='major', ls='-', color='k', alpha=0.3)
        ax.grid(True, which='minor', ls=':', color='k', alpha=0.3)

        # Slider
        if slider:
            self.initial_x = XLim[0]
            self.vline = ax.axvline(self.initial_x, color='red', linestyle='--')
            ext = ax.get_position()
            self.slider_ax = fig.add_axes([ext.x0, 0.02, ext.width, 0.03])
            self.x_slider = Slider(self.slider_ax, "", XLim[0], XLim[1],
                                   valinit=self.initial_x,
                                   valfmt=f"%.2f {Xunit}")
            self.x_slider.on_changed(self.update)

        self.Axis = {}
        self.Curve = {}
        self.fig.canvas.mpl_connect('resize_event', self._on_resize)

    def updateDataBoxes(self, t):
        for data in self.Curve.values():
            x, y = data['Curve'].get_xdata(), data['Curve'].get_ydata()
            yv = np.interp(t, x, y)
            vb = data['ValueBox']
            vb.set_position((t, yv))
            vb.set_text(f"{yv:.3f}" if self.XLim[0] <= t <= self.XLim[1] else "")

    def update(self, val):
        self.vline.set_xdata([val, val])
        self.updateDataBoxes(val)
        self.fig.canvas.draw_idle()

    def AddAxis(self, Name, GridHeight=1, GridPos=0, Unit="", Position="Left", YLims="Auto", offset=0.02):
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

        # initial limits
        lims = [-1,1] if YLims=='Auto' else YLims
        ax2.set_ylim(lims)
        ax2.set_autoscaley_on(True)
        lo, hi = lims
        ax2.set_yticks([lo, (lo+hi)/2, hi])

        ax2.set_xlim(self.XLim)
        ax2.set_xticks([]); ax2.set_xticklabels([])
        ax2.set_ylabel(f"{Name}\n{Unit}", labelpad=0)
        for s in ('top','bottom'): ax2.spines[s].set_visible(False)
        ax2.patch.set_visible(False)
        if Position.lower()=='left':
            ax2.spines['right'].set_visible(False)
            ax2.spines['left'].set_position(('axes',-offset))
        else:
            ax2.yaxis.tick_right(); ax2.yaxis.set_label_position('right')
            ax2.spines['left'].set_visible(False)
            ax2.spines['right'].set_position(('axes',1+offset))

        self.Axis[Name] = {'ax':ax2,'AutoScale':True,'GridHeight':GridHeight,'GridPos':GridPos,'Position':Position}

    def AddCurve(self, Name, Axis, Xdata, Ydata, **kwargs):
        ax2 = self.Axis[Axis]['ax']
        line,_ = ax2.plot(Xdata, Ydata, **kwargs)
        idx = len(Xdata)//2
        ax2.text(Xdata[idx],Ydata[idx],Name,ha='center',va='center',color=line.get_color(),bbox=dict(facecolor='white',edgecolor='white',boxstyle='round,pad=0.1'))
        vb = ax2.text(Xdata[0],Ydata[0],'',ha='center',va='center',color=line.get_color(),bbox=dict(facecolor='white',edgecolor='white',boxstyle='round,pad=0.1',alpha=1.0))
        self.Curve[Name]={'Curve':line,'ValueBox':vb}

        info=self.Axis[Axis]
        if info['AutoScale']:
            dmin,dmax=Ydata.min(),Ydata.max()
            raw_lo,raw_hi=np.floor(dmin),np.ceil(dmax)
            lo=raw_lo-1 if np.isclose(dmin,raw_lo) else raw_lo
            hi=raw_hi+1 if np.isclose(dmax,raw_hi) else raw_hi
            ax2.set_ylim(lo,hi)
            ax2.set_yticks([lo,(lo+hi)/2,hi])

    def enable_autoscale(self,Name):
        if Name not in self.Axis: return
        info=self.Axis[Name]
        info['AutoScale']=True
        ax2=info['ax']
        ys=np.hstack([ln.get_ydata() for ln in ax2.get_lines()])
        raw_lo,raw_hi=np.floor(ys.min()),np.ceil(ys.max())
        lo=raw_lo-1 if np.isclose(ys.min(),raw_lo) else raw_lo
        hi=raw_hi+1 if np.isclose(ys.max(),raw_hi) else raw_hi
        ax2.set_ylim(lo,hi)
        ax2.set_yticks([lo,(lo+hi)/2,hi])

    def _on_resize(self,event):
        for info in self.Axis.values():
            if not info['AutoScale']: continue
            ax2=info['ax']
            ys=np.hstack([ln.get_ydata() for ln in ax2.get_lines()])
            raw_lo,raw_hi=np.floor(ys.min()),np.ceil(ys.max())
            lo=raw_lo-1 if np.isclose(ys.min(),raw_lo) else raw_lo
            hi=raw_hi+1 if np.isclose(ys.max(),raw_hi) else raw_hi
            ax2.set_ylim(lo,hi)
            ax2.set_yticks([lo,(lo+hi)/2,hi])



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

