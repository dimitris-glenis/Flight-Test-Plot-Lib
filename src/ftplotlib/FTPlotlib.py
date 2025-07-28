# -*- coding: utf-8 -*-
"""
Created on Tue May 20 11:51:17 2025

@author: aantonak
"""
import textwrap

import matplotlib
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.ticker import AutoMinorLocator, MaxNLocator


# pip install --upgrade --force-reinstall -e "git+ssh://git@github.com/dimitris-glenis/Flight-Test-Plot-Lib.git@pip-installable#egg=ftplotlib"

class FTPlot:
    """
    Aircraft Simulator for Spin
    """
    def __init__(self, fig, ax,
                 Ngridx=10, Ngridy=13,
                 XLim=(0, 1), slider=False, Xunit="s"):
        # core references
        self.fig = fig
        self.ax = ax
        self.Ngridx = Ngridx
        self.Ngridy = Ngridy
        self.XLim = XLim
        self.Xunit = Xunit
        self.slider = slider

        # adjust main axes
        plt.subplots_adjust(left=0.15, right=0.85, top=0.95)
        if slider:
            ext = ax.get_position()
            plt.subplots_adjust(bottom=0.15, top=ext.y1 + 0.15 - ext.y0)

        # set up main grid
        self.ax.set_xlim(XLim)
        self.ax.set_ylim(0, 1)
        self.ax.set_xticks(np.linspace(XLim[0], XLim[1], Ngridx + 1))
        self.ax.set_yticks(np.linspace(0, 1, Ngridy + 1))
        self.ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.set_yticklabels([])
        self.ax.set_xlabel(f"Time ({Xunit})")
        self.ax.grid(True, which='major', ls='-', color='k', alpha=0.3)
        self.ax.grid(True, which='minor', ls=':', color='k', alpha=0.3)

        # optional slider
        if slider:
            self.initial_x = XLim[0]
            self.vline = ax.axvline(self.initial_x, color='red', linestyle='--')
            ext = ax.get_position()
            self.slider_ax = fig.add_axes([ext.x0, 0.02, ext.width, 0.03])
            self.x_slider = Slider(self.slider_ax, "",
                                   XLim[0], XLim[1],
                                   valinit=self.initial_x,
                                   valfmt=f"%.2f {Xunit}")
            self.x_slider.on_changed(self.update)

        # containers
        self.Axis = {}
        self.Curve = {}
        self.fig.canvas.mpl_connect('resize_event', self._on_resize)
        # setup interactive label dragging
        self._dragging = None
        self.fig.canvas.mpl_connect('pick_event', self._on_pick)
        self.fig.canvas.mpl_connect('motion_notify_event', self._on_motion)
        self.fig.canvas.mpl_connect('button_release_event', self._on_release)

    def updateDataBoxes(self, t):
        for entry in self.Curve.values():
            x, y = entry['Curve'].get_xdata(), entry['Curve'].get_ydata()
            yv = np.interp(t, x, y)
            vb = entry['ValueBox']
            vb.set_position((t, yv))
            vb.set_text(f"{yv:.3f}" if self.XLim[0] <= t <= self.XLim[1] else "")

    def update(self, val):
        self.vline.set_xdata([val, val])
        self.updateDataBoxes(val)
        self.fig.canvas.draw_idle()

    def AddAxis(self, Name, GridHeight=1, GridPos=0,
                Unit="", Position="Left", YLims="Auto", offset=0.02):
        """
        Create a sub-axis spanning GridHeight cells,
        centered on the GridPos-th grid cell.
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

        # initial limits
        lims = (0, 1) if YLims == 'Auto' else YLims
        ax2.set_ylim(lims)
        ax2.set_autoscaley_on(True)
        lo, hi = lims
        ax2.set_yticks([lo, (lo + hi) / 2, hi])

        # share x
        ax2.set_xlim(self.XLim)
        ax2.set_xticks([]); ax2.set_xticklabels([])
        ax2.set_ylabel(f"{Name}\n{Unit}", labelpad=0)
        for s in ('top', 'bottom'): ax2.spines[s].set_visible(False)
        ax2.patch.set_visible(False)
        if Position.lower() == 'left':
            ax2.spines['right'].set_visible(False)
            ax2.spines['left'].set_position(('axes', -offset))
        else:
            ax2.yaxis.tick_right()
            ax2.yaxis.set_label_position('right')
            ax2.spines['left'].set_visible(False)
            ax2.spines['right'].set_position(('axes', 1 + offset))

        self.Axis[Name] = dict(ax=ax2, GridHeight=GridHeight,
                               GridPos=GridPos, Position=Position,
                               AutoScale=True)

    def AddCurve(self, Name, Axis, Xdata, Ydata, **kwargs):
        """
        Plot on sub-axis; auto-rescale y-limits to include all curves on that axis,
        but keep axis height constant.
        """
        ax2 = self.Axis[Axis]['ax']
        line, = ax2.plot(Xdata, Ydata, **kwargs)
        lbl = ax2.text(
            Xdata[len(Xdata)//2], Ydata[len(Ydata)//2], Name,
            ha='center', va='center', color=line.get_color(),
            bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1')
        )
        lbl.set_picker(True)
        vb = ax2.text(
            Xdata[0], Ydata[0], '',
            ha='center', va='center', color=line.get_color(),
            bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1', alpha=1.0)
        )
        self.Curve[Name] = dict(Curve=line, Label=lbl, ValueBox=vb)

        # recompute y-limits across all curves on this axis
        info = self.Axis[Axis]
        if info['AutoScale']:
            y_all = np.hstack([ln.get_ydata() for ln in ax2.get_lines()])
            dmin, dmax = y_all.min(), y_all.max()
            lo, hi = np.floor(dmin), np.ceil(dmax)
            ax2.set_ylim(lo, hi)
            ax2.set_yticks([lo, (lo + hi) / 2, hi])

    def RemoveCurve(self, Name):
        """
        Remove a previously added curve, its label, and its value box.
        """
        if Name not in self.Curve:
            print(f"No such curve '{Name}'")
            return
        entry = self.Curve.pop(Name)
        for art in entry.values():
            try:
                art.remove()
            except Exception:
                pass
        self.fig.canvas.draw_idle()

    def RemoveAxis(self, Name):
        if Name not in self.Axis: return
        ax2 = self.Axis[Name]['ax']
        for c, ent in list(self.Curve.items()):
            if ent['Curve'].axes is ax2:
                self.RemoveCurve(c)
        ax2.remove()
        del self.Axis[Name]
        self.fig.canvas.draw_idle()

    def SetAxisLimits(self, Name, lo, hi):
        """
        Manually set the y-axis limits for sub-axis 'Name'.
        """
        if Name not in self.Axis:
            print(f"No such axis '{Name}'")
            return
        ax2 = self.Axis[Name]['ax']
        ax2.set_ylim(lo, hi)
        ax2.set_yticks([lo, (lo + hi) / 2, hi])
        self.fig.canvas.draw_idle()

    def SetAxisHeight(self, Name, GridHeight):
        """
        Manually adjust the height (number of grid cells) of sub-axis 'Name'.
        """
        if Name not in self.Axis:
            print(f"No such axis '{Name}'")
            return
        info = self.Axis[Name]
        info['GridHeight'] = GridHeight
        self._resize_axis(Name)
        self.fig.canvas.draw_idle()

    def SetLabelPosition(self, Name, x, y):
        """
        Manually reposition the text label for curve 'Name'.
        """
        if Name not in self.Curve:
            print(f"No such curve '{Name}'")
            return
        lbl = self.Curve[Name]['Label']
        lbl.set_position((x, y))
        self.fig.canvas.draw_idle()

    def enable_autoscale(self, Name):
        if Name not in self.Axis: return
        self.Axis[Name]['AutoScale'] = True
        ax2 = self.Axis[Name]['ax']
        y_all = np.hstack([ln.get_ydata() for ln in ax2.get_lines()])
        lo, hi = np.floor(y_all.min()), np.ceil(y_all.max())
        ax2.set_ylim(lo, hi)

    def AddVerticalLine(self, Name, Xpos, YLims=None, linestyle='--', marker=None, **kwargs):
        if Name in self.Curve:
            print(f"Vertical line '{Name}' already exists")
            return
        if YLims is None:
            YLims = self.ax.get_ylim()
        if 'linestyle' not in kwargs and 'ls' not in kwargs:
            kwargs['linestyle'] = linestyle
        line = self.ax.axvline(Xpos, **kwargs)
        vb = self.ax.text(
            Xpos, YLims[0] - 0.05, f"{Xpos:.3f}",
            ha='center', va='top', color=line.get_color(),
            bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1', alpha=1.0)
        )
        intersections = {}
        self.Curve[Name] = {'Curve': line, 'ValueBox': vb, 'IsVerticalLine': True, 'XPosition': Xpos, 'Intersections': intersections}
        for curve_name, data in self.Curve.items():
            if curve_name == Name or data.get('IsVerticalLine', False):
                continue
            try:
                curve = data['Curve']
                x_data = curve.get_xdata()
                y_data = curve.get_ydata()
                if Xpos >= x_data.min() and Xpos <= x_data.max():
                    y_val = np.interp(Xpos, x_data, y_data)
                    mark = None
                    if marker is not None:
                        mark = curve.axes.plot([Xpos], [y_val], marker, color=curve.get_color(), markersize=6, zorder=10)[0]
                    ib = curve.axes.text(
                        Xpos, y_val, f"{y_val:.3f}",
                        ha='left', va='bottom', color=curve.get_color(),
                        bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1', alpha=1.0)
                    )
                    intersections[curve_name] = {'marker': mark, 'text': ib}
            except:
                pass
        self.fig.canvas.draw_idle()

    def _resize_axis(self, Name):
        info = self.Axis[Name]
        ext = self.ax.get_position()
        cell_h = ext.height / self.Ngridy
        center = ext.y0 + (info['GridPos'] + 0.5)*cell_h
        height = info['GridHeight'] * cell_h
        bottom = center - height/2
        info['ax'].set_position([ext.x0, bottom, ext.width, height])

    def _on_pick(self, event):
        """Start dragging a picked label."""
        if isinstance(event.artist, plt.Text):
            self._dragging = event.artist

    def _on_motion(self, event):
        """Drag the picked label with the mouse."""
        if self._dragging and event.inaxes == self._dragging.axes:
            self._dragging.set_position((event.xdata, event.ydata))
            self.fig.canvas.draw_idle()

    def _on_release(self, event):
        """Release the dragged label."""
        self._dragging = None
        for Name, info in self.Axis.items():
            if info['AutoScale']:
                ax2 = info['ax']
                y_all = np.hstack([ln.get_ydata() for ln in ax2.get_lines()])
                lo, hi = np.floor(y_all.min()), np.ceil(y_all.max())
                ax2.set_ylim(lo, hi)
                ax2.set_yticks([lo, (lo + hi) / 2, hi])

    def AddFooter(self, info_dict, fontsize=8, color='gray', pad=0.5, wrap_width=80):
        """
        Place a multi‑line footer immediately under the x‑axis LABEL.
        - info_dict: metadata dict (one key per line)
        - fontsize, color: styling
        - pad: fraction of the label's height to separate the footer
        - wrap_width: target wrap chars × fraction of axis width
        """
        # Store footer info for resize events
        self._footer_metadata = (info_dict, fontsize, color, pad, wrap_width)

        # Remove any old footer texts
        for txt in getattr(self, '_footer_texts', []):
            txt.remove()
        self._footer_texts = []

        # Force draw to get accurate measurements
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()

        # Get x-axis label position
        xlabel = self.ax.xaxis.label
        lbl_box = xlabel.get_window_extent(renderer)
        lbl_box_fig = lbl_box.transformed(self.fig.transFigure.inverted())

        # Format and wrap metadata text
        lines = [f"{k}: {v}" for k, v in info_dict.items()]
        axis_w = self.ax.get_position().width
        wrap_chars = max(10, int(wrap_width * axis_w))
        wrapped = "\n".join(textwrap.fill(ln, wrap_chars) for ln in lines)

        # Calculate vertical position with small padding
        pad_fig = pad * lbl_box_fig.height
        x_fig = 0.5 * (lbl_box_fig.x0 + lbl_box_fig.x1)
        y_fig = lbl_box_fig.y0 - pad_fig

        # Add footer text
        footer = self.fig.text(
            x_fig, y_fig, wrapped,
            transform=self.fig.transFigure,
            ha='center', va='top',
            fontsize=fontsize,
            color=color,
            family='monospace',
            linespacing=1.2
        )
        footer.footer_tag = True
        self._footer_texts.append(footer)

        # Get footer dimensions
        self.fig.canvas.draw()
        fbbox = footer.get_window_extent(renderer)
        fbbox_fig = fbbox.transformed(self.fig.transFigure.inverted())

        # Calculate slider position
        slider_height = 0.03
        slider_padding = 0.01
        min_bottom = 0.01

        # Calculate total required space
        if self.slider:
            min_required_y = min_bottom + slider_height + slider_padding
        else:
            min_required_y = min_bottom

        # If footer extends below minimum height, adjust layout
        if fbbox_fig.y0 < min_required_y:
            # Calculate needed adjustment
            needed = min_required_y - fbbox_fig.y0

            # Store original positions
            orig_positions = {}
            for name, axis_info in self.Axis.items():
                orig_positions[name] = axis_info['ax'].get_position().bounds
            main_pos = self.ax.get_position().bounds

            # Adjust figure's bottom margin
            old_bottom = self.fig.subplotpars.bottom
            new_bottom = old_bottom + needed
            plt.subplots_adjust(bottom=new_bottom)

            # Update main axis position
            scale_factor = (1.0 - new_bottom) / (1.0 - old_bottom)
            new_height = main_pos[3] * scale_factor
            self.ax.set_position([main_pos[0], new_bottom, main_pos[2], new_height])

            # Update all sub-axes proportionally
            for name, (x0, y0, width, height) in orig_positions.items():
                rel_y = (y0 - main_pos[1]) / main_pos[3]  # Relative position
                rel_height = height / main_pos[3]  # Relative height

                new_y = new_bottom + (rel_y * new_height)
                new_height_sub = rel_height * new_height

                self.Axis[name]['ax'].set_position([x0, new_y, width, new_height_sub])

            # Redraw and reposition footer with new coordinates
            self.fig.canvas.draw()
            lbl_box = xlabel.get_window_extent(renderer)
            lbl_box_fig = lbl_box.transformed(self.fig.transFigure.inverted())
            y_fig = lbl_box_fig.y0 - pad_fig
            footer.set_position((x_fig, y_fig))

        # Position slider if present
        if self.slider:
            # Position slider just below footer
            fbbox = footer.get_window_extent(renderer)
            fbbox_fig = fbbox.transformed(self.fig.transFigure.inverted())

            slider_y = fbbox_fig.y0 - slider_padding - slider_height
            slider_y = max(min_bottom, slider_y)  # Ensure minimum bottom margin

            ext = self.ax.get_position()
            self.slider_ax.set_position([ext.x0, slider_y, ext.width, slider_height])

        self.fig.canvas.draw_idle()

    def _on_resize(self, event):
        """Handle resize events by updating all axes and footer"""
        # Store current size
        old_width = getattr(self, '_last_width', event.width)
        old_height = getattr(self, '_last_height', event.height)

        # Force initial draw to get accurate measurements
        self.fig.canvas.draw()

        # First, update all axis positions
        for name, info in self.Axis.items():
            self._resize_axis(name)

        # Re-autoscale axes that need it
        for info in self.Axis.values():
            if info.get('AutoScale', False):
                ax2 = info['ax']
                all_y = np.hstack([ln.get_ydata() for ln in ax2.get_lines()])
                if len(all_y) > 0:  # Only adjust if there's data
                    lo, hi = np.floor(all_y.min()), np.ceil(all_y.max())
                    ax2.set_ylim(lo, hi)
                    ax2.set_yticks([lo, 0.5 * (lo + hi), hi])

        # Reset margins when window size changes significantly
        size_changed = (abs(event.width - old_width) > 10 or
                       abs(event.height - old_height) > 10)

        if size_changed:
            # Reset to default margins
            plt.subplots_adjust(bottom=0.1, top=0.9)

        # Store current size for next comparison
        self._last_width = event.width
        self._last_height = event.height

        # Restore title if it exists
        if hasattr(self, '_title_metadata'):
            title, fontsize, pad, kwargs = self._title_metadata
            self.ax.set_title(title, fontsize=fontsize, pad=pad, **kwargs)

            # Force draw to get accurate measurements
            self.fig.canvas.draw()
            renderer = self.fig.canvas.get_renderer()

            # Get title dimensions and adjust top margin
            bbox = self.ax.title.get_window_extent(renderer)
            bbox_fig = bbox.transformed(self.fig.transFigure.inverted())

            # Add extra padding to ensure title is visible
            extra_padding = 0.03
            required_top = 1.0 - bbox_fig.height - extra_padding

            # Always adjust top margin for title
            plt.subplots_adjust(top=required_top)

        # Force draw before recalculating footer position
        self.fig.canvas.draw()

        # Update the footer with proper layout (do this last)
        if hasattr(self, '_footer_metadata'):
            self.AddFooter(*self._footer_metadata)


    def AddTitle(self, title, fontsize=12, pad=10, **kwargs):
        """
        Add a title to the plot.

        Parameters:
        -----------
        title : str
            The title text
        fontsize : int
            Font size for the title
        pad : float
            Padding between the title and the plot
        **kwargs : dict
            Additional keyword arguments passed to plt.title()
        """
        # Set default values for some parameters if not provided
        if 'fontweight' not in kwargs:
            kwargs['fontweight'] = 'bold'

        # Remove any existing title
        if hasattr(self.ax, 'title') and self.ax.title:
            self.ax.title.set_text("")

        # Add title to the main axis
        title_obj = self.ax.set_title(title, fontsize=fontsize, pad=pad, **kwargs)

        # Force draw to get accurate measurements
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()

        # Get title dimensions
        bbox = title_obj.get_window_extent(renderer)
        bbox_fig = bbox.transformed(self.fig.transFigure.inverted())

        # Calculate required top margin with extra padding
        extra_padding = 0.03  # Add more space to ensure visibility
        required_top = 1.0 - bbox_fig.height - extra_padding

        # Ensure top margin is sufficient
        plt.subplots_adjust(top=required_top)

        # Store title for resize events
        self._title_metadata = (title, fontsize, pad, kwargs)

        # Redraw the figure
        self.fig.canvas.draw_idle()


if __name__ == "__main__":

    # matplotlib.use('Qt5Agg')
    plt.ion()
    plt.close('all')
    fig, ax = plt.subplots()

    metadata = {'aero_db': 'PDR_CAND4_Wg009_Fl001_ST007_TT064_Fu005_Bf017_Vt004_Ht004_Na003_Py016_HQ',
     'thrust_model': '20250528_throttle_to_thrust',
     'aircraft_characteristics': 'PDR_cand4_v1.yaml',
     'version': '1.5.0',
     'branch': '100-new-main',
     'commit': '634dae3',
     'run_date': '2025-07-22T15:03:05.206297'}

    Fdr = FTPlot(fig, ax,slider=True)
    Fdr.AddTitle("Flight Test Data", fontsize=14, color='navy')
    Fdr.AddFooter(metadata)
    Xdata = np.linspace(0, 1)
    Ydata = 5 * np.cos(13 * Xdata)
    Fdr.AddAxis(Name='Axis 1',GridHeight=1,GridPos=1,Unit='m/s')
    Fdr.AddCurve('C1', 'Axis 1', Xdata, Ydata)
    Fdr.AddAxis(Name='Axis 2',GridPos=3,Unit='m/s',Position='Right')
    Fdr.AddCurve('C2', 'Axis 2', Xdata, Ydata, color='r')
    Fdr.AddAxis(Name='Axis 3',GridPos=5,Unit='deg',offset=.1)
    Fdr.AddCurve('C3', 'Axis 3', Xdata, Ydata, color='m')

    Fdr.RemoveCurve('C1')
    Fdr.RemoveAxis('Axis 2')
    Fdr.SetAxisLimits('Axis 3', -6, 6)
    Fdr.SetAxisHeight('Axis 3', 3)
    Fdr.AddVerticalLine(Name='0.5', Xpos=0.5, linestyle='--')
    Fdr.AddVerticalLine(Name='0.7', Xpos=0.7, linestyle=':', color='green')
    Fdr.AddVerticalLine(Name='0.3', Xpos=0.3, linestyle='-.', marker='s', color='blue')
    Fdr.SetLabelPosition('C2', 0.5, 0.8)

