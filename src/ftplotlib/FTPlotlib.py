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
        fig.set_constrained_layout(False)
        self.fig = fig
        self.ax = ax
        self.Ngridx = Ngridx
        self.Ngridy = Ngridy
        self.XLim = XLim
        self.Xunit = Xunit
        self.slider = slider

        # adjust main axes
        # plt.subplots_adjust(left=0.15, right=0.85, top=0.95)
        # if slider:
        #     ext = ax.get_position()
        #     plt.subplots_adjust(bottom=0.15, top=ext.y1 + 0.15 - ext.y0)

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
          # 1) draw the slider line
            self.vline = ax.axvline(self.initial_x, color='red', linestyle='--')

          # 2) create its ValueBox up front (hidden at t=0)
            y_min = self.ax.get_ylim()[0]
            self.slider_vb = ax.text(
                self.initial_x, y_min, "",
                ha = "center", va = "top", color = "red",
                bbox = dict(facecolor="white", edgecolor="white",
                boxstyle = "round,pad=0.1", alpha = 1.0),
                visible = False
            )

          # 3) add the Matplotlib slider widget
            ext = ax.get_position()
            self.slider_ax = fig.add_axes([ext.x0, 0.02, ext.width, 0.03])
            self.x_slider = Slider(
                self.slider_ax, "",
                XLim[0], XLim[1],
                valinit = self.initial_x,
                valfmt = f"%.2f {Xunit}"
            )
            self.x_slider.on_changed(self.update)

            # Initialize container for slider intersections with curves
            self.slider_intersections = {}

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
        for name, entry in self.Curve.items():
            # Skip vertical lines - they have their own text that should stay visible
            if entry.get('IsVerticalLine', False):
                continue

            x, y = entry['Curve'].get_xdata(), entry['Curve'].get_ydata()
            yv = np.interp(t, x, y)
            vb = entry['ValueBox']
            vb.set_position((t, yv))
            vb.set_text(f"{yv:.3f}" if self.XLim[0] <= t <= self.XLim[1] else "")

        # Handle slider intersections with curves
        if hasattr(self, 'slider_intersections') and self.slider_intersections:
            # Update or hide all intersection points for the slider
            for curve_name, intersection in self.slider_intersections.items():
                if t == self.XLim[0]:
                    # Hide at t=0
                    intersection['text'].set_visible(False)
                    if intersection.get('marker'):
                        intersection['marker'].set_visible(False)
                else:
                    # Show and update at other positions
                    curve_data = self.Curve.get(curve_name)
                    if curve_data and not curve_data.get('IsVerticalLine', False):
                        curve = curve_data['Curve']
                        x_data = curve.get_xdata()
                        y_data = curve.get_ydata()
                        if t >= min(x_data) and t <= max(x_data):
                            y_val = np.interp(t, x_data, y_data)
                            intersection['text'].set_position((t, y_val))
                            intersection['text'].set_text(f"{y_val:.3f}")
                            intersection['text'].set_visible(True)
                            if intersection.get('marker'):
                                intersection['marker'].set_visible(True)
                                intersection['marker'].set_data([t], [y_val])
                        else:
                            # Hide if outside curve range
                            intersection['text'].set_visible(False)
                            if intersection.get('marker'):
                                intersection['marker'].set_visible(False)

    def update(self, val):
        """
        Called when the slider moves. Updates slider line, shows/hides only the slider's value box and intersections at t=0,
        updates all curve value boxes, but leaves other vertical lines and their intersections untouched.
        """
        self.vline.set_xdata([val, val])

        # Show/hide the slider's own ValueBox
        if val == self.XLim[0]:
            # Hide only the slider's own value box at t=0
            self.slider_vb.set_visible(False)
        else:
            # Show and update text + position for slider's value box
            y_min = self.ax.get_ylim()[0]
            self.slider_vb.set_text(f"{val:.3f}")
            self.slider_vb.set_position((val, y_min))
            self.slider_vb.set_visible(True)

        # Always update the data boxes (including intersections) - this will handle hiding/showing
        # the intersection value boxes based on slider position
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
                               AutoScale=True, offset=offset)

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

        # Add slider intersection point if slider exists
        if hasattr(self, 'slider') and self.slider and hasattr(self, 'slider_intersections'):
            slider_pos = self.x_slider.val
            if slider_pos > self.XLim[0]:  # Only create if slider not at 0
                y_val = np.interp(slider_pos, Xdata, Ydata)
                ib = ax2.text(
                    slider_pos, y_val, f"{y_val:.3f}",
                    ha='left', va='bottom', color=line.get_color(),
                    bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1', alpha=1.0)
                )
                self.slider_intersections[Name] = {'text': ib}

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

        # Value box showing the x-position (this is the only label we'll keep)
        vb = self.ax.text(
            Xpos, YLims[0] - 0.05, f"{Xpos:.3f}",
            ha='center', va='top', color=line.get_color(),
            bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.1', alpha=1.0)
        )
        # Make the value box draggable
        vb.set_picker(True)

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
                    # Make the intersection value boxes draggable too
                    ib.set_picker(True)

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
        if self._dragging is None or event.x is None or event.y is None:
            return

        txt = self._dragging
        ax = txt.axes

        # convert mouse pixel → data coords in the text's own axes
        inv = ax.transData.inverted()
        xdata, ydata = inv.transform((event.x, event.y))

        txt.set_position((xdata, ydata))
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
                ax2.set_yticks([lo, 0.5*(lo+hi), hi])

    def AddFooter(self, info_dict, fontsize=8, color='gray',
                      pad=0.5, wrap_width=150):
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

            # force draw for accurate coords
            self.fig.canvas.draw()
            render = self.fig.canvas.get_renderer()

            # locate the xlabel in figure coords
            xlabel = self.ax.xaxis.label
            xb = xlabel.get_window_extent(render)
            xf = xb.transformed(self.fig.transFigure.inverted())

            # wrap each dict entry on its own line
            lines = [f"{k}: {v}" for k, v in info_dict.items()]
            max_chars = max(10, int(wrap_width * self.ax.get_position().width))
            text = "\n".join(textwrap.fill(l, max_chars) for l in lines)

            # position just below the label
            pad_fig = pad * xf.height
            x_fig = 0.5 * (xf.x0 + xf.x1)
            y_fig = xf.y0 - pad_fig

            footer = self.fig.text(
                x_fig, y_fig, text,
                transform=self.fig.transFigure,
                ha='center', va='top',
                fontsize=fontsize, color=color,
                family='monospace', linespacing=1.2,
            )
            footer.footer_tag = True
            self._footer_texts.append(footer)

            # Get footer dimensions
            self.fig.canvas.draw()
            fbbox = footer.get_window_extent(render)
            fbbox_fig = fbbox.transformed(self.fig.transFigure.inverted())

            # Calculate slider position
            slider_height = 0.03
            slider_padding = 0.01
            min_bottom = 0.01

            # If footer extends below minimum height, adjust layout
            if fbbox_fig.y0 < min_bottom:
                old_bottom = plt.rcParams['figure.subplot.bottom']
                required_space = min_bottom - fbbox_fig.y0
                new_bottom = old_bottom + required_space

                # Store original positions
                orig_positions = {name: info['ax'].get_position().bounds for name, info in self.Axis.items()}
                main_pos = self.ax.get_position()

                # Update main axis position
                scale_factor = (1.0 - new_bottom) / (1.0 - old_bottom)
                new_height = main_pos.height * scale_factor
                self.ax.set_position([main_pos.x0, new_bottom, main_pos.width, new_height])

                # Update all sub-axes proportionally
                for name, (x0, y0, width, height) in orig_positions.items():
                    rel_y = (y0 - main_pos.y0) / main_pos.height  # Relative position
                    rel_height = height / main_pos.height  # Relative height

                    new_y = new_bottom + (rel_y * new_height)
                    new_height_sub = rel_height * new_height

                    self.Axis[name]['ax'].set_position([x0, new_y, width, new_height_sub])

                # Redraw and reposition footer with new coordinates
                self.fig.canvas.draw()
                lbl_box = xlabel.get_window_extent(render)
                lbl_box_fig = lbl_box.transformed(self.fig.transFigure.inverted())
                y_fig = lbl_box_fig.y0 - pad_fig
                footer.set_position((x_fig, y_fig))
                fbbox = footer.get_window_extent(render)
                fbbox_fig = fbbox.transformed(self.fig.transFigure.inverted())

            # Calculate slider position
            slider_height = 0.03
            slider_padding = 0.01
            if self.slider:
                slider_y = fbbox_fig.y0 - slider_height - slider_padding
            else:
                slider_y = fbbox_fig.y0 - slider_padding

            # Position slider if present
            if self.slider:
                ext = self.ax.get_position()
                self.slider_ax.set_position([ext.x0, slider_y, ext.width, slider_height])
                self.fig.canvas.draw_idle()

    def _on_resize(self, event):
        """Handle resize events by properly positioning all elements with correct spacing"""
        # Disable conflicting layout systems
        self.fig.set_constrained_layout(False)
        self.fig.set_tight_layout(False)

        # Get figure dimensions
        fig_width, fig_height = self.fig.get_size_inches() * self.fig.dpi

        # Get minimum width needed for axis labels
        min_left_margin = 0.05
        min_right_margin = 0.05

        # Calculate minimum required margins by measuring axis labels
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()

        # Check if we have any axes to determine required margins
        left_offset = 0
        right_offset = 0

        for name, info in self.Axis.items():
            ax2 = info['ax']
            position = info['Position'].lower()
            offset = info.get('offset', 0.02)

            # Get label width to ensure it's visible
            if position == 'left':
                if ax2.yaxis.label.get_text():
                    bbox = ax2.yaxis.label.get_window_extent(renderer)
                    bbox_fig = bbox.transformed(self.fig.transFigure.inverted())
                    left_offset = max(left_offset, bbox_fig.width * (1 + offset))
            else:  # right position
                if ax2.yaxis.label.get_text():
                    bbox = ax2.yaxis.label.get_window_extent(renderer)
                    bbox_fig = bbox.transformed(self.fig.transFigure.inverted())
                    right_offset = max(right_offset, bbox_fig.width * (1 + offset))

        # Add minimum padding for left/right margins based on actual axis labels
        min_left_margin = max(min_left_margin, left_offset + 0.05)
        min_right_margin = max(min_right_margin, right_offset + 0.05)

        # Dynamic margins - adjust margins based on window size but ensure axis labels are visible
        # For larger windows, use smaller relative margins to maximize plot area
        # For smaller windows, ensure all axis labels remain visible
        h_scale_factor = min(1.0, max(0.5, fig_width / 1000))  # Scale based on width
        v_scale_factor = min(1.0, max(0.5, fig_height / 800))  # Scale based on height

        # Calculate horizontal margins - ensure they're large enough for axis labels
        left_margin = max(min_left_margin, 0.12 * h_scale_factor)
        right_margin = min(1.0 - min_left_margin, 1.0 - (0.12 * h_scale_factor))

        # Calculate vertical margins
        top_margin = min(0.98, 0.95 + (0.03 * (1 - v_scale_factor)))
        bottom_margin = max(0.05, 0.15 * v_scale_factor)

        # First apply dynamic margins
        plt.subplots_adjust(left=left_margin, right=right_margin,
                            bottom=bottom_margin, top=top_margin)

        # Force draw to get accurate positions
        self.fig.canvas.draw()
        main = self.ax.get_position()

        # 1. Position all sub-axes based on main axis
        for name, info in self.Axis.items():
            ax2 = info['ax']
            cell_h = main.height / self.Ngridy
            center = main.y0 + (info['GridPos'] + 0.5) * cell_h
            height = info['GridHeight'] * cell_h
            bottom = center - height/2
            ax2.set_position([main.x0, bottom, main.width, height])
            position = info['Position'].lower()
            offset = info.get('offset', 0.02)  # Get stored offset
            if position == 'left':
                ax2.spines['right'].set_visible(False)
                ax2.spines['left'].set_position(('axes', -offset))
            else:
                ax2.yaxis.tick_right()
                ax2.yaxis.set_label_position('right')
                ax2.spines['left'].set_visible(False)
                ax2.spines['right'].set_position(('axes', 1 + offset))
            ax2.set_xlim(self.ax.get_xlim())

        # 2. Apply autoscaling to axes
        for info in self.Axis.values():
            if info.get('AutoScale', False):
                ax2 = info['ax']
                lines = ax2.get_lines()
                if lines:
                    ys = np.hstack([ln.get_ydata() for ln in lines])
                    if ys.size:
                        lo, hi = np.floor(ys.min()), np.ceil(ys.max())
                        ax2.set_ylim(lo, hi)
                        ax2.set_yticks([lo, 0.5*(lo+hi), hi])

        # 3. Calculate space needed for all elements
        # Get space for title if present
        title_height = 0
        if hasattr(self, '_title_metadata'):
            renderer = self.fig.canvas.get_renderer()
            if self.ax.title and self.ax.title.get_text():
                title_box = self.ax.title.get_window_extent(renderer)
                title_box_fig = title_box.transformed(self.fig.transFigure.inverted())
                title_height = title_box_fig.height + 0.02  # Add padding

        # Calculate space needed for slider
        slider_height = 0
        slider_padding = 0
        if self.slider:
            slider_height = 0.03  # Fixed height for slider
            slider_padding = 0.08  # Increased padding above slider for footer

        # 4. Calculate space needed for footer
        footer_height = 0
        if hasattr(self, '_footer_metadata'):
            # Clear old footer
            for txt in getattr(self, '_footer_texts', []):
                txt.remove()
            self._footer_texts = []

            # Temporarily position footer to measure it
            info_dict, fontsize, color, pad, wrap_width = self._footer_metadata
            lines = [f"{k}: {v}" for k, v in info_dict.items()]
            max_chars = max(10, int(wrap_width * main.width))
            text = "\n".join(textwrap.fill(l, max_chars) for l in lines)

            # Create temp footer for measurement
            temp_footer = self.fig.text(
                0.5, 0.1, text,  # Temporary position
                transform=self.fig.transFigure,
                ha='center', va='top',
                fontsize=fontsize, color=color,
                family='monospace', linespacing=1.2,
            )

            # Measure footer height
            self.fig.canvas.draw()
            renderer = self.fig.canvas.get_renderer()
            fb = temp_footer.get_window_extent(renderer)
            ff = fb.transformed(self.fig.transFigure.inverted())
            footer_height = ff.height + 0.03  # Add padding
            temp_footer.remove()

        # 5. Adjust figure layout to accommodate all elements
        # Calculate minimum required bottom margin
        required_bottom = 0.02  # Minimum padding at bottom

        # Calculate exact space needed for slider if present
        if self.slider:
            required_bottom += slider_height  # Slider height

        # Calculate exact space needed for footer if present
        if footer_height > 0:
            # If slider is present, place footer right above it
            if self.slider:
                # We need to ensure the bottom margin includes BOTH the slider and footer heights
                required_bottom = max(required_bottom, 0.02 + slider_height + footer_height + 0.01)
            else:
                # Just space for footer
                required_bottom = max(required_bottom, 0.02 + footer_height)

        # Calculate required top margin with minimal padding
        required_top = 1.0 - title_height - 0.01  # Minimal padding at top

        # Ensure our margins don't go below the minimum required
        bottom_margin = max(bottom_margin, required_bottom)
        top_margin = min(top_margin, required_top)

        # Apply final margins to maximize plot space while keeping all elements visible
        plt.subplots_adjust(bottom=bottom_margin, top=top_margin)

        # Redraw to update positions
        self.fig.canvas.draw()
        main = self.ax.get_position()

        # 6. Position the footer if it exists
        if hasattr(self, '_footer_metadata'):
            self.fig.canvas.draw()
            renderer = self.fig.canvas.get_renderer()

            # Get the position of the x-axis label
            xlabel = self.ax.xaxis.label
            xb = xlabel.get_window_extent(renderer)
            xf = xb.transformed(self.fig.transFigure.inverted())

            info_dict, fontsize, color, pad, wrap_width = self._footer_metadata
            lines = [f"{k}: {v}" for k, v in info_dict.items()]
            max_chars = max(10, int(wrap_width * main.width))
            text = "\n".join(textwrap.fill(l, max_chars) for l in lines)

            # Position footer directly below the x-axis label with proper spacing
            x_fig = 0.5 * (xf.x0 + xf.x1)  # Center horizontally based on label

            # Calculate y position based on the bottom of the xlabel with padding
            pad_fig = pad * xf.height  # Scale padding based on label height
            y_fig = xf.y0 - pad_fig  # Position just below the label

            # Add footer with top alignment to grow downward from this position
            footer = self.fig.text(
                x_fig, y_fig, text,
                transform=self.fig.transFigure,
                ha='center', va='top',  # Top alignment so text grows downward
                fontsize=fontsize, color=color,
                family='monospace', linespacing=1.2,
            )
            self._footer_texts.append(footer)

            # Get footer dimensions to position slider if needed
            self.fig.canvas.draw()
            fbbox = footer.get_window_extent(renderer)
            fbbox_fig = fbbox.transformed(self.fig.transFigure.inverted())

            # Ensure footer is within window bounds
            # If it would go out of bounds, adjust the plot accordingly
            if fbbox_fig.y0 < 0.02 + (slider_height if self.slider else 0):
                # Calculate how much the plot needs to shrink
                adjustment = 0.02 + (slider_height if self.slider else 0) - fbbox_fig.y0

                # Adjust main plot position to make room
                current_pos = self.ax.get_position()
                new_bottom = current_pos.y0 + adjustment
                new_height = current_pos.height - adjustment

                # Apply the adjustment
                self.ax.set_position([current_pos.x0, new_bottom, current_pos.width, new_height])

                # Redraw and reposition the x-axis label and footer
                self.fig.canvas.draw()
                xb = xlabel.get_window_extent(renderer)
                xf = xb.transformed(self.fig.transFigure.inverted())

                # Recalculate footer position
                y_fig = xf.y0 - pad_fig
                footer.set_position((x_fig, y_fig))

                # Update footer dimensions
                self.fig.canvas.draw()
                fb = footer.get_window_extent(renderer)

            # Reposition all sub-axes after main axis adjustment
            self._update_subaxes_positions()

        # 7. Position slider at the bottom with minimal space
        if self.slider:
            slider_y = 0.02  # Fixed position at bottom with minimal padding
            self.slider_ax.set_position([main.x0, slider_y, main.width, slider_height])

        # 8. Update vertical lines and intersections
        for name, data in self.Curve.items():
            if data.get('IsVerticalLine', False):
                xpos = data.get('XPosition')
                line = data['Curve']
                line.set_xdata([xpos, xpos])

        # Final redraw
        self.fig.canvas.draw_idle()

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
        # plt.subplots_adjust(top=required_top)

        # Store title for resize events
        self._title_metadata = (title, fontsize, pad, kwargs)

        # Redraw the figure
        self.fig.canvas.draw_idle()

    def _update_subaxes_positions(self):
        """Update the positions of all sub-axes based on the main axis position"""
        main = self.ax.get_position()

        for name, info in self.Axis.items():
            ax2 = info['ax']
            cell_h = main.height / self.Ngridy
            center = main.y0 + (info['GridPos'] + 0.5) * cell_h
            height = info['GridHeight'] * cell_h
            bottom = center - height/2
            ax2.set_position([main.x0, bottom, main.width, height])

            # Get the position and offset from axis configuration
            position = info['Position'].lower()
            offset = info.get('offset', 0.02)  # Default to 0.02 if not stored

            # Apply the correct spine position with stored offset
            if position == 'left':
                ax2.spines['right'].set_visible(False)
                ax2.spines['left'].set_position(('axes', -offset))
            else:
                ax2.yaxis.tick_right()
                ax2.yaxis.set_label_position('right')
                ax2.spines['left'].set_visible(False)
                ax2.spines['right'].set_position(('axes', 1 + offset))

            ax2.set_xlim(self.ax.get_xlim())


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
    Fdr.AddAxis(Name='Axis 3',GridPos=5,Unit='deg',offset=0.1)
    Fdr.AddCurve('C3', 'Axis 3', Xdata, Ydata, color='m')

    Fdr.RemoveCurve('C1')
    Fdr.RemoveAxis('Axis 2')
    Fdr.SetAxisLimits('Axis 3', -6, 6)
    Fdr.SetAxisHeight('Axis 3', 3)
    Fdr.AddVerticalLine(Name='0.5', Xpos=0.5, linestyle='--')
    Fdr.AddVerticalLine(Name='0.7', Xpos=0.7, linestyle=':', color='green')
    Fdr.AddVerticalLine(Name='0.3', Xpos=0.3, linestyle='-.', marker='s', color='blue')
    Fdr.SetLabelPosition('C2', 0.5, 0.8)
