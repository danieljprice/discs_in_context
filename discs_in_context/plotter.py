"""
Main plotting class for star-forming regions with extinction maps and disc positions.
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import time

import astropy.units as units
from astropy.coordinates import SkyCoord

from dustmaps.sfd import SFDQuery
from dustmaps.planck import PlanckQuery
from dustmaps.bayestar import BayestarQuery

from .regions import REGIONS

class plotcloud:
    """
    A class for plotting extinction maps of star-forming regions with
    protostars and protoplanetary discs overlaid.

    Can be initialized with:
    - Preset region labels (e.g., 'taurus', 'lupus', 'orion')
    - RA/Dec coordinate ranges
    - Galactic l/b coordinate ranges
    """

    def __init__(self, region=None, ra_range=None, dec_range=None,
                 galactic_l=None, galactic_b=None, coord_system='galactic',
                 num_points=2048, object=None, csvfile=None, image_size=None,
                 image_size_unit='arcsec'):
        """
        Initialize the plotter with a region specification.

        Parameters
        ----------
        region : str, optional
            Preset region name (e.g., 'taurus', 'lupus', 'orion').
            If provided, other parameters are ignored.
        ra_range : tuple, optional
            (ra_min, ra_max) in hours. Requires dec_range.
        dec_range : tuple, optional
            (dec_min, dec_max) in degrees. Requires ra_range.
        galactic_l : dict, optional
            {'l_min': min_l, 'l_max': max_l} in degrees.
            Requires galactic_b.
        galactic_b : dict, optional
            {'b_min': min_b, 'b_max': max_b} in degrees.
            Requires galactic_l.
        coord_system : str, default 'galactic'
            Coordinate system to use: 'galactic' or 'icrs'.
            Only used when region is None.
        num_points : int, default 2048
            Number of grid points for the extinction map.
        object : str, optional
            Object name to look up in the discs CSV file (e.g., 'HD 142527').
            If provided, initializes the map centred on this object using
            coordinates from the discs catalogue.
        csvfile : str, optional
            Path to discs CSV file. If None, uses the copy bundled in the
            package data directory.
        image_size : float, optional
            Size of the image in arcsec or degrees around the object.
        image_size_unit : str, default 'arcsec'
            Unit for image_size: 'arcsec' or 'degrees'.
        """
        self.num_points = num_points
        self.coord_system = coord_system
        self.prev_positions = []
        self.overlap_range = 0.15
        self.overlap_range_ra = 0.1
        self.region = region  # Store region name for title

        # Initialize coordinates based on input
        if region is not None:
            self._init_from_region(region)
        elif object is not None:
            if image_size is None:
                raise ValueError(
                    "If 'object' is provided, 'image_size' is also required"
                )
            self._init_from_object(object, csvfile, image_size, image_size_unit)
        elif ra_range is not None and dec_range is not None:
            self._init_from_ra_dec(ra_range, dec_range)
        elif galactic_l is not None and galactic_b is not None:
            self._init_from_galactic(galactic_l, galactic_b)
        else:
            raise ValueError(
                "Must provide either 'region', 'object' (with csvfile and image_size), "
                "('ra_range', 'dec_range'), or ('galactic_l', 'galactic_b')"
            )

        # Generate coordinate grid
        self._generate_coord_grid()

    def _init_from_region(self, region):
        """Initialize from a preset region name."""
        region_lower = region.lower()
        if region_lower not in REGIONS:
            available = ', '.join(REGIONS.keys())
            raise ValueError(
                f"Unknown region '{region}'. Available regions: {available}"
            )

        region_config = REGIONS[region_lower]
        if self.coord_system == 'galactic' and 'galactic' in region_config:
            config = region_config['galactic']
            self._init_from_galactic(
                {'l_min': config['l_min'], 'l_max': config['l_max']},
                {'b_min': config['b_min'], 'b_max': config['b_max']}
            )
        elif self.coord_system == 'icrs' and 'icrs' in region_config:
            config = region_config['icrs']
            self._init_from_ra_dec(
                (config['ra_min'], config['ra_max']),
                (config['dec_min'], config['dec_max'])
            )
        else:
            # Default to galactic if available, otherwise icrs
            if 'galactic' in region_config:
                config = region_config['galactic']
                self.coord_system = 'galactic'
                self._init_from_galactic(
                    {'l_min': config['l_min'], 'l_max': config['l_max']},
                    {'b_min': config['b_min'], 'b_max': config['b_max']}
                )
            else:
                config = region_config['icrs']
                self.coord_system = 'icrs'
                self._init_from_ra_dec(
                    (config['ra_min'], config['ra_max']),
                    (config['dec_min'], config['dec_max'])
                )

    def _init_from_ra_dec(self, ra_range, dec_range):
        """Initialize from RA/Dec ranges."""
        self.coord_system = 'icrs'
        self.ra_min, self.ra_max = ra_range
        self.dec_min, self.dec_max = dec_range
        self.ra_span = self.ra_max - self.ra_min
        self.dec_span = self.dec_max - self.dec_min
        self.max_span = max(self.ra_span, self.dec_span / 15.0)

    def _init_from_galactic(self, galactic_l, galactic_b):
        """Initialize from Galactic l/b ranges."""
        self.coord_system = 'galactic'
        self.l_min = galactic_l['l_min']
        self.l_max = galactic_l['l_max']
        self.b_min = galactic_b['b_min']
        self.b_max = galactic_b['b_max']
        self.l_span = self.l_max - self.l_min
        self.b_span = self.b_max - self.b_min
        self.max_span = max(self.l_span, self.b_span)

    def _init_from_object(self, object, csvfile, image_size, image_size_unit='arcsec'):
        """
        Initialize from an object name looked up in a CSV file.

        Parameters
        ----------
        object : str
            Object name to search for in the CSV file (e.g., 'HD 142527').
        csvfile : str or None
            Path to CSV file containing disc data. If None, uses the
            bundled discs catalogue in the package data directory.
        image_size : float
            Size of the image in arcsec or degrees.
        image_size_unit : str, default 'arcsec'
            Unit for image_size: 'arcsec' or 'degrees'.
        """
        import os
        import pandas as pd
        from astropy.coordinates import SkyCoord
        import astropy.units as u

        # If no CSV path is provided, use the bundled discs catalogue
        if csvfile is None:
            from . import data_utils

            csvfile = data_utils.get_discs_catalog_path()

        # Expand user path and read CSV file
        csvfile = os.path.expanduser(str(csvfile))
        try:
            df = pd.read_csv(csvfile, usecols=['target_id', 'l', 'b'])
        except KeyError as exc:
            raise ValueError(
                "Discs CSV file must contain columns: 'target_id', 'l', 'b'"
            ) from exc

        # Find the object
        matches = df[df['target_id'] == object]
        if len(matches) == 0:
            raise ValueError(f"Object '{object}' not found in CSV file '{csvfile}'")
        if len(matches) > 1:
            raise ValueError(f"Multiple entries found for object '{object}' in CSV file")

        # Get coordinates
        l = matches.iloc[0]['l']
        b = matches.iloc[0]['b']

        # Convert image size to degrees
        if image_size_unit == 'arcsec':
            size_deg = image_size / 3600.0
        elif image_size_unit == 'degrees':
            size_deg = image_size
        else:
            raise ValueError(f"image_size_unit must be 'arcsec' or 'degrees', got '{image_size_unit}'")

        # Create coordinate range centered on the source
        half_size = size_deg / 2.0

        # Convert to the desired coordinate system
        if self.coord_system == 'galactic':
            # Use galactic coordinates directly
            self.l_min = l - half_size
            self.l_max = l + half_size
            self.b_min = b - half_size
            self.b_max = b + half_size
            self.l_span = self.l_max - self.l_min
            self.b_span = self.b_max - self.b_min
            self.max_span = max(self.l_span, self.b_span)
        else:
            # Convert to ICRS (RA/Dec)
            coord = SkyCoord(l=l*u.deg, b=b*u.deg, frame='galactic')
            ra_center = coord.icrs.ra.hourangle
            dec_center = coord.icrs.dec.deg

            # For RA, account for the fact that 1 degree of RA = 1/15 hours
            # For small fields, we can approximate the size in RA hours
            ra_half_size = half_size / 15.0
            dec_half_size = half_size

            self.ra_min = ra_center - ra_half_size
            self.ra_max = ra_center + ra_half_size
            self.dec_min = dec_center - dec_half_size
            self.dec_max = dec_center + dec_half_size
            self.ra_span = self.ra_max - self.ra_min
            self.dec_span = self.dec_max - self.dec_min
            self.max_span = max(self.ra_span, self.dec_span / 15.0)

        # Store object name for title
        self.object = object

    def _generate_coord_grid(self):
        """Generate the coordinate grid for extinction map calculation."""
        if self.coord_system == 'icrs':
            ra_step = self.max_span / self.num_points
            dec_step = self.max_span / (self.num_points * 15)
            ra = np.linspace(self.ra_min,self.ra_max,self.num_points)
            dec = np.linspace(self.dec_min,self.dec_max,self.num_points)
            ra, dec = np.meshgrid(ra, dec)
            self.coords = SkyCoord(
                ra * units.hourangle,
                dec * units.deg,
                frame='icrs',
                obstime="J2000"
            )
        else:  # galactic
            l = np.linspace(self.l_min,self.l_max,self.num_points)
            b = np.linspace(self.b_min,self.b_max,self.num_points)
            l, b = np.meshgrid(l, b)
            self.coords = SkyCoord(l * units.deg,b * units.deg,frame='galactic')

    def _apply_scaling(self, ax, scale_factor, default_fontsize):
        """
        Apply scaling to all plot elements for interactive display.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to scale.
        scale_factor : float
            Scaling factor (typically 0.2 for interactive plots).
        default_fontsize : float
            Default font size to scale from.

        Returns
        -------
        scaled_fontsize : float
            Scaled font size.
        labelpad : float
            Scaled label padding.
        """
        scaled_fontsize = default_fontsize * scale_factor
        tick_pad = 2 * scale_factor
        title_pad = 4 * scale_factor
        labelpad = 3 * scale_factor

        # Scale fonts
        ax.tick_params(labelsize=scaled_fontsize, width=scale_factor,
                      length=4 * scale_factor, pad=tick_pad)
        ax.xaxis.label.set_fontsize(scaled_fontsize)
        ax.yaxis.label.set_fontsize(scaled_fontsize)
        ax.title.set_fontsize(scaled_fontsize)

        # Scale line widths
        for spine in ax.spines.values():
            spine.set_linewidth(0.5 * scale_factor)

        # Scale grid line widths (if grid exists)
        for line in ax.xaxis.get_gridlines():
            line.set_linewidth(0.3 * scale_factor)
        for line in ax.yaxis.get_gridlines():
            line.set_linewidth(0.3 * scale_factor)

        # Store for later use
        self._title_pad = title_pad
        self._scaled_fontsize = scaled_fontsize
        self._interactive_scale = scale_factor

        return scaled_fontsize, labelpad

    def _restore_scaling(self, ax, default_fontsize):
        """
        Restore original (non-scaled) sizes for saving.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to restore.
        default_fontsize : float
            Default font size to restore to.
        """
        ax.tick_params(labelsize=default_fontsize, width=1.0, length=4.0, pad=2.0)
        ax.xaxis.label.set_fontsize(default_fontsize)
        ax.yaxis.label.set_fontsize(default_fontsize)
        ax.title.set_fontsize(default_fontsize)

        # Restore line widths
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)
        for line in ax.xaxis.get_gridlines():
            line.set_linewidth(0.3)
        for line in ax.yaxis.get_gridlines():
            line.set_linewidth(0.3)

        # Restore padding
        ax.set_xlabel(ax.get_xlabel(), labelpad=3)
        ax.set_ylabel(ax.get_ylabel(), labelpad=3)
        ax.set_title(ax.get_title(), pad=4.0, fontsize=default_fontsize)

    def _place_non_overlapping_label(self, ax, x, y, label, font_scale,
                                     color='cyan', ha='left',
                                     offset_factor=0.05):
        """
        Place a text label near (x, y) while trying to avoid overlap with
        previously drawn labels.

        This works for both RA/Dec and Galactic coordinates since it only
        uses the plotted x/y positions and the shared overlap ranges.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to draw on.
        x, y : float
            Data coordinates of the associated point.
        label : str
            Text label to draw.
        font_scale : float
            Scale factor applied to the base font size and overlap ranges.
        color : str, optional
            Text colour.
        ha : str, optional
            Horizontal alignment ('left' or 'right').
        offset_factor : float, optional
            Horizontal offset (in data units) as a multiple of font_scale.
        """
        # Horizontal offset from the point
        offset_x = offset_factor * font_scale
        label_x = x - offset_x if ha == 'right' else x + offset_x
        label_y = y

        # Scale overlap ranges with font size
        scaled_overlap_x = self.overlap_range_ra * font_scale
        scaled_overlap_y = self.overlap_range * font_scale

        # Shift the label upwards until it no longer overlaps previous ones
        while any(
            abs(label_x - x_prev) <= scaled_overlap_x and
            abs(label_y - y_prev) <= scaled_overlap_y
            for x_prev, y_prev in self.prev_positions
        ):
            label_y += scaled_overlap_y

        ax.text(
            label_x,
            label_y,
            label,
            color=color,
            va='top',
            ha=ha,
            size=10 * font_scale,
            clip_on=True,
            rotation=0.0,
        )
        self.prev_positions.append((label_x, label_y))

    def _setup_interactive_labels(self, ax, scatter_points, labels_data, scale_factor=1.0):
        """
        Set up interactive label display on hover/click for scatter points.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes containing the scatter points.
        scatter_points : list
            List of scatter plot objects.
        labels_data : list
            List of (label, x, y) tuples for each point.
        scale_factor : float, default 1.0
            Scale factor for annotation size (for interactive plots).
        """
        # Create annotation that will show labels
        # Scale annotation size for interactive plots
        annot = ax.annotate('', xy=(0, 0), xytext=(20 * scale_factor, 20 * scale_factor),
                           textcoords="offset points",
                           bbox=dict(boxstyle="round", fc="white", alpha=0.8,linewidth=0.5 * scale_factor),
                           arrowprops=dict(arrowstyle="->", linewidth=0.5 * scale_factor),
                           fontsize=10 * scale_factor)
        annot.set_visible(False)

        # Create a mapping from scatter plot to label indices
        scatter_to_labels = {}
        label_idx = 0
        for scatter in scatter_points:
            offsets = scatter.get_offsets()
            if len(offsets.shape) == 2:
                num_points = len(offsets)
            else:
                num_points = 1
            scatter_to_labels[scatter] = list(range(label_idx, label_idx + num_points))
            label_idx += num_points

        def update_annot(scatter, idx):
            """Update annotation for a given scatter point and index."""
            label_list_idx = scatter_to_labels[scatter][idx]
            if label_list_idx < len(labels_data):
                label, x, y = labels_data[label_list_idx]
                annot.xy = (x, y)
                annot.set_text(label)
                annot.set_visible(True)
                ax.figure.canvas.draw_idle()

        def hover(event):
            """Handle mouse hover events."""
            if event.inaxes == ax:
                vis = annot.get_visible()
                found = False
                for scatter in scatter_points:
                    cont, ind = scatter.contains(event)
                    if cont:
                        idx = ind["ind"][0]
                        update_annot(scatter, idx)
                        found = True
                        break

                if not found and vis:
                    annot.set_visible(False)
                    ax.figure.canvas.draw_idle()

        def click(event):
            """Handle mouse click events."""
            if event.inaxes == ax:
                for scatter in scatter_points:
                    cont, ind = scatter.contains(event)
                    if cont:
                        idx = ind["ind"][0]
                        update_annot(scatter, idx)
                        break

        # Connect event handlers
        ax.figure.canvas.mpl_connect("motion_notify_event", hover)
        ax.figure.canvas.mpl_connect("button_press_event", click)

    def plot_all_discs(self, ax, csvfile=None, interactive=False, font_scale=1.0, marker_scale=1.0):
        """
        Plot all protoplanetary discs from a CSV file.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to plot on.
        csvfile : str, optional
            Path to CSV file with disc data. If None, uses default path.
        interactive : bool, default False
            If True, labels are shown on hover/click instead of always visible.
        """
        import os
        import pandas as pd

        # If no CSV path is provided, use the bundled discs catalogue
        if csvfile is None:
            from . import data_utils

            csvfile = data_utils.get_discs_catalog_path()

        csvfile = os.path.expanduser(str(csvfile))
        df = pd.read_csv(csvfile)

        # Ensure required columns exist
        required_cols = {'target_id', 'l', 'b'}
        missing = required_cols.difference(df.columns)
        if missing:
            missing_str = ', '.join(sorted(missing))
            raise ValueError(f"Discs CSV file is missing required columns: {missing_str}")

        # Optional distance column (specific to discs catalogue)
        disc_distance_col = None
        for cand in ('Distance', 'distance_pc', 'dist_pc', 'distance', 'dist'):
            if cand in df.columns:
                disc_distance_col = cand
                break
        nd = df.shape[0]
        print(f"got {nd} discs")

        # Get plot limits
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        x_min, x_max = min(xlim), max(xlim)
        y_min, y_max = min(ylim), max(ylim)

        # Store scatter points for interactive mode
        scatter_points = []
        labels_data = []
        disc_label_set = set()

        for i in range(nd):
            row = df.loc[i, :]
            label = row['target_id']
            l = row['l']
            b = row['b']

            # Build display label for interactive mode (optionally including distance)
            if disc_distance_col is not None:
                dist_val = row[disc_distance_col]
                try:
                    if np.isfinite(dist_val):
                        display_label = f"{label} ({dist_val:.0f} pc)"
                    else:
                        display_label = label
                except TypeError:
                    display_label = label
            else:
                display_label = label

            # Transform to the plotting coordinate system
            if self.coord_system == 'icrs':
                # Convert galactic (l, b) to ICRS (RA, Dec)
                co = SkyCoord(l * units.deg, b * units.deg, frame='galactic')
                x = co.icrs.ra.degree
                y = co.icrs.dec.degree
                # Account for potentially reversed x-limits in RA (degrees)
                in_x_range = (x_min <= x <= x_max) if x_min < x_max else (x_max <= x <= x_min)
            else:
                # Plot directly in galactic coordinates
                x = l
                y = b
                in_x_range = x_min <= x <= x_max

            # Only plot if within plot limits
            if in_x_range and y_min <= y <= y_max:
                # Discs: cyan markers, drawn above PMS using higher zorder.
                scatter = ax.scatter(x, y, marker='*', s=5 * marker_scale,
                                     color='cyan', zorder=3)
                scatter_points.append(scatter)
                # Use display_label for interactive annotations
                labels_data.append((display_label, x, y))
                disc_label_set.add(label)
                if not interactive:
                    # Use shared helper to avoid overlapping labels
                    self._place_non_overlapping_label(
                        ax,
                        x,
                        y,
                        label,
                        font_scale,
                        color='cyan',
                        ha='right',
                        offset_factor=0.02,
                    )

        # Remember which labels correspond to discs so that PMS labels
        # for the same objects can be suppressed to avoid duplicates.
        self._disc_label_set = disc_label_set

        # For interactive labels, defer setting up the event handlers until
        # PMS sources have also been added, so that there is a single shared
        # annotation object for both discs and PMS.
        if interactive and scatter_points:
            pending = getattr(self, '_pending_interactive', {'scatter': [], 'labels': []})
            pending['scatter'].extend(scatter_points)
            pending['labels'].extend(labels_data)
            self._pending_interactive = pending

        # Set up interactive label display
        if interactive and scatter_points:
            scale = getattr(self, '_interactive_scale', 1.0)
            self._setup_interactive_labels(ax, scatter_points, labels_data, scale_factor=scale)

    def plot_kenyon08_pms(self, ax, csvfile=None, only_label_famous=False,
                          interactive=False, font_scale=1.0, marker_scale=1.0):
        """
        Plot PMS sources from Kenyon 2008 catalog.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to plot on.
        csvfile : str, optional
            Path to CSV file with PMS data. If None, uses default path.
        only_label_famous : bool, default False
            If True, only label "famous" objects (non-catalogue names)
            based on a built-in heuristic. If False, labels are plotted
            for all sources within the plot region.
        interactive : bool, default False
            If True, labels are shown on hover/click instead of always visible.
        """
        import os
        import pandas as pd

        # Load PMS catalogue:
        # - If csvfile is None, look up the copy bundled inside the
        #   discs_in_context/data directory.
        # - If csvfile is provided, fall back to reading from the user-specified path.
        if csvfile is None:
            from . import data_utils

            csvfile = data_utils.get_tau_sources_path()

        csvfile = os.path.expanduser(str(csvfile))
        df = pd.read_csv(csvfile)

        # Ensure required columns exist
        required_cols = {'PMS', 'RA', 'Dec'}
        missing = required_cols.difference(df.columns)
        if missing:
            missing_str = ', '.join(sorted(missing))
            raise ValueError(f"PMS CSV file is missing required columns: {missing_str}")

        nd = df.shape[0]
        print(f"got {nd} PMS sources")

        def _is_famous_label(label: str) -> bool:
            """Return True if label is considered 'famous' (not a catalogue ID)."""
            return not (label.startswith("J0") or
                        label.startswith("I0") or
                        label.startswith("IC") or
                        label.startswith("J1") or
                        label.startswith("JH") or
                        label.startswith("ITG") or
                        label.startswith("Haro") or
                        label.startswith("XEST") or
                        label.startswith("CFHT") or
                        label.startswith("V410") or
                        label.startswith("MHO") or
                        label.startswith("LR1") or
                        label.startswith("Hubble") or
                        label.startswith("Anon") or
                        label.startswith("KPNO"))

        # Get plot limits
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        x_min, x_max = min(xlim), max(xlim)
        y_min, y_max = min(ylim), max(ylim)

        # Store scatter points and labels for interactive mode (PMS sources
        # use white markers; interactive labels show the PMS name only).
        scatter_points = []
        labels_data = []

        # If discs have been plotted, avoid duplicating labels for objects
        # that are present in both the discs catalogue and the PMS catalogue.
        disc_label_set = getattr(self, '_disc_label_set', set())

        for i in range(nd):
            row = df.loc[i, :]
            label = row['PMS']
            ra_str = row['RA']
            dec_str = row['Dec']
            co1 = SkyCoord(
                [ra_str + dec_str],
                frame='icrs',
                unit=(units.hourangle, units.deg),
                obstime="J2000"
            )

            if self.coord_system == 'icrs':
                ra1, dec1 = (co1.ra.degree[0], co1.dec.degree[0])
                # RA is already in degrees, x-axis is also in degrees
                # Account for potentially reversed xlim
                ra1_deg = ra1
                in_x_range = (x_min <= ra1_deg <= x_max) if x_min < x_max else (x_max <= ra1_deg <= x_min)
                # Only plot if within plot limits
                if in_x_range and y_min <= dec1 <= y_max:
                    plot_label = _is_famous_label(label) if only_label_famous else True
                    if label in disc_label_set:
                        plot_label = False

                    scatter = ax.scatter(ra1, dec1, marker='*', s=10 * marker_scale,
                                         color='white', zorder=2)
                    if interactive and plot_label:
                        scatter_points.append(scatter)
                        labels_data.append((label, ra1, dec1))
                    if not interactive and plot_label:
                        # Use shared helper to avoid overlapping static labels
                        self._place_non_overlapping_label(ax,ra1,dec1,label,font_scale,color='cyan',
                                                          ha='left',offset_factor=0.05)
            else:  # galactic
                galactic_coords = co1.galactic
                l, b = (galactic_coords.l.degree[0], galactic_coords.b.degree[0])
                # Only plot if within plot limits
                if x_min <= l <= x_max and y_min <= b <= y_max:
                    plot_label = _is_famous_label(label) if only_label_famous else True
                    if label in disc_label_set:
                        plot_label = False

                    scatter = ax.scatter(l, b, marker='*', s=10 * marker_scale,
                                         color='white', zorder=2)
                    if interactive and plot_label:
                        scatter_points.append(scatter)
                        labels_data.append((label, l, b))
                    if not interactive and plot_label:
                        # For galactic axes, treat l like RA and b like Dec
                        self._place_non_overlapping_label(ax,l,b,label,font_scale,color='cyan',
                                                          ha='left',offset_factor=0.05)

        # Set up interactive label display for discs+PMS together. Discs may have
        # stored "pending" interactive points; combine them here so there is only
        # one shared annotation object and event handler.
        if interactive:
            pending = getattr(self, '_pending_interactive', {'scatter': [], 'labels': []})
            pending['scatter'].extend(scatter_points)
            pending['labels'].extend(labels_data)
            if pending['scatter']:
                scale = getattr(self, '_interactive_scale', 1.0)
                self._setup_interactive_labels(
                    ax,
                    pending['scatter'],
                    pending['labels'],
                    scale_factor=scale,
                )
            # Clear pending so subsequent plots start fresh
            if hasattr(self, '_pending_interactive'):
                del self._pending_interactive

    def plot(self, dustmap='planck', figsize=(15, 12), dpi=300,
             vmin=0.0, vmax=2.0, cmap=None, plot_discs=False,
             plot_pms=False, pms_csvfile=None, discs_csvfile=None,
             only_label_famous=False, bayestar_distance=None,
             save_path=None, interactive=False):
        """
        Create the extinction map plot.

        Parameters
        ----------
        dustmap : str, default 'planck'
            Dust map to use: 'planck', 'sfd', or 'bayestar'.
        figsize : tuple, default (18, 10)
            Figure size in inches.
        dpi : int, default 300
            Figure resolution.
        vmin : float, default 0.0
            Minimum value for colormap.
        vmax : float, default 2.0
            Maximum value for colormap.
        cmap : str, optional
            Colormap name. If None, uses 'binary'.
        plot_discs : bool, default False
            Whether to plot all discs.
        plot_pms : bool, default True
            Whether to plot PMS sources.
        pms_csvfile : str, optional
            Path to PMS CSV file.
        discs_csvfile : str, optional
            Path to discs CSV file.
        only_label_famous : bool, default False
            If True, only label "famous" PMS objects (non-catalogue names),
            using the same heuristic as in plot_kenyon08_pms.
        bayestar_distance : float or astropy.units.Quantity, optional
            Distance at which to evaluate the Bayestar dust map when
            ``dustmap='bayestar'``. If given as a float, it is interpreted
            as kiloparsecs. If None, the map is integrated along the entire
            line of sight using the maximum E(B−V) in each direction.
        save_path : str, optional
            Path to save figure. If None and interactive=False, doesn't save.
        interactive : bool, default False
            If True, display the figure interactively. If False and
            save_path is provided, saves the figure.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The created figure.
        ax : matplotlib.axes.Axes
            The axes object.
        """
        # Get extinction map
        if dustmap == 'planck':
            planck = PlanckQuery()
            av = 3.1 * planck(self.coords)
            dustmap_name = 'Planck'
        elif dustmap == 'sfd':
            sfd = SFDQuery()
            av = 2.742 * sfd(self.coords)
            dustmap_name = 'SFD'
        elif dustmap == 'bayestar':
            bayestar = BayestarQuery(max_samples=1)

            # If a distance is provided, evaluate Bayestar at that distance.
            # Otherwise, integrate along the line of sight by taking the
            # maximum E(B−V) value along the distance axis.
            if bayestar_distance is not None:
                from astropy import units as u

                dist = bayestar_distance
                if not hasattr(dist, 'unit'):
                    dist = dist * u.kpc

                if self.coord_system == 'icrs':
                    coords_dist = SkyCoord(
                        ra=self.coords.ra,
                        dec=self.coords.dec,
                        distance=dist,
                        frame='icrs'
                    )
                else:
                    coords_dist = SkyCoord(
                        l=self.coords.l,
                        b=self.coords.b,
                        distance=dist,
                        frame='galactic'
                    )
                ebv = bayestar(coords_dist)
            else:
                ebv = bayestar(self.coords)
                # Bayestar returns E(B−V) vs distance: collapse distance axis.
                if ebv.ndim == 3:
                    ebv = np.nanmax(ebv, axis=-1)

            av = 2.742 * ebv
            dustmap_name = 'Bayestar'
        else:
            raise ValueError(f"Unknown dustmap: {dustmap}")

        # Format title: include region name if using a preset region
        if self.region is not None and self.region != 'allsky':
            region_name = self.region.capitalize()
            title = f'{dustmap_name} extinction map of {region_name}'
        else:
            title = f'{dustmap_name} extinction map'

        # For interactive plots, use smaller figure size for screen display
        # but maintain the same proportions
        if interactive:
            # Scale down to fit screen better (e.g., 60% of original size)
            scale_factor = 0.2
            display_figsize = (figsize[0] * scale_factor, figsize[1] * scale_factor)
            display_dpi = dpi
        else:
            scale_factor = 1.0
            display_figsize = figsize
            display_dpi = dpi

        # Create figure
        fig = plt.figure(figsize=display_figsize, dpi=display_dpi)
        ax = fig.add_subplot(1, 1, 1)

        # Scale font sizes, line widths, and other elements for interactive display
        import matplotlib as mpl
        default_fontsize = mpl.rcParams['font.size']
        if interactive:
            scaled_fontsize, labelpad = self._apply_scaling(ax, scale_factor, default_fontsize)
        else:
            scaled_fontsize = default_fontsize
            labelpad = 3

        # Set colormap
        if cmap is None:
            cmap = 'binary'

        # Plot extinction map
        if self.coord_system == 'icrs':
            ax.imshow(
                np.sqrt(av)[::, ::-1],
                vmin=vmin,
                vmax=vmax,
                origin='lower',
                interpolation='bilinear',
                cmap=cmap,
                aspect='equal',
                extent=[
                    self.ra_max * 15.0,
                    self.ra_min * 15.0,
                    self.dec_min,
                    self.dec_max
                ],
            )
            ax.set_xlim(self.ra_max * 15.0, self.ra_min * 15.0)
            ax.set_ylim(self.dec_min, self.dec_max)
            # Label offsets are handled in _apply_scaling, use the returned value
            ax.set_xlabel('RA (J2000)', labelpad=labelpad)
            ax.set_ylabel('Dec (J2000)', labelpad=labelpad)

            formatter = matplotlib.ticker.FuncFormatter(
                lambda ra, x: time.strftime(
                    '%Hh %M',
                    time.gmtime((ra / 15.0) * 3600.0)
                )
            )
            ax.xaxis.set_major_formatter(formatter)
            # Grid linewidth will be scaled in _apply_scaling if interactive
            plt.grid(axis='x', color='0.3', linestyle='dashed', alpha=0.3, linewidth=0.3)
            plt.grid(axis='y', color='0.3', linestyle='dashed', alpha=0.3, linewidth=0.3)
        else:  # galactic
            ax.imshow(
                np.sqrt(av)[::, ::-1],
                vmin=vmin,
                vmax=vmax,
                origin='lower',
                interpolation='bilinear',
                cmap=cmap,
                aspect='equal',
                extent=[
                    self.l_max,
                    self.l_min,
                    self.b_min,
                    self.b_max
                ]
            )
            ax.set_xlim(self.l_min, self.l_max)
            ax.set_ylim(self.b_min, self.b_max)
            # Label offsets are handled in _apply_scaling, use the returned value
            ax.set_xlabel('Galactic longitude', labelpad=labelpad)
            ax.set_ylabel('Galactic latitude', labelpad=labelpad)
            ax.set_aspect('equal')

        # Plot sources
        # Scale markers by the same factor as font size for interactive plots
        marker_scale = scale_factor**2 if interactive else 1.0

        # Plot discs first (stored for interactive labels), then PMS so that
        # discs can still be drawn above PMS via zorder and interactive labels
        # are set up once for both.
        if plot_discs:
            self.plot_all_discs(
                ax,
                csvfile=discs_csvfile,
                interactive=interactive,
                font_scale=scale_factor,
                marker_scale=marker_scale,
            )

        if plot_pms:
            self.plot_kenyon08_pms(
                ax,
                csvfile=pms_csvfile,
                only_label_famous=only_label_famous,
                interactive=interactive,
                font_scale=scale_factor,
                marker_scale=marker_scale,
            )

        # Set title with scaled padding and font size for interactive plots
        if interactive:
            title_pad = getattr(self, '_title_pad', 4.0)
            scaled_fontsize = getattr(self, '_scaled_fontsize', None)
            if scaled_fontsize:
                ax.set_title(title, pad=title_pad, fontsize=scaled_fontsize)
            else:
                ax.set_title(title, pad=title_pad)
        else:
            ax.set_title(title, pad=4.0)

        # Adjust layout to ensure title and labels are visible
        # For interactive plots, use tight_layout with minimal padding
        if interactive:
            fig.tight_layout(pad=0.3, rect=[0, 0, 1, 0.98])  # Minimal padding, leave small space for title
        else:
            fig.subplots_adjust(wspace=0.0, hspace=0.0, top=0.95, bottom=0.1, left=0.1, right=0.95)

        if save_path:
            # Save with original figure size and DPI for high quality
            # Temporarily resize figure and restore font sizes for saving
            original_figsize = fig.get_size_inches()
            if interactive:
                # Restore original font sizes and line widths for saving
                import matplotlib as mpl
                default_fontsize = mpl.rcParams['font.size']
                ax.tick_params(labelsize=default_fontsize, width=1.0, length=4.0)
                ax.xaxis.label.set_fontsize(default_fontsize)
                ax.yaxis.label.set_fontsize(default_fontsize)
                ax.title.set_fontsize(default_fontsize)
                # Restore label offsets
                ax.set_xlabel(ax.get_xlabel(), labelpad=3)
                ax.set_ylabel(ax.get_ylabel(), labelpad=3)
                # Restore grid line widths
                for line in ax.xaxis.get_gridlines():
                    line.set_linewidth(0.3)
                for line in ax.yaxis.get_gridlines():
                    line.set_linewidth(0.3)
                # Restore axis spine line widths
                for spine in ax.spines.values():
                    spine.set_linewidth(0.5)
                # Restore subplots_adjust for saved plots
                fig.subplots_adjust(wspace=0.0, hspace=0.0)

            fig.set_size_inches(figsize)
            plt.savefig(save_path, bbox_inches='tight', dpi=dpi)
            if save_path.endswith('.pdf'):
                # Also save PNG version
                png_path = save_path.replace('.pdf', '.png')
                plt.savefig(png_path, dpi=dpi, bbox_inches='tight')

            # Restore display size and scaled fonts/line widths if interactive
            if interactive:
                fig.set_size_inches(display_figsize)
                # Re-apply scaling for interactive display
                scaled_fontsize, labelpad = self._apply_scaling(ax, scale_factor, default_fontsize)
                ax.set_xlabel(ax.get_xlabel(), labelpad=labelpad)
                ax.set_ylabel(ax.get_ylabel(), labelpad=labelpad)
                title_pad = getattr(self, '_title_pad', 4.0)
                ax.set_title(ax.get_title(), pad=title_pad, fontsize=scaled_fontsize)
                # Restore grid line widths
                grid_linewidth = 0.3 * scale_factor
                for line in ax.xaxis.get_gridlines():
                    line.set_linewidth(grid_linewidth)
                for line in ax.yaxis.get_gridlines():
                    line.set_linewidth(grid_linewidth)
                # Restore tight layout for interactive
                fig.tight_layout(pad=0.3, rect=[0, 0, 1, 0.98])

        if interactive:
            plt.show()

        return fig, ax

