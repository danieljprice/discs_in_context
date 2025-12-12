# Discs in Context

A Python wrapper for plotting extinction maps of star-forming regions via the dustmaps package, with protostars and protoplanetary discs overlaid.

## Installation

### Install the Package

Install in editable mode (recommended for development):

```bash
pip install -e .
```

Or install normally:

```bash
pip install .
```

### Dependencies

This package requires the following dependencies (automatically installed with the package):
- matplotlib
- numpy
- astropy
- dustmaps (with Planck, SFD, and/or Bayestar maps downloaded)
- pandas

### Downloading Dust Maps

After installing, you'll need to download the dust map data using the `dustmaps` command-line tools:

```bash
dustmaps_download planck
dustmaps_download sfd
dustmaps_download bayestar
```

## Usage

### Basic Usage with Preset Regions

```python
from discs_in_context import plotcloud as pc

# Plot Taurus region in galactic coordinates
plotter = pc(region='taurus', coord_system='galactic')
fig, ax = plotter.plot(save_path='taurus.pdf', plot_pms=True)
```

### Using RA/Dec Coordinates

```python
plotter = pc(
    ra_range=(4.0, 5.3),  # RA in hours
    dec_range=(17.5, 31.5),  # Dec in degrees
    coord_system='icrs'
)
fig, ax = plotter.plot(save_path='region_icrs.pdf')
```

### Using Galactic Coordinates

```python
plotter = pc(
    galactic_l={'l0': 173.1565, 'lsize': 7.5},
    galactic_b={'b0': -15.9071, 'bsize': 7.5},
    coord_system='galactic'
)
fig, ax = plotter.plot(save_path='region_galactic.pdf')
```

### Available Preset Regions

- `taurus` - Taurus-Auriga molecular cloud
- `lupus` - Lupus star-forming region
- `orion` - Orion molecular cloud
- `ophiuchus` - Ophiuchus star-forming region
- `chamaeleon` - Chamaeleon molecular cloud
- `allsky` - Full sky view (0-360° galactic longitude, -90 to 90° latitude)

### Plotting Options

The `plot()` method accepts many customization options:

```python
plotter.plot(
    dustmap='planck',  # or 'sfd', 'bayestar'
    figsize=(18, 10),
    dpi=300,
    vmin=0.0,
    vmax=2.0,
    cmap='inferno',  # colormap name
    plot_discs=False,  # plot all discs from CSV
    plot_pms=True,  # plot PMS sources
    pms_csvfile='tau-sources.csv',  # custom PMS file
    discs_csvfile='discs.csv',  # custom discs file
    save_path='output.pdf',
    show=False  # display interactively
)
```

## Class API

### plotcloud

Main plotting class.

#### Initialization

- `region` (str): Preset region name
- `ra_range` (tuple): (ra_min, ra_max) in hours
- `dec_range` (tuple): (dec_min, dec_max) in degrees
- `galactic_l` (dict): {'l_min': min_l, 'l_max': max_l}
- `galactic_b` (dict): {'b_min': min_b, 'b_max': max_b}
- `coord_system` (str): 'galactic' or 'icrs'
- `num_points` (int): Grid resolution (default 2048)

#### Methods

- `plot()`: Create the extinction map plot
- `plot_all_discs()`: Plot protoplanetary discs from CSV
- `plot_kenyon08_pms()`: Plot PMS sources from Kenyon 2008 catalog

