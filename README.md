# Discs in Context

A Python wrapper for plotting extinction maps of star-forming regions via the dustmaps package, with protostars and protoplanetary discs overlaid, e.g.:

<img width="990" height="917" alt="image" src="https://github.com/user-attachments/assets/a5c97092-34d3-476c-84bf-06d920843b80" />

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
- reproject (for ``dustmap='herschel'``; SciPy fallback if missing)
- astroquery (ESASky SPIRE download for Herschel)
- scipy

### Downloading Dust Maps

After installing, you'll need to download the dust map data using a short python script:

```python
from dustmaps.planck import fetch as fp
from dustmaps.sfd import fetch as fs
from dustmaps.bayestar import fetch as fb

fp()
fs()
fb()
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
    dustmap='planck',  # or 'sfd', 'bayestar', 'herschel', 'herschel_hips'
    figsize=(18, 10),
    dpi=300,
    vmin=0.0,
    vmax=4.0,
    cmap='inferno',  # colormap name
    colorbar=True,  # draw colour bar (default True)
    stretch='linear',  # or 'sqrt' / 'log' (Herschel defaults to log)
    plot_discs=False,  # plot all discs from CSV
    plot_pms=True,  # plot PMS sources
    pms_csvfile='tau-sources.csv',  # custom PMS file
    discs_csvfile='discs.csv',  # custom discs file
    save_path='output.pdf',
    show=False  # display interactively
)
```

### Herschel SPIRE background

**Pointed maps** (`dustmap='herschel'`): science-grade Level-2 SPIRE via ESASky (or a local FITS path).

```python
obs = pc(
    object='IRAS 18148-0440',
    image_size=1.0,
    image_size_unit='degrees',
    coord_system='icrs',
)
obs.plot(
    dustmap='herschel',
    herschel_band='psw',          # 250 um; or 'pmw' / 'plw'
    herschel_cache_dir='herschel_l483',  # optional cache path
    # herschel_fits='path/to/spire.fits.gz',  # skip download
    vmin=60,
    vmax=4000,
    interactive=True,
    plot_discs=True,
    plot_halpha=True,
)
```

**HiPS mosaic** (`dustmap='herschel_hips'`): ESA all-observed-sky SPIRE HiPS cutout via hips2fits. Convenient anywhere Herschel looked; lower fidelity than pointed downloads.

```python
obs.plot(
    dustmap='herschel_hips',
    herschel_band='psw',
    vmin=60,
    vmax=4000,
    interactive=True,
    plot_discs=True,
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

