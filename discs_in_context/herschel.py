"""
Herschel SPIRE helpers for discs_in_context.

Two backends:
- Pointed Level-2 maps via ESASky (``dustmap='herschel'``)
- ESA SPIRE HiPS cutouts via hips2fits (``dustmap='herschel_hips'``)

Both return surface brightness on a plotcloud coordinate grid (MJy/sr).
"""

from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

# SPIRE band short names used in ESA product filenames / HiPS ids
VALID_BANDS = ('psw', 'pmw', 'plw')  # 250, 350, 500 um

# ESA ESDC Herschel SPIRE HiPS (visualisation mosaic of public maps)
HIPS_IDS = {
    'psw': 'ESAVO/P/HERSCHEL/SPIRE-250',
    'pmw': 'ESAVO/P/HERSCHEL/SPIRE-350',
    'plw': 'ESAVO/P/HERSCHEL/SPIRE-500',
}


def default_cache_dir():
    """Default on-disk cache for ESASky SPIRE downloads."""
    return Path.home() / '.cache' / 'discs_in_context' / 'herschel'


def open_spire_image(fits_path):
    """
    Return image data, WCS, and BUNIT for a SPIRE map product.

    Extension 1 is the science image in standard ESASky/HSA SPIRE products.
    """
    fits_path = Path(fits_path)
    with fits.open(fits_path) as hdul:
        hdu = hdul[1]
        data = np.array(hdu.data, dtype=float)
        wcs = WCS(hdu.header)
        bunit = hdu.header.get('BUNIT', 'MJy/sr')
    return data, wcs, bunit


def ensure_spire_download(position, download_dir=None, band='psw',
                          radius=30 * u.arcmin):
    """
    Query ESASky for SPIRE maps near ``position`` and cache them on disk.

    Parameters
    ----------
    position : str or SkyCoord
        Object name or sky position passed to ESASky.
    download_dir : path-like, optional
        Cache directory. Defaults to ``~/.cache/discs_in_context/herschel``.
    band : str, default 'psw'
        SPIRE band: 'psw' (250 um), 'pmw' (350 um), or 'plw' (500 um).
    radius : Quantity, default 30 arcmin
        Cone-search radius for ESASky.

    Returns
    -------
    Path
        Path to the cached SPIRE FITS (gzipped) for the requested band.
    """
    band = str(band).lower()
    if band not in VALID_BANDS:
        raise ValueError(
            f"herschel_band must be one of {VALID_BANDS}, got '{band}'"
        )

    download_dir = Path(download_dir) if download_dir is not None else default_cache_dir()
    download_dir.mkdir(parents=True, exist_ok=True)

    # Reuse any previously downloaded product for this band
    existing = list(download_dir.glob(f'**/hspire{band}*.fits.gz'))
    if existing:
        return existing[0]

    from astroquery.esasky import ESASky
    from astroquery.utils import TableList

    if isinstance(position, SkyCoord):
        query_pos = position
    else:
        query_pos = str(position)

    print(f"Querying ESASky for Herschel SPIRE near {query_pos}...")
    maps = ESASky.query_region_maps(
        position=query_pos,
        radius=radius,
        missions=['HERSCHEL'],
    )
    if 'HERSCHEL' not in maps.keys():
        raise ValueError(f"No Herschel maps found near {query_pos}")

    table = maps['HERSCHEL']
    spire = table[table['instrument'] == 'SPIRE']
    if len(spire) == 0:
        raise ValueError(f"No Herschel SPIRE maps found near {query_pos}")

    print(f"Downloading {len(spire)} SPIRE observation(s) to {download_dir}...")
    ESASky.get_maps(
        TableList({'HERSCHEL': spire}),
        missions=['HERSCHEL'],
        download_dir=str(download_dir),
    )

    existing = list(download_dir.glob(f'**/hspire{band}*.fits.gz'))
    if not existing:
        raise ValueError(
            f"Download finished but no hspire{band}*.fits.gz found under {download_dir}"
        )
    return existing[0]


def _wcs_from_skycoord_grid(coords):
    """
    Build a simple TAN WCS matching a regular SkyCoord meshgrid.

    Pixel (1, 1) is the centre of ``coords[0, 0]``; axis 1 follows RA/l,
    axis 2 follows Dec/b.
    """
    ny, nx = coords.shape
    header = fits.Header()
    header['NAXIS'] = 2
    header['NAXIS1'] = nx
    header['NAXIS2'] = ny
    header['CRPIX1'] = 1.0
    header['CRPIX2'] = 1.0
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'

    if hasattr(coords, 'ra'):
        lon = coords.ra.degree
        lat = coords.dec.degree
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
    else:
        lon = coords.l.degree
        lat = coords.b.degree
        header['CTYPE1'] = 'GLON-TAN'
        header['CTYPE2'] = 'GLAT-TAN'

    header['CRVAL1'] = float(lon[0, 0])
    header['CRVAL2'] = float(lat[0, 0])
    if nx > 1:
        header['CDELT1'] = float(lon[0, -1] - lon[0, 0]) / (nx - 1)
    else:
        header['CDELT1'] = 1.0 / 3600.0
    if ny > 1:
        header['CDELT2'] = float(lat[-1, 0] - lat[0, 0]) / (ny - 1)
    else:
        header['CDELT2'] = 1.0 / 3600.0

    return WCS(header)


def _fov_radius(coords):
    """Half-diagonal of the FOV as an angular Quantity (floor 30 arcmin)."""
    if hasattr(coords, 'ra'):
        c00 = coords[0, 0]
        c11 = coords[-1, -1]
    else:
        c00 = coords[0, 0]
        c11 = coords[-1, -1]
    sep = c00.separation(c11) / 2.0
    return max(sep, 30 * u.arcmin)


def herschel_on_grid(coords, fits_path=None, position=None, band='psw',
                     cache_dir=None):
    """
    Return SPIRE surface brightness on the same grid as ``coords``.

    Parameters
    ----------
    coords : SkyCoord array
        2D coordinate grid from plotcloud (ICRS or galactic).
    fits_path : path-like, optional
        Local SPIRE FITS. If None, download via ESASky.
    position : str or SkyCoord, optional
        Query centre for ESASky when ``fits_path`` is None.
        Defaults to the FOV centre.
    band : str, default 'psw'
        SPIRE band short name when downloading.
    cache_dir : path-like, optional
        ESASky download cache directory.

    Returns
    -------
    map_2d : ndarray
        Reprojected intensity (MJy/sr); NaN outside the SPIRE footprint.
    bunit : str
        Unit string from the FITS header.
    """
    if fits_path is None:
        if position is None:
            # FOV centre as ICRS for the ESASky cone search
            mid = coords[coords.shape[0] // 2, coords.shape[1] // 2]
            position = mid.icrs if hasattr(mid, 'icrs') else mid
        radius = _fov_radius(coords)
        fits_path = ensure_spire_download(
            position,
            download_dir=cache_dir,
            band=band,
            radius=radius,
        )

    data, src_wcs, bunit = open_spire_image(fits_path)
    map_2d = _reproject_to_coords(data, src_wcs, coords)
    return map_2d, bunit


def herschel_hips_on_grid(coords, band='psw', hips_id=None, min_coverage=1e-4):
    """
    Return ESA Herschel SPIRE HiPS intensity on the ``coords`` grid.

    Uses CDS hips2fits with a WCS matched to the plotcloud meshgrid.
    Empty / unobserved sky is NaN. Raises ``ValueError`` if essentially
    no HiPS coverage falls in the FOV.

    Parameters
    ----------
    coords : SkyCoord array
        2D coordinate grid from plotcloud.
    band : str, default 'psw'
        SPIRE band: 'psw' (250 um), 'pmw' (350), 'plw' (500).
    hips_id : str, optional
        Override HiPS identifier (default from ``HIPS_IDS[band]``).
    min_coverage : float, default 1e-4
        Minimum fraction of finite pixels required; below this, raise.

    Returns
    -------
    map_2d : ndarray
        HiPS intensity cutout (typically MJy/sr); NaN outside coverage.
    bunit : str
        Unit label for the colour bar.
    """
    from astroquery.hips2fits import hips2fits

    band = str(band).lower()
    if hips_id is None:
        if band not in HIPS_IDS:
            raise ValueError(
                f"herschel_band must be one of {VALID_BANDS}, got '{band}'"
            )
        hips_id = HIPS_IDS[band]

    target_wcs = _wcs_from_skycoord_grid(coords)
    print(f"Fetching Herschel HiPS cutout ({hips_id})...")
    hdul = hips2fits.query_with_wcs(
        hips=hips_id,
        wcs=target_wcs,
        format='fits',
    )
    data = np.array(hdul[0].data, dtype=float)
    if data.shape != coords.shape:
        # Rare mismatch: reproject returned cutout onto the plot grid
        map_2d = _reproject_to_coords(data, WCS(hdul[0].header), coords)
    else:
        map_2d = data

    finite_frac = float(np.isfinite(map_2d).mean()) if map_2d.size else 0.0
    if finite_frac < min_coverage:
        raise ValueError(
            f"No Herschel HiPS coverage in this FOV for '{hips_id}' "
            f"(finite fraction={finite_frac:.2e})."
        )

    bunit = hdul[0].header.get('BUNIT', 'MJy/sr')
    return map_2d, bunit


def _reproject_to_coords(data, src_wcs, coords):
    """Reproject image onto coords grid (reproject if available, else scipy)."""
    try:
        from reproject import reproject_interp

        target_wcs = _wcs_from_skycoord_grid(coords)
        map_2d, _footprint = reproject_interp(
            (data, src_wcs),
            target_wcs,
            shape_out=coords.shape,
        )
        return map_2d
    except ImportError:
        pass

    # Fallback: sample the SPIRE image at each plotcloud sky position
    from scipy.ndimage import map_coordinates

    world = coords.icrs if hasattr(coords, 'ra') else coords
    x, y = src_wcs.world_to_pixel(world)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ny, nx = data.shape
    outside = (
        (x < 0) | (x > nx - 1) | (y < 0) | (y > ny - 1)
        | ~np.isfinite(x) | ~np.isfinite(y)
    )

    data_f = np.array(data, dtype=float)
    valid = np.isfinite(data_f)
    data_f = np.where(valid, data_f, 0.0)

    map_2d = map_coordinates(
        data_f, [y, x], order=1, mode='constant', cval=np.nan, prefilter=False,
    )
    footprint = map_coordinates(
        valid.astype(float), [y, x], order=0, mode='constant', cval=0.0,
    )
    map_2d = np.asarray(map_2d, dtype=float)
    map_2d[footprint < 0.5] = np.nan
    map_2d[outside] = np.nan
    return map_2d
