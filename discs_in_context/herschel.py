"""
Far-infrared imaging helpers for discs_in_context.

Backends:
- Pointed Herschel SPIRE Level-2 maps via ESASky (``dustmap='herschel'``)
- ESA SPIRE HiPS cutouts via hips2fits (``dustmap='herschel_hips'``)
- AKARI FIS HiPS (``dustmap='akari_hips'``)
- IRAS 100 um via SkyView IRIS (``dustmap='iras'``)

All return surface brightness on a plotcloud coordinate grid (MJy/sr).
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

# AKARI FIS all-sky HiPS (useful near ~160 um when Herschel has no coverage)
AKARI_VALID_BANDS = ('n160', 'widel')
AKARI_HIPS_IDS = {
    'n160': 'CDS/P/AKARI/FIS/N160',    # narrow 160 um
    'widel': 'CDS/P/AKARI/FIS/WideL',  # Wide-L (~140-180 um)
}

# NASA SkyView survey names for IRIS (improved IRAS; all-sky, holes filled)
IRAS_VALID_BANDS = ('12', '25', '60', '100')
IRAS_SKYVIEW_SURVEYS = {
    '12': 'IRIS  12',
    '25': 'IRIS  25',
    '60': 'IRIS  60',
    '100': 'IRIS 100',
}


def default_cache_dir():
    """Default on-disk cache for ESASky SPIRE downloads."""
    return Path.home() / '.cache' / 'discs_in_context' / 'herschel'


def default_hips_cache_dir():
    """Default on-disk cache for HiPS cutouts (Herschel / AKARI)."""
    return Path.home() / '.cache' / 'discs_in_context' / 'hips_cutouts'


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


def _as_skycoord(position):
    """Resolve a name or SkyCoord to an ICRS SkyCoord."""
    if isinstance(position, SkyCoord):
        return position.icrs if hasattr(position, 'icrs') else position
    return SkyCoord.from_name(str(position))


def _spire_covers_position(fits_path, position, margin_pix=1.0):
    """
    Return True if ``position`` falls inside the SPIRE image footprint.

    Used so a shared download cache never reuses a map of the wrong field
    (e.g. V380 Ori for V883 Ori).
    """
    try:
        data, wcs, _bunit = open_spire_image(fits_path)
    except Exception:
        return False
    pos = _as_skycoord(position)
    try:
        x, y = wcs.world_to_pixel(pos)
    except Exception:
        return False
    x = float(np.asarray(x).reshape(-1)[0])
    y = float(np.asarray(y).reshape(-1)[0])
    if not (np.isfinite(x) and np.isfinite(y)):
        return False
    ny, nx = data.shape
    return (
        margin_pix <= x < (nx - 1 - margin_pix)
        and margin_pix <= y < (ny - 1 - margin_pix)
    )


def _pick_covering_spire(paths, position):
    """
    Prefer a cached SPIRE FITS that covers ``position``.

    Among covering maps, prefer the largest file (often a mosaic / large scan).
    """
    covering = [p for p in paths if _spire_covers_position(p, position)]
    if not covering:
        return None
    return max(covering, key=lambda p: p.stat().st_size)


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
        Cone-search radius for ESASky. Capped at 45 arcmin so a wide plot
        FOV does not download every SPIRE map in several degrees.

    Returns
    -------
    Path
        Path to a cached SPIRE FITS (gzipped) for the requested band that
        covers ``position``.
    """
    band = str(band).lower()
    if band not in VALID_BANDS:
        raise ValueError(
            f"herschel_band must be one of {VALID_BANDS}, got '{band}'"
        )

    download_dir = Path(download_dir) if download_dir is not None else default_cache_dir()
    download_dir.mkdir(parents=True, exist_ok=True)

    position_coord = _as_skycoord(position)

    # Reuse a cached product only if it actually covers this target
    existing = list(download_dir.glob(f'**/hspire{band}*.fits.gz'))
    picked = _pick_covering_spire(existing, position_coord)
    if picked is not None:
        print(f"Using cached SPIRE {band} map covering target: {picked.name}")
        return picked

    from astroquery.esasky import ESASky
    from astroquery.utils import TableList

    # Limit the cone so wide FOVs do not trigger multi-degree mass downloads
    query_radius = min(radius, 45 * u.arcmin)

    print(
        f"Querying ESASky for Herschel SPIRE near {position_coord.icrs.to_string('hmsdms')} "
        f"(r={query_radius})..."
    )
    maps = ESASky.query_region_maps(
        position=position_coord,
        radius=query_radius,
        missions=['HERSCHEL'],
    )
    if 'HERSCHEL' not in maps.keys():
        raise ValueError(f"No Herschel maps found near {position}")

    table = maps['HERSCHEL']
    spire = table[table['instrument'] == 'SPIRE']
    if len(spire) == 0:
        raise ValueError(f"No Herschel SPIRE maps found near {position}")

    # Prefer maps whose centres are closest to the target
    ra = np.asarray(spire['ra_deg'], dtype=float)
    dec = np.asarray(spire['dec_deg'], dtype=float)
    centres = SkyCoord(ra * u.deg, dec * u.deg, frame='icrs')
    sep = position_coord.separation(centres)
    spire = spire[np.argsort(sep)]

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

    picked = _pick_covering_spire(existing, position_coord)
    if picked is None:
        raise ValueError(
            f"Downloaded SPIRE {band} maps near {position}, but none cover the "
            f"target position. Try dustmap='herschel_hips' for mosaic coverage."
        )
    print(f"Selected SPIRE {band} map: {picked.name}")
    return picked


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
                     cache_dir=None, min_coverage=0.05):
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
    min_coverage : float, default 0.05
        If a pointed map covers less than this fraction of the FOV,
        fall back to the ESA SPIRE HiPS mosaic (needed for wide fields
        such as V883 Ori at a few degrees).

    Returns
    -------
    map_2d : ndarray
        Reprojected intensity (MJy/sr); NaN outside the SPIRE footprint.
    bunit : str
        Unit string from the FITS header.
    """
    if position is None:
        mid = coords[coords.shape[0] // 2, coords.shape[1] // 2]
        position = mid.icrs if hasattr(mid, 'icrs') else mid

    if fits_path is None:
        radius = _fov_radius(coords)
        try:
            fits_path = ensure_spire_download(
                position,
                download_dir=cache_dir,
                band=band,
                radius=radius,
            )
        except ValueError as exc:
            print(
                f"Pointed SPIRE download failed ({exc}); "
                f"falling back to Herschel SPIRE HiPS."
            )
            return herschel_hips_on_grid(coords, band=band)

    data, src_wcs, bunit = open_spire_image(fits_path)
    map_2d = _reproject_to_coords(data, src_wcs, coords)
    finite_frac = float(np.isfinite(map_2d).mean()) if map_2d.size else 0.0
    if finite_frac < min_coverage:
        print(
            f"Pointed SPIRE covers only {finite_frac:.1%} of the FOV "
            f"({Path(fits_path).name}); falling back to Herschel SPIRE HiPS."
        )
        return herschel_hips_on_grid(coords, band=band)
    return map_2d, bunit


def _hips_on_grid(coords, hips_id, min_coverage=1e-4, label='HiPS',
                  cache_dir=None, use_cache=True):
    """
    Fetch a HiPS FITS cutout via hips2fits onto the ``coords`` grid.

    Empty / unobserved sky is NaN. Raises ``ValueError`` if essentially
    no HiPS coverage falls in the FOV.

    Cutouts are cached under ``~/.cache/discs_in_context/hips_cutouts/``
    keyed by HiPS id and WCS (FOV / pixel grid), so replotting the same
    field does not re-hit the network.
    """
    import hashlib
    from astroquery.hips2fits import hips2fits

    target_wcs = _wcs_from_skycoord_grid(coords)
    hdr = target_wcs.to_header()
    # Stable key from survey + FOV geometry (not the full float WCS dump)
    key_parts = [
        str(hips_id),
        str(coords.shape[0]), str(coords.shape[1]),
        f"{hdr.get('CRVAL1', 0):.6f}", f"{hdr.get('CRVAL2', 0):.6f}",
        f"{hdr.get('CDELT1', 0):.8e}", f"{hdr.get('CDELT2', 0):.8e}",
        str(hdr.get('CTYPE1', '')), str(hdr.get('CTYPE2', '')),
    ]
    digest = hashlib.sha1('|'.join(key_parts).encode('utf-8')).hexdigest()[:16]
    safe_id = str(hips_id).replace('/', '_').replace(' ', '_')
    cache_root = Path(cache_dir) if cache_dir is not None else default_hips_cache_dir()
    cache_path = cache_root / f'{safe_id}_{digest}.fits.gz'

    if use_cache and cache_path.is_file():
        print(f"Using cached {label} cutout: {cache_path.name}")
        with fits.open(cache_path) as hdul:
            data = np.array(hdul[0].data, dtype=float)
            src_wcs = WCS(hdul[0].header)
            bunit = hdul[0].header.get('BUNIT', 'MJy/sr')
    else:
        print(f"Fetching {label} cutout ({hips_id})...")
        hdul = hips2fits.query_with_wcs(
            hips=hips_id,
            wcs=target_wcs,
            format='fits',
        )
        data = np.array(hdul[0].data, dtype=float)
        src_wcs = WCS(hdul[0].header)
        bunit = hdul[0].header.get('BUNIT', 'MJy/sr')
        if use_cache:
            cache_root.mkdir(parents=True, exist_ok=True)
            out_hdr = fits.Header(hdul[0].header)
            out_hdr['BUNIT'] = bunit
            out_hdr['HIPSID'] = str(hips_id)
            fits.PrimaryHDU(data=data, header=out_hdr).writeto(
                cache_path, overwrite=True,
            )
            print(f"Cached {label} cutout -> {cache_path}")

    if data.shape != coords.shape:
        map_2d = _reproject_to_coords(data, src_wcs, coords)
    else:
        map_2d = data

    finite_frac = float(np.isfinite(map_2d).mean()) if map_2d.size else 0.0
    if finite_frac < min_coverage:
        raise ValueError(
            f"No {label} coverage in this FOV for '{hips_id}' "
            f"(finite fraction={finite_frac:.2e})."
        )

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
    band = str(band).lower()
    if hips_id is None:
        if band not in HIPS_IDS:
            raise ValueError(
                f"herschel_band must be one of {VALID_BANDS}, got '{band}'"
            )
        hips_id = HIPS_IDS[band]

    return _hips_on_grid(
        coords, hips_id, min_coverage=min_coverage, label='Herschel HiPS',
    )


def akari_hips_on_grid(coords, band='n160', hips_id=None, min_coverage=1e-4):
    """
    Return AKARI FIS HiPS intensity on the ``coords`` grid.

    AKARI provides all-sky ~160 um imaging, useful when Herschel PACS/SPIRE
    never observed the field (e.g. GW Ori / lambda Ori).

    Parameters
    ----------
    coords : SkyCoord array
        2D coordinate grid from plotcloud.
    band : str, default 'n160'
        'n160' (narrow 160 um) or 'widel' (Wide-L).
    hips_id : str, optional
        Override HiPS identifier (default from ``AKARI_HIPS_IDS[band]``).
    min_coverage : float, default 1e-4
        Minimum fraction of finite pixels required; below this, raise.

    Returns
    -------
    map_2d : ndarray
        HiPS intensity cutout (typically MJy/sr); NaN outside coverage.
    bunit : str
        Unit label for the colour bar.
    """
    band = str(band).lower()
    if hips_id is None:
        if band not in AKARI_HIPS_IDS:
            raise ValueError(
                f"akari_band must be one of {AKARI_VALID_BANDS}, got '{band}'"
            )
        hips_id = AKARI_HIPS_IDS[band]

    return _hips_on_grid(
        coords, hips_id, min_coverage=min_coverage, label='AKARI HiPS',
    )


def iras_on_grid(coords, band='100', pixels=512):
    """
    Return IRAS/IRIS surface brightness on the ``coords`` grid.

    Uses NASA SkyView's IRIS maps (Miville-Deschenes & Lagache 2005):
    all-sky IRAS photometry with improved calibration and DIRBE fill
    in unobserved strips.

    Parameters
    ----------
    coords : SkyCoord array
        2D coordinate grid from plotcloud.
    band : str, default '100'
        IRAS band in microns: '12', '25', '60', or '100'.
    pixels : int, default 512
        SkyView cutout size (native IRAS resolution is ~4 arcmin).

    Returns
    -------
    map_2d : ndarray
        Intensity on the plot grid (MJy/sr).
    bunit : str
        Unit label for the colour bar.
    """
    from astroquery.skyview import SkyView

    band = str(band).lower().replace('um', '').replace('micron', '').strip()
    if band not in IRAS_SKYVIEW_SURVEYS:
        raise ValueError(
            f"iras_band must be one of {IRAS_VALID_BANDS}, got '{band}'"
        )

    ny, nx = coords.shape
    mid = coords[ny // 2, nx // 2]
    centre = mid.icrs if hasattr(mid, 'icrs') else mid
    radius = _fov_radius(coords)

    print(f"Fetching IRAS/IRIS {band} um from SkyView (bilinear)...")
    imgs = SkyView.get_images(
        position=centre,
        survey=IRAS_SKYVIEW_SURVEYS[band],
        radius=radius,
        pixels=int(pixels),
        sampler='LI',  # bilinear resampling of the native IRAS pixels
    )
    hdu = imgs[0][0]
    data = np.array(hdu.data, dtype=float)
    src_wcs = WCS(hdu.header).celestial
    bunit = hdu.header.get('BUNIT', 'MJy/sr')
    if isinstance(bunit, str) and bunit.upper() == 'MJY/SR':
        bunit = 'MJy/sr'
    map_2d = _reproject_to_coords(data, src_wcs, coords)
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
            order=1,  # bilinear
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
