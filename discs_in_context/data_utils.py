from importlib import resources


def _resource_path(filename: str) -> str:
    """Return a filesystem path to a data file bundled with the package.

    Parameters
    ----------
    filename : str
        File name inside the ``discs_in_context/data`` directory.
    """
    csv_res = resources.files("discs_in_context") / "data" / filename  # type: ignore[operator]
    return str(csv_res)


def get_tau_sources_path() -> str:
    """Return a filesystem path to the bundled ``tau-sources.csv``."""
    return _resource_path("tau-sources.csv")


def load_tau_sources_text() -> str:
    """Return the contents of the bundled ``tau-sources.csv`` as text."""
    with open(get_tau_sources_path(), "r", encoding="utf-8") as f:
        return f.read()


def get_discs_catalog_path() -> str:
    """Return a filesystem path to the bundled discs catalogue CSV.

    By default this looks for ``circumstellar-disk-catalog-with-gaiadr3.csv``
    inside the ``discs_in_context/data`` directory. To use a different
    catalogue, either pass ``csvfile=...`` to the relevant plotting
    methods or replace this file in the data directory.
    """
    return _resource_path("circumstellar-disk-catalog-with-gaiadr3.csv")


def get_scocen_catalog_path() -> str:
    """Return a filesystem path to the bundled Sco-Cen catalogue CSV.

    By default this looks for ``Sco_Cen_all.csv`` inside the
    ``discs_in_context/data`` directory. To use a different catalogue,
    pass ``scocen_csvfile=...`` to ``plot()``.
    """
    return _resource_path("Sco_Cen_all.csv")


def get_halpha_catalog_path() -> str:
    """Return a filesystem path to the bundled Halpha sources catalogue CSV.

    By default this looks for ``Delfini_strong_Halpha_subset.csv`` inside the
    ``discs_in_context/data`` directory. To use a different catalogue,
    pass ``halpha_csvfile=...`` to ``plot()``.
    """
    return _resource_path("Delfini_strong_Halpha_subset.csv")

