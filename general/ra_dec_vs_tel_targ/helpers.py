import os
from concurrent.futures import ThreadPoolExecutor
from astropy.table import Table, vstack
from matplotlib import pyplot as plt
from astropy.coordinates import FK5, SkyCoord
import astropy.units as u
import numpy as np
import fnmatch
from astroquery.eso import Eso # import the ESO module from astroquery
from astroquery.simbad import Simbad
from tqdm import tqdm

# Global plot style settings
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.weight"] = "bold"
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in" 
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.linewidth'] = 0.5
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['axes.titleweight'] = 'bold'   # normal, bold, or a numeric weight
plt.rcParams['axes.titlesize']   = 12       # or 'large', 'medium', etc.
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.linestyle'] = ':'

def build_main_teltarg_table(collection_name, maxrec=1, max_workers=8, overwrite=False):
    eso = Eso()
    eso.ROW_LIMIT = maxrec

    output_file = f"main_teltarg_table_{collection_name}.fits"

    if os.path.exists(output_file) and not overwrite:
        print(f"{output_file} already exists, skipping query.")
        final_table = Table.read(output_file)
        return final_table, len(final_table)

    if overwrite and os.path.exists(output_file):
        print(f"{output_file} already exists, overwriting.")
    else:
        print(f"{output_file} does not exist, performing query.")

    table = eso.query_surveys(collection_name)
    N = len(table)

    header_vals = fetch_header_values_parallel(
        eso,
        table["dp_id"],
        desc="Fetching headers",
        max_workers=max_workers,
        pattern_keys=[
            ("ra_teltarg_raw", "HIERARCH ESO TEL* TARG ALPHA"),
            ("dec_teltarg_raw", "HIERARCH ESO TEL* TARG DELTA"),
        ],
    )

    coord_teltarg = packed_ra_dec_to_skycoord(
        header_vals["ra_teltarg_raw"],
        header_vals["dec_teltarg_raw"],
    )

    ra_main = np.asarray(table["s_ra"], dtype=float)
    dec_main = np.asarray(table["s_dec"], dtype=float)
    ra_teltarg = np.asarray(coord_teltarg.ra.deg, dtype=float)
    dec_teltarg = np.asarray(coord_teltarg.dec.deg, dtype=float)

    final_table = Table()
    final_table["dp_id"] = table["dp_id"]
    final_table["target_name"] = table["target_name"]
    final_table["collection"] = table["obs_collection"]
    final_table["ra_main"] = ra_main
    final_table["dec_main"] = dec_main
    final_table["ra_teltarg"] = ra_teltarg
    final_table["dec_teltarg"] = dec_teltarg
    final_table["ra_main_teltarg_delta"] = wrapped_ra_delta_deg(ra_main, ra_teltarg)
    final_table["dec_main_teltarg_delta"] = dec_main - dec_teltarg

    ra_simbad, dec_simbad = fetch_simbad_coords_parallel(
        final_table["target_name"],
        desc="Resolving SIMBAD coordinates",
        max_workers=max_workers,
    )
    ra_simbad = np.asarray(ra_simbad, dtype=float)
    dec_simbad = np.asarray(dec_simbad, dtype=float)

    final_table["ra_simbad"] = ra_simbad
    final_table["dec_simbad"] = dec_simbad
    final_table["ra_main_simbad_delta"] = wrapped_ra_delta_deg(ra_main, ra_simbad)
    final_table["dec_main_simbad_delta"] = dec_main - dec_simbad
    final_table["ra_teltarg_simbad_delta"] = wrapped_ra_delta_deg(ra_teltarg, ra_simbad)
    final_table["dec_teltarg_simbad_delta"] = dec_teltarg - dec_simbad

    for col in ["dp_id", "collection", "target_name"]:
        force_fixed_str_column(final_table, col)

    final_table.write(output_file, overwrite=overwrite)

    return final_table, N


def main(collection_name, maxrec=1, max_workers=8, overwrite=False):
    return build_main_teltarg_table(
        collection_name,
        maxrec=maxrec,
        max_workers=max_workers,
        overwrite=overwrite,
    )

def force_fixed_str_column(tbl, colname):
    """Convert an object/string column to a fixed-width string dtype."""
    vals = ["" if v is None else str(v) for v in tbl[colname]]
    maxlen = max((len(v) for v in vals), default=1)
    tbl[colname] = np.array(vals, dtype=f"U{maxlen}")  # fixed-width unicode
    return maxlen


def fetch_headers_parallel(eso, dp_ids, desc, max_workers=8):
    dp_ids = list(dp_ids)
    if not dp_ids:
        return Table()

    def _get_header(dp_id):
        return eso.get_headers([dp_id])

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        headers = list(
            tqdm(
                executor.map(_get_header, dp_ids),
                total=len(dp_ids),
                desc=desc,
            )
        )
    return vstack(headers)


def fetch_header_values_parallel(
    eso,
    dp_ids,
    desc,
    max_workers=8,
    direct_keys=None,
    pattern_keys=None,
):
    direct_keys = direct_keys or []
    pattern_keys = pattern_keys or []
    out_keys = [name for name, _ in direct_keys + pattern_keys]
    dp_ids = list(dp_ids)
    if not dp_ids:
        return {name: [] for name in out_keys}

    def _extract_values(header):
        values = {}
        if header is None or len(header) == 0:
            for name in out_keys:
                values[name] = np.nan
            return values

        row = header[0]
        colnames = header.colnames

        for out_name, key in direct_keys:
            value = row[key] if key in colnames else np.nan
            values[out_name] = np.nan if _is_missing_value(value) else value

        for out_name, pattern in pattern_keys:
            matches = sorted(k for k in colnames if fnmatch.fnmatch(k, pattern))
            values[out_name] = coalesce_row_values(row, matches)

        return values

    def _get_values(dp_id):
        if _is_missing_value(dp_id):
            return {name: np.nan for name in out_keys}
        header = eso.get_headers([dp_id])
        return _extract_values(header)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(
            tqdm(
                executor.map(_get_values, dp_ids),
                total=len(dp_ids),
                desc=desc,
            )
        )

    out = {name: [] for name in out_keys}
    for row in results:
        for name in out_keys:
            out[name].append(row.get(name, np.nan))
    return out


def fetch_simbad_coords_parallel(target_names, desc, max_workers=8):
    target_names = list(target_names)
    if not target_names:
        return [], []

    cleaned_names = [
        "" if _is_missing_value(name) else str(name).strip()
        for name in target_names
    ]
    unique_names = list(dict.fromkeys(name for name in cleaned_names if name))

    def _resolve_target(name):
        try:
            simbad = Simbad()
            simbad.add_votable_fields("ra(d)", "dec(d)")
            result = simbad.query_object(name)
            return extract_simbad_j2000_coords(result)
        except Exception:
            return np.nan, np.nan

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(
            tqdm(
                executor.map(_resolve_target, unique_names),
                total=len(unique_names),
                desc=desc,
            )
        )

    resolved_coords = dict(zip(unique_names, results))
    ra_vals = []
    dec_vals = []
    for name in cleaned_names:
        ra_deg, dec_deg = resolved_coords.get(name, (np.nan, np.nan))
        ra_vals.append(ra_deg)
        dec_vals.append(dec_deg)

    return ra_vals, dec_vals


def extract_simbad_j2000_coords(result_table):
    if result_table is None or len(result_table) == 0:
        return np.nan, np.nan

    row = result_table[0]
    colnames = {name.lower(): name for name in result_table.colnames}

    for ra_key, dec_key in [("ra", "dec"), ("ra_d", "dec_d"), ("ra_deg", "dec_deg")]:
        if ra_key in colnames and dec_key in colnames:
            return simbad_coords_to_j2000(
                row[colnames[ra_key]],
                row[colnames[dec_key]],
            )

    if "ra" in result_table.colnames and "dec" in result_table.colnames:
        return simbad_coords_to_j2000(row["ra"], row["dec"])

    if "RA" in result_table.colnames and "DEC" in result_table.colnames:
        return simbad_coords_to_j2000(row["RA"], row["DEC"])

    return np.nan, np.nan


def simbad_coords_to_j2000(ra_value, dec_value):
    if _is_missing_value(ra_value) or _is_missing_value(dec_value):
        return np.nan, np.nan

    try:
        coord = SkyCoord(
            ra=float(ra_value) * u.deg,
            dec=float(dec_value) * u.deg,
            frame="icrs",
        )
    except (TypeError, ValueError):
        ra_text = value_to_clean_text(ra_value)
        dec_text = value_to_clean_text(dec_value)
        if not ra_text or not dec_text:
            return np.nan, np.nan
        coord = SkyCoord(
            ra_text,
            dec_text,
            unit=(u.hourangle, u.deg),
            frame="icrs",
        )

    coord_j2000 = coord.transform_to(FK5(equinox="J2000"))
    return coord_j2000.ra.deg, coord_j2000.dec.deg


def packed_ra_dec_to_skycoord(ra_packed, dec_packed, equinox="J2000"):
    """
    Convert ESO packed sexagesimal RA/Dec to SkyCoord.

    Parameters
    ----------
    ra_packed : float, str, or array-like
        RA in HHMMSS.SSS format (e.g. 200544.317)
    dec_packed : float, str, or array-like
        Dec in DDMMSS.SSS format (e.g. 40252.812 or -40252.812)
    equinox : str
        Equinox string (default J2000)

    Returns
    -------
    SkyCoord
    """
    ra_packed = np.asarray(ra_packed, dtype=float)
    dec_packed = np.asarray(dec_packed, dtype=float)

    # --- RA ---
    ra_h = (ra_packed // 10000).astype(int)
    ra_m = ((ra_packed % 10000) // 100).astype(int)
    ra_s = ra_packed % 100

    # --- Dec ---
    sign = np.where(dec_packed < 0, -1, 1)
    dec_packed = np.abs(dec_packed)

    dec_d = (dec_packed // 10000).astype(int)
    dec_m = ((dec_packed % 10000) // 100).astype(int)
    dec_s = dec_packed % 100

    ra = (ra_h + ra_m / 60 + ra_s / 3600) * 15.0  # hours → degrees
    dec = sign * (dec_d + dec_m / 60 + dec_s / 3600)

    return SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="fk5", equinox=equinox)


def wrapped_ra_delta_deg(ra_main_deg, ra_teltarg_deg):
    ra_main_deg = np.asarray(ra_main_deg, dtype=float)
    ra_teltarg_deg = np.asarray(ra_teltarg_deg, dtype=float)
    return (ra_main_deg - ra_teltarg_deg + 180.0) % 360.0 - 180.0

def find_header_key(headers, pattern):
    return find_header_keys(headers, pattern)[0]


def find_header_keys(headers, pattern):
    matches = sorted(k for k in headers.keys() if fnmatch.fnmatch(k, pattern))
    if not matches:
        raise KeyError(f"No header matches pattern: {pattern}")
    return matches


def _is_missing_value(value):
    if np.ma.is_masked(value):
        return True
    if value is None:
        return True
    if isinstance(value, (bytes, np.bytes_)):
        value = value.decode(errors="ignore")
    if isinstance(value, str) and value.strip() == "":
        return True
    try:
        return np.isnan(value)
    except Exception:
        return False


def value_to_clean_text(value):
    if _is_missing_value(value):
        return ""
    if isinstance(value, (bytes, np.bytes_)):
        value = value.decode(errors="ignore")
    return str(value).strip()


def coalesce_columns(table, keys, fill_value=np.nan):
    if not keys:
        return np.array([], dtype=float)
    nrows = len(table)
    result = np.empty(nrows, dtype=object)
    for i in range(nrows):
        chosen = None
        for key in keys:
            value = table[key][i]
            if not _is_missing_value(value):
                chosen = value
                break
        result[i] = fill_value if chosen is None else chosen
    return result


def coalesce_row_values(row, keys, fill_value=np.nan):
    for key in keys:
        value = row[key]
        if not _is_missing_value(value):
            return value
    return fill_value

def remove_suffix_from_column(col, suffix=".fits"):
    """
    Remove a trailing suffix from all entries in a column-like object.

    Parameters
    ----------
    col : array-like
        Iterable of values (Astropy Column, pandas Series, list, ndarray).
    suffix : str, optional
        Suffix to remove (default: ".fits").

    Returns
    -------
    list
        List of strings with suffix removed (only if present at end).
    """
    suffix_len = len(suffix)

    values = []
    for v in col:
        if _is_missing_value(v):
            values.append("")
            continue
        text = str(v)
        values.append(text[:-suffix_len] if text.endswith(suffix) else text)
    return values


def plot_pointing_offsets(
    table,
    collection_name,
    p3orig,
    x_col="ra_main_teltarg_delta",
    y_col="dec_main_teltarg_delta",
    outdir="./figures",
    bins=100,
    lim=60,
    figsize=(18, 6),
    panel_specs=None,
):
    """
    Plot RA/Dec pointing offsets with scatter + marginal histograms.

    Parameters
    ----------
    table : astropy.table.Table or pandas.DataFrame
        Input table containing offset columns.
    collection_name : str
        Name used for plot annotation and filename.
    p3orig : str
        Subdirectory name for output.
    x_col, y_col : str
        Column names for the default first RA and Dec offset panels (degrees).
    outdir : str
        Base output directory.
    bins : int
        Number of histogram bins.
    lim : float
        Axis limits (+/- lim).
    figsize : tuple
        Figure size.
    panel_specs : list of tuple, optional
        List of (x_col, y_col, title) panel definitions. Values are assumed to
        be in degrees and converted to arcsec for plotting.
    """

    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    import os

    if panel_specs is None:
        panel_specs = [
            (x_col, y_col, "Main - Tel Targ"),
            ("ra_main_simbad_delta", "dec_main_simbad_delta", "Main - SIMBAD"),
            ("ra_teltarg_simbad_delta", "dec_teltarg_simbad_delta", "Tel Targ - SIMBAD"),
        ]

    fig = plt.figure(figsize=figsize)
    outer_gs = GridSpec(1, len(panel_specs), figure=fig, wspace=0.28)

    for i, (panel_x_col, panel_y_col, panel_title) in enumerate(panel_specs):
        x = np.asarray(table[panel_x_col], dtype=float) * 3600.0
        y = np.asarray(table[panel_y_col], dtype=float) * 3600.0
        valid = np.isfinite(x) & np.isfinite(y)
        x = x[valid]
        y = y[valid]

        inner_gs = outer_gs[i].subgridspec(
            2,
            2,
            height_ratios=[1, 4],
            width_ratios=[4, 1],
            hspace=0.0,
            wspace=0.0,
        )

        ax_histx = fig.add_subplot(inner_gs[0, 0])
        ax_scatter = fig.add_subplot(inner_gs[1, 0], sharex=ax_histx)
        ax_histy = fig.add_subplot(inner_gs[1, 1], sharey=ax_scatter)
        ax_corner = fig.add_subplot(inner_gs[0, 1])
        ax_corner.axis("off")

        ax_scatter.scatter(x, y, s=10, alpha=0.5)
        ax_scatter.axhline(0.0, color="0.5", lw=0.8, zorder=0)
        ax_scatter.axvline(0.0, color="0.5", lw=0.8, zorder=0)
        ax_scatter.set_xlim(-lim, lim)
        ax_scatter.set_ylim(-lim, lim)
        ax_scatter.set_title(panel_title)
        ax_scatter.set_xlabel(f"RA delta ({panel_title}) [arcsec]")
        if i == 0:
            ax_scatter.set_ylabel("Dec delta [arcsec]")
        else:
            ax_scatter.tick_params(labelleft=False)

        ax_histx.hist(x, bins=bins, range=(-lim, lim), ec="black")
        ax_histy.hist(
            y,
            bins=bins,
            range=(-lim, lim),
            orientation="horizontal",
            ec="black",
            fc="C1",
        )

        ax_histx.tick_params(labelbottom=False)
        ax_histy.tick_params(labelleft=False)
        ax_histx.grid(False)
        ax_histy.grid(False)

    fig.suptitle(collection_name, fontsize=16, y=0.98)
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.12, top=0.88, wspace=0.28)

    # Output
    os.makedirs(f"{outdir}/{p3orig}", exist_ok=True)
    outfile = f"{outdir}/{p3orig}/{collection_name}_ra_dec_deltas.png"
    fig.savefig(outfile, dpi=300)

    return fig, outfile
