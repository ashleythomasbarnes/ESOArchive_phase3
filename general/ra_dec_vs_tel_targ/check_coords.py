from astropy.coordinates import SkyCoord
import astropy.units as u


def packed_ra_dec_to_skycoord(ra_packed, dec_packed, equinox="J2000"):
    """
    Convert ESO packed sexagesimal RA/Dec to SkyCoord.

    Parameters
    ----------
    ra_packed : float or str
        RA in HHMMSS.SSS format (e.g. 200544.317)
    dec_packed : float or str
        Dec in DDMMSS.SSS format (e.g. 40252.812 or -40252.812)
    equinox : str
        Equinox string (default J2000)

    Returns
    -------
    SkyCoord
    """
    ra_packed = float(ra_packed)
    dec_packed = float(dec_packed)

    # --- RA ---
    ra_h = int(ra_packed // 10000)
    ra_m = int((ra_packed % 10000) // 100)
    ra_s = ra_packed % 100

    # --- Dec ---
    sign = -1 if dec_packed < 0 else 1
    dec_packed = abs(dec_packed)

    dec_d = int(dec_packed // 10000)
    dec_m = int((dec_packed % 10000) // 100)
    dec_s = dec_packed % 100

    ra = (ra_h + ra_m / 60 + ra_s / 3600) * 15.0  # hours → degrees
    dec = sign * (dec_d + dec_m / 60 + dec_s / 3600)

    return SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="fk5", equinox=equinox)


# Example using your header values
c = packed_ra_dec_to_skycoord("200544.317", "40252.812")

print("RA  [deg]:", c.ra.deg)
print("Dec [deg]:", c.dec.deg)