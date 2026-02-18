"""Generate absorbed surface flux files for the heat1d thermal model.

Computes the absorbed solar flux at a given latitude using the standard
heat1d insolation physics (angle-dependent albedo, orbital geometry)
and writes it in the standard flux file format readable by both the
Python Model and C backend.

Supports optional eclipse/shadow events with configurable duration,
start time, and obscuration fraction.

Epoch convention (consistent with heat1d):
    Local time 0 = local noon (hour angle = 0).
    flux[0] corresponds to model time t = 0.

Note on time units:
    Local time arguments (--t-start, --t-stop, --eclipse-start) are in
    *planetary* hours on a 0-24 scale, where 24 = one full diurnal cycle.
    On the Moon, 1 local hour = planet.day/24 ~ 29.5 Earth hours.
    Eclipse duration (--eclipse-duration) is in SI seconds.

Usage::

    python -m heat1d.generate_flux --help
    generate-flux --lat 0 -o flux_equator.txt
    generate-flux --lat 0 --eclipse-start 2.0 --eclipse-duration 3600 -o flux.txt
"""

import argparse
import sys

import numpy as np
import planets

from . import orbits
from .flux import write_flux_file
from .properties import albedoVar


def compute_flux_array(planet, lat_rad, nsteps, t_start_hr=0.0,
                       t_stop_hr=24.0, dec_rad=0.0, r=None,
                       lon=0.0):
    """Compute absorbed solar flux over a local-time window.

    Parameters
    ----------
    planet : object
        Planet object with S, albedo, albedoCoef, day, rAU, year,
        eccentricity attributes.
    lat_rad : float
        Latitude [rad].
    nsteps : int
        Total number of time steps in the output array.
    t_start_hr : float
        Start local time [decimal hours, 0 = noon].
    t_stop_hr : float
        Stop local time [decimal hours].
    dec_rad : float
        Solar declination [rad]. Used only for circular orbit (ecc <= 0.01).
    r : float or None
        Heliocentric distance [AU]. Used only for circular orbit (ecc <= 0.01).
    lon : float
        Observer longitude [rad]. Default 0.

    Returns
    -------
    flux : np.ndarray
        Absorbed surface flux [W/m^2], length *nsteps*.
    dt : float
        Time step [s].
    """
    duration_s = (t_stop_hr - t_start_hr) / 24.0 * planet.day
    if duration_s <= 0:
        raise ValueError("t_stop must be greater than t_start")
    dt = duration_s / nsteps

    t_start_s = t_start_hr * planet.day / 24.0
    t = t_start_s + np.arange(nsteps) * dt

    Sabs = planet.S * (1.0 - planet.albedo)
    a_coef = planet.albedoCoef[0]
    b_coef = planet.albedoCoef[1]

    ecc = planet.eccentricity
    if ecc > 0.05:
        # General eccentric orbit
        P_sid = orbits.siderealPeriod(planet.day, planet.year)
        M0 = orbits.meanFromTrue(lon, ecc)
        M = M0 + orbits.TWOPI * t / planet.year
        nu = orbits.trueFromMean(M, ecc)

        h = orbits.TWOPI * t / P_sid + lon - nu

        a_au = planet.rAU
        x = a_au * (1 - ecc**2)
        r_t = x / (1 + ecc * np.cos(nu))

        obliq = planet.obliquity
        Lp = planet.Lp if planet.Lp is not None else 0.0
        dec_t = np.arcsin(np.sin(obliq) * np.sin(nu + Lp))

        c = orbits.cosSolarZenith(lat_rad, dec_t, h)
        inc = np.arccos(c)
        A_var = albedoVar(planet.albedo, a_coef, b_coef, inc)
        f = (1.0 - A_var) / (1.0 - planet.albedo)
        flux = f * Sabs * (r_t / planet.rAU) ** -2 * c
    else:
        # Circular orbit approximation
        if r is None:
            r = planet.rAU
        h = orbits.hourAngle(t, planet.day)
        c = orbits.cosSolarZenith(lat_rad, dec_rad, h)
        inc = np.arccos(c)
        A_var = albedoVar(planet.albedo, a_coef, b_coef, inc)
        f = (1.0 - A_var) / (1.0 - planet.albedo)
        flux = f * Sabs * (r / planet.rAU) ** -2 * c

    return flux, dt


def apply_eclipse(flux, dt, planet, t_start_hr, eclipse_start_hr,
                  eclipse_duration_s, eclipse_fraction=1.0):
    """Zero out or reduce flux during an eclipse window.

    Parameters
    ----------
    flux : np.ndarray
        Absorbed flux array (modified in-place).
    dt : float
        Time step [s].
    planet : object
        Planet object (for day length).
    t_start_hr : float
        Local time corresponding to flux[0] [decimal hours].
    eclipse_start_hr : float
        Eclipse start in local time [decimal hours, 0 = noon].
    eclipse_duration_s : float
        Eclipse duration [SI seconds].
    eclipse_fraction : float
        Fraction of flux blocked (1.0 = total shadow, 0.5 = 50%).

    Returns
    -------
    n_affected : int
        Number of flux samples modified.
    n_daytime : int
        Number of those that had non-zero flux (daytime samples).
    """
    t_start_s = t_start_hr * planet.day / 24.0
    t = t_start_s + np.arange(len(flux)) * dt

    eclipse_t0 = eclipse_start_hr * planet.day / 24.0
    eclipse_t1 = eclipse_t0 + eclipse_duration_s

    mask = (t >= eclipse_t0) & (t < eclipse_t1)
    n_affected = int(np.sum(mask))
    n_daytime = int(np.sum(flux[mask] > 0)) if n_affected > 0 else 0
    flux[mask] *= (1.0 - eclipse_fraction)
    return n_affected, n_daytime


def _required_nsteps(planet, t_start_hr, t_stop_hr, eclipse_duration_s,
                     min_samples_per_eclipse=10):
    """Compute minimum total nsteps to resolve an eclipse.

    Ensures at least *min_samples_per_eclipse* flux samples fall within
    the eclipse window.
    """
    duration_s = (t_stop_hr - t_start_hr) / 24.0 * planet.day
    dt_max = eclipse_duration_s / min_samples_per_eclipse
    return int(np.ceil(duration_s / dt_max))


def plot_flux(flux, dt, planet, lat_deg, t_start_hr=0.0,
              eclipse_info=None, save_path=None):
    """Plot the flux time series.

    Parameters
    ----------
    flux : np.ndarray
        Absorbed flux [W/m^2].
    dt : float
        Time step [s].
    planet : object
        Planet object.
    lat_deg : float
        Latitude [degrees] (for title).
    t_start_hr : float
        Local time of flux[0] [decimal hours].
    eclipse_info : dict or None
        Keys: 'start_hr', 'duration_s', 'fraction'.
    save_path : str or None
        If given, save plot to this file instead of showing.
    """
    import matplotlib.pyplot as plt

    t_hr = t_start_hr + np.arange(len(flux)) * dt / planet.day * 24.0

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t_hr, flux, 'k-', linewidth=1)
    ax.set_xlabel("Local time [hours past noon]")
    ax.set_ylabel("Absorbed flux [W/m$^2$]")
    ax.set_title(f"Surface flux — lat={lat_deg:.1f}°, {planet.name}")
    ax.set_xlim(t_hr[0], t_hr[-1])
    ax.set_ylim(bottom=0)

    if eclipse_info is not None:
        e_start = eclipse_info['start_hr']
        e_end = e_start + eclipse_info['duration_s'] / planet.day * 24.0
        ax.axvspan(e_start, e_end, alpha=0.25, color='gray',
                   label=f"Eclipse ({eclipse_info['fraction']*100:.0f}%)")
        ax.legend()

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"Plot saved: {save_path}")
    else:
        plt.show()
    plt.close(fig)


def main(argv=None):
    """CLI entry point for generating flux files."""
    parser = argparse.ArgumentParser(
        description="Generate absorbed surface flux files for heat1d.",
        epilog=(
            "Local time convention: 0 = local noon (hour angle = 0). "
            "One local 'hour' = planet.day/24 in SI seconds "
            "(~29.5 Earth hours on the Moon)."
        ),
    )
    parser.add_argument("--lat", type=float, default=0.0,
                        help="Latitude [degrees] (default: 0)")
    parser.add_argument("--t-start", type=float, default=0.0,
                        help="Start local time [planetary hours, 0=noon] (default: 0)")
    parser.add_argument("--t-stop", type=float, default=24.0,
                        help="Stop local time [planetary hours] (default: 24)")
    parser.add_argument("--nsteps", type=int, default=None,
                        help="Total number of flux samples (default: auto, "
                             "at least 480/day; increased to resolve eclipses)")
    parser.add_argument("--planet", type=str, default="Moon",
                        help="Planet name from planets package (default: Moon)")
    parser.add_argument("--albedo", type=float, default=None,
                        help="Override normal-incidence albedo A0")
    parser.add_argument("--declination", type=float, default=0.0,
                        help="Solar declination [degrees] (default: 0)")
    parser.add_argument("--distance", type=float, default=None,
                        help="Heliocentric distance [AU] (default: planet.rAU)")

    eclipse = parser.add_argument_group("eclipse/shadow options")
    eclipse.add_argument("--eclipse-start", type=float, default=None,
                         help="Eclipse start [planetary local time hours, 0=noon]")
    eclipse.add_argument("--eclipse-duration", type=float, default=None,
                         help="Eclipse duration [SI seconds]")
    eclipse.add_argument("--eclipse-fraction", type=float, default=1.0,
                         help="Fraction of flux blocked, 0-1 (default: 1.0 = total)")

    output = parser.add_argument_group("output options")
    output.add_argument("-o", "--output", type=str, default="flux_output.txt",
                        help="Output file path (default: flux_output.txt)")
    output.add_argument("--plot", action="store_true",
                        help="Show a plot of the generated flux")
    output.add_argument("--plot-file", type=str, default=None,
                        help="Save plot to file (implies --plot)")

    args = parser.parse_args(argv)

    # Resolve planet
    try:
        planet = getattr(planets, args.planet)
    except AttributeError:
        parser.error(f"Unknown planet: {args.planet}")

    # Apply albedo override
    if args.albedo is not None:
        planet.albedo = args.albedo

    lat_rad = np.deg2rad(args.lat)
    dec_rad = np.deg2rad(args.declination)

    # Determine nsteps: default 480/day, auto-increase to resolve eclipses
    ndays = (args.t_stop - args.t_start) / 24.0
    default_nsteps = int(round(480 * ndays))

    if args.nsteps is not None:
        nsteps = args.nsteps
    elif args.eclipse_duration is not None:
        nsteps_eclipse = _required_nsteps(
            planet, args.t_start, args.t_stop, args.eclipse_duration)
        nsteps = max(default_nsteps, nsteps_eclipse)
    else:
        nsteps = default_nsteps

    # Compute flux
    flux, dt = compute_flux_array(
        planet, lat_rad, nsteps,
        t_start_hr=args.t_start, t_stop_hr=args.t_stop,
        dec_rad=dec_rad, r=args.distance,
    )

    # Apply eclipse if specified
    eclipse_info = None
    if args.eclipse_start is not None and args.eclipse_duration is not None:
        n_affected, n_daytime = apply_eclipse(
            flux, dt, planet, args.t_start,
            args.eclipse_start, args.eclipse_duration,
            args.eclipse_fraction)
        eclipse_info = {
            'start_hr': args.eclipse_start,
            'duration_s': args.eclipse_duration,
            'fraction': args.eclipse_fraction,
        }
    elif args.eclipse_start is not None or args.eclipse_duration is not None:
        parser.error("Both --eclipse-start and --eclipse-duration are required")

    # Write flux file
    write_flux_file(args.output, flux, dt)

    # Summary
    print(f"Wrote {len(flux)} flux samples (dt={dt:.2f} s) to {args.output}")
    print(f"  Planet:      {planet.name}")
    print(f"  Latitude:    {args.lat:.2f} deg")
    print(f"  Local time:  {args.t_start:.2f} - {args.t_stop:.2f} planetary hr")
    print(f"  Declination: {args.declination:.2f} deg")
    print(f"  1 local hr = {planet.day/24.0:.1f} s = {planet.day/24.0/3600:.2f} Earth hr")
    if args.albedo is not None:
        print(f"  Albedo (A0): {args.albedo:.4f}")
    if eclipse_info:
        ecl_dur_lt = args.eclipse_duration / planet.day * 24.0
        print(f"  Eclipse:     start={args.eclipse_start:.2f} local hr, "
              f"duration={args.eclipse_duration:.0f} s "
              f"({ecl_dur_lt:.4f} local hr), "
              f"fraction={args.eclipse_fraction:.2f}")
        print(f"  Eclipse:     {n_affected} samples affected, "
              f"{n_daytime} during daytime")
        if n_daytime == 0:
            print(f"  WARNING: eclipse falls entirely during nighttime "
                  f"(flux already zero) -- no effect on temperatures")
    print(f"  Flux max:    {flux.max():.2f} W/m^2")
    print(f"  Flux mean:   {flux.mean():.2f} W/m^2")

    # Plot
    if args.plot or args.plot_file:
        plot_flux(flux, dt, planet, args.lat, t_start_hr=args.t_start,
                  eclipse_info=eclipse_info, save_path=args.plot_file)


if __name__ == "__main__":
    main()
