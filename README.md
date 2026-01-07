# Readme

This library provides interaction with multiple subjects relating to aerodynamics and atmospherics.

## Atmospherics

This crate seeks, in part, to provide access to some of the standard atmospheres used in aerospace and atmospheric science. The following atmospheres are planned:

- International Standard Atmosphere (ISA) (ICAO Doc 7488/3, 1993)
   - Work in progress. See [ISA.md](ISA.md) for details.
- US Standard Atmosphere, 1976
    - Not yet started
- ISO Standard Atmosphere (ISO 2533, 1975)
    - Not yet started

At altitudes near the earth's surface, these models provide similar data. However, the altitude ranges, some constants, and other details differ. Difrent applications require different atmospheres, so it is useful to have access to multiple models.

## Testing

Anywhere where it is practical, calculations are tested to a precision of 0.0005, that is 0.05%. The constant `PRECISION` is devined in `lib.rs`.

## Long-Term Goals

**Aerodynamics** — In addition to atmospherics, this crate will eventually provide functions and data structures useful for aerodynamic calculations. This may include airfoil analysis, including calculating $C_L$ and $C_D$, and other related functionality. The Open-source project [XFOIL](http://web.mit.edu/drela/Public/web/xfoil/), Distributed under the GNU General Public License, provides Fortran code for airfoil analysis that may be useful as a reference. This crate may eventually provide Rust implementations of some of XFOIL's functionality.

**GUI App** — This crate, in and of itself, is useful for incorporation into larger codebases. However, many end-users may not wish to interact with the crate at this level. Therefore, long-term goals include building a stand-along GUI application. This application will likely be written Rust using the [Dioxus](https://dioxuslabs.com/) GUI framework.
