# Readme

This library provides interaction with multiple subjects relating to aerodynamics and atmospherics.

## atmospherics

The primary atmospheric model here is the standard atmosphere. Two standard atmospheres are provided:

- 1976 US Standard Atmosphere
- 1993 International Standard Atmosphere (ISA) (ICAO Doc 7488/3)

## Implimentations

See [ISA.md](ISA.md) for constants and symbols used, and for the progress of implimenting computations.

## Testing

Anywhere where it is practical, calculations are tested to a precision of 0.0005, that is 0.05%. The constant `PRECISION` is devined in `lib.rs`.
