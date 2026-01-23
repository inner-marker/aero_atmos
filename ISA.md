# International Standard Atmosphere

## Implimentations

### Geopotential altitude ($H$) __to__

- [x] Temperature ($T$)
- [x] Temperature Ratio ($T/T_0$)
- [x] Pressure ($P$)
- [x] Pressure Ratio ($p/P_0$)
- [x] Density ($\rho$, rho)
- [x] Density Ratio ($\rho/\rho_0$)
- [x] Square Root of Density Ratio ($\sqrt{\rho/\rho_0}$)
- [x] Speed of Sound ($a$)
- [x] Dynamic Viscosity ($\mu$, mu)
- [x] Kinematic Viscosity ($ν$) ($ν = μ/ρ$)
- [x] Thermal Conductivity ($λ$, lambda)
- [x] Number Density ($n$) ($N = \rho/M$)
- [x] Mean particle speed ($\bar{v}$)
- [x] Mean free path ($l$)
- [x] Collision Frequency ($\omega$)
- [x] Gravity acceleration, local $g$
- [x] Specific Weight ($\gamma$, gamma) $\gamma = \rho g$

Generally, __to__ functions will return a `Result<uom::si::f64::XXX>`, where `XXX` is the appropriate quantity.

### Geopotential altitude ($H$) __from__

For the __from__ functions, where appropriate, a `Result<Vec<uom::si::f64::Length>, IsaError>` will be returned. For situations where multiple altitudes correspond to the input value (for example, temperature in the stratosphere), all possible geopotential altitudes will be returned. For all other situations, the function returns a `Result<uom::si::f64::Length, IsaError>`.

- [x] Pressure ($P$)
- [ ] Pressure Ratio ($P/P_0$)
- [x] Density ($\rho$, rho)
- [ ] Density Ratio ($\rho/\rho_0$)

For most item marked complete, there are extensive tests in the `tests` module to verify the calculations to a precision of 0.0005 (0.05%).
For some smaller functions, such as square root of density ratio, only a `doctest` may be provided.

### Constants

Functions are provided which return either an `f64` or a `uom::si::f64::XXX` for the following constants:

- [x] Earth's Radius ($R_e$)
- [x] Gravitational acceleration at sea level ($G$)
- [x] Molar mass of Earth's air ($M_0$)
- [x] Avogadro's number ($N_A$)
- [x] Specific gas constant for air ($R$)
- [x] Universal gas constant ($R*$)
- [x] Sutherland's Constant ($β_S$)
- [x] Sutherland constant ($S$)
- [x] Ratio of specific heats (κ, kappa)
- [x] Effective collision diameter ($σ$)
- [x] Standard pressure at sea level ($P_0$)

## Layers

According to ICAO Doc 7488/3, the atmosphere is divided into layers. Each layer has a base geopotential altitude ($H_b$), base temperature ($T_b$), base pressure ($P_b$), and temperature lapse rate ($β$). The layers are as follows:

| Layer Name   | Geopotential Altitude Range         |
|--------------|-------------------------------------|
|              | -5 km to 0 km                       |
| Troposphere  | 0 km to 11 km                       |
| Tropopause   | 11 km to 20 km                      |
| Stratosphere | 20 km to 32 km                      |
| Stratopause  | 32 km to 47 km                      |
| Mesosphere   | 47 km to 51 km                      |
| Mesopause    | 51 km to 71 km                      |
| Thermosphere | 71 km to 84.852 km                  |

## Symbols and Units

THis table represents the symbols, units, and values (if applicable) used in the International Standard Atmosphere (ISA) model, ICAO Doc 7488/3.

| Symbol        | Name                               | Unit             | Value (if applicable)           | Notes or Derrivation     |
|---------------|------------------------------------|------------------|---------------------------------|--------------------------|
|  $a$          | Speed of sound                     | $m/s$              |                                 |                          |
|  $g$          | acceleration due to gravity        | $m/s^2$             |                                 |                          |
| $G$           | Standard gravity at sea level      | $m/s^2$             | 9.806_65                        |                          |
|  $h$          | Geometric altitude                 | $m$                |                                 |                          |
| $H$           | Geopotential altitude              | $m$                |                                 |                          |
| $H_b$         | Base geopotential altitude         | $m$                | Varies by layer, 0 @ SL         |                          |
| $l$           | Mean free path                     |                  |                                 |                          |
| $M_0$         | Molar mass of Earth's air          | $g/mol$            | 28_964.420                      |                          |
| $n$         | Number density                     |                  |                                 | $N = \frac{\rho}{M}$   |
| $N_A$       | Avogadro's number                  | $kmol^{-1}$      | 602.257e24                      |                          |
| $p$         | Pressure                           | $Pa$               |                                 |                          |
| $P_0$       | pressure @SL                       | $Pa$               | 101_325 @ SL                    |                          |
| $P_b$       | Base pressure at layer base        | $Pa$               | Varies by layer, 101_325 @ SL   |                          |
| $R_e$       | Earth's radius                     | $m$                | 6_356_766                       |                          |
| $R$         | Specific gas constant for air      | $J/(kg·K)$         | 287.052_87                      | $R = \frac{R^*}{M_0}$  |
| $R^*$        | Universal gas constant             | $J/(kmol·K)$       | 831.432                         |                          |
| $S$         | Sutherland constant                | $K$                | 110.4                           |                          |
| $T$         | Temperature                        | $K$                |                                 |                          |
| $T_0$       | Standard temperature @SL           | $K$                | 288.15                          |                          |
| $t_0$       | Standard temperature @SL           | $°C$               | 15                              | $t_0 = T_0 - 273.15$  |
| $T_i$       | Kelvin ice point temperature @SL   | $K$                | 273.15                          |                          |
| $t_i$       | Celcius ice point temperature @SL  | $°C$               | 0                               | $t_i = T_i - 273.15$  |
| $T_b$       | Base temperature at layer base     | $K$                | Varies by layer    |                          |
| $\bar{v}$  | Mean particle speed                | $m/s$              |                                 |                          |
| $\rho$  | Air density                        | $kg/m³$            |                                 |                          |
| $\rho_0$    | Standard air density @SL           | $kg/m³$            | 1.225                           |                          |
| $\beta$ | Temperature lapse rate             | $K/km$             | Varies by layer      |                          |
| $\beta_S$   | Sutherland's Constant               | $kg/m•s•K^{1/2}$ | 1.458e-6                        |                          |
| $\gamma$ | Specific weight                    | $N/m³$             |                                 | $\gamma = \rho g$      |
| $\kappa$ | Ratio of specific heats            |                  | 1.4                             | Dimensionels            |
| $\lambda$ | Thermal conductivity               | $W/(m·K)$          |                                 |                          |
| $\mu$    | Dynamic viscosity                  | $Pa·s$             |                                 |                          |
| $\nu$    | Kinematic viscosity                | $m²/s$             |                                 |                          |
| $\sigma$    | Effective collision diameter       | $m$                | 0.365e-9                        |                          |
| $\omega$   | Collision frequency                | $s^{-1}$       |                                 |                          |

[Greek letters.](GREEK.md)

## Pressure

See [PRESSURE.md](PRESSURE.md) for the formulas and constants used in pressure calculations.

