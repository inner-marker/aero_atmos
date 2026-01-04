# International Standard Atmosphere

## Implimentations

Geopotential altitude (H) __to__ ...

- [x] Temperature ($T$)
- [ ] Pressure ($P$)
- [ ] Pressure Ratio ($p/P_0$)
- [ ] Density ($\rho$, rho)
- [ ] Density Ratio ($ \rho/\rho_0 $)
- [ ] Square Root of Density Ratio ($ \sqrt{\rho/\rho_0} $)
- [ ] Speed of Sound ($a$)
- [ ] Dynamic Viscosity (μ, mu)
- [ ] Kinematic Viscosity ($ν$) (ν = μ/ρ)
- [ ] Thermal Conductivity ($λ$, lambda)
- [ ] Spcific Weights ($γ$, gamma) (γ = g)
- [ ] Number Density ($n$) (N = ρ/M)
- [ ] Mean particle speed ($\bar{v}$)
- [ ] Collision Frequency ($\omega$)
- [ ] Mean free path ($ l $)

Generally, __to__ functions will return a `Result<uom::si::f64::XXX>`, where `XXX` is the appropriate quantity.

Geopotential altitude (H) __from__ ...

- [ ] Temperature ($T$)
- [ ] Pressure ($P$)
- [ ] Density ($\rho$, rho)

For the __from__ functions, a `Result<Vec<uom::si::f64::Length>>` will be returned.

## Symbols and Units

| Symbol        | Name                               | Unit             | Value (if applicable)           | Notes or Derrivation             | IMpl'd? |
|---------------|------------------------------------|------------------|---------------------------------|----------------------------------|---------|
| $ a $         | Speed of sound                     | m/s              |                                 |                                  |         |
| $ h $         | Geometric altitude                 | m                |                                 |                                  |         |
| $ g $         | acceleration due to gravity        | m/s²             |                                 |                                  |         |
| $ G $         | Standard gravity at sea level      | m/s²             | 9.806_65                        |                                  | [x]     |
| $ H $         | Geopotential altitude              | m                |                                 |                                  |         |
| $ H_b $       | Base geopotential altitude         | m                | Varies by layer, 0 @ SL         |                                  |         |
| $ l $         | Mean free path                     |                  |                                 |                                  |         |
| $ M_0 $       | Molar mass of Earth's air          | g/mol            | 28_964.420                      |                                  |         |
| $ n $         | Number density                     |                  |                                 | $ N = \frac{\rho}{M} $           |         |
| $ N_A $       | Avogadro's number                  | $kmol^{-1}$      | 602.257e24                      |                                  |         |
| $ p $         | Pressure                           | Pa               |                                 |                                  |         |
| $ P_0 $       | pressure @SL                       | Pa               | 101_325 @ SL                    |                                  |         |
| $ P_b $       | Base pressure at layer base        | Pa               | Varies by layer, 101_325 @ SL   |                                  |         |
| $ R_e $       | Earth's radius                     | m                | 6_356_766                       |                                  |         |
| $ R $         | Specific gas constant for air      | J/(kg·K)         | 287.052_87                      | $ R = \frac{R^*}{M_0} $          |         |
| $ R* $        | Universal gas constant             | J/(kmol·K)       | 831.432                         |                                  |         |
| $ S $         | Sutherland constant                | K                | 110.4                           |                                  |         |
| $ T $         | Temperature                        | K                |                                 |                                  |         |
| $ T_0 $       | Standard temperature @SL           | K                | 288.15                          |                                  |         |
| $ t_0 $       | Standard temperature @SL           | °C               | 15                              | $ t_0 = T_0 - 273.15 $          |         |
| $ T_i $       | Kelvin ice point temperature @SL   | K                | 273.15                          |                                  |         |
| $ t_i $       | Celcius ice point temperature @SL  | °C               | 0                               | $ t_i = T_i - 273.15 $          |         |
| $ T_b $       | Base temperature at layer base     | K                | Varies by layer, 288.15 @ SL    |                                  |         |
| $ \bar{v} $   | Mean particle speed                | m/s              |                                 |                                  |         |
| $ \rho $ rho  | Air density                        | kg/m³            |                                 |                                  |         |
| $ \rho_0 $    | Standard air density @SL           | kg/m³            | 1.225                           |                                  |         |
| $ \beta $ beta | Temperature lapse rate             | K/km             | Varies by layer, -6.5 @ SL      |                                  |         |
| $ \beta_S $   |                                    | $kg/m•s•K^{1/2}$ | 1.458e-6                        |                                  |         |
| $ \gamma $ gamma | Specific weight                    | N/m³             |                                 | $ \gamma = \rho g $              |         |
| $ \kappa $ kappa | Ratio of specific heats            |                  | 1.4                             | Dimensionelss                    |         |
| $ \lambda $ lambda | Thermal conductivity               | W/(m·K)          |                                 |                                  |         |
| $ \mu $ mu    | Dynamic viscosity                  | Pa·s             |                                 |                                  |         |
| $ \nu $ nu    | Kinematic viscosity                | m²/s             |                                 |                                  |         |
| $ \sigma $    | Effective collision diameter       | m                | 0.365e-9                        |                                  |         |
| $ \omega $    | Collision frequency                | $ s^{-1} $       |                                 |                                  |         |
|               |                                    |                  |                                 |                                  |         |
|               |                                    |                  |                                 |                                  |         |
|               |                                    |                  |                                 |                                  |         |
|               |                                    |                  |                                 |                                  |         |

[Greek letters.](GREEK.md)

## Pressure

See [PRESSURE.md](PRESSURE.md) for the formulas and constants used in pressure calculations.
