# Pressure

Here are the formulas for pressure with height.

The relationship between geopotential altitude and pressure is derived from the hydrostatic equation and the ideal gas law, using geopotential altitude $ H $ as the vertical coordinate.
 The hydrostatic equation in terms of geopotential altitude is:

$$ \frac{dP}{dH} = - \rho g $$

where $ P $ is pressure, $ \rho $ is air density, and $ G $ is the acceleration due to gravity at sea level.
 Using the ideal gas law, $ \rho = \frac{MP}{RT} $, where $ M $ is the molar mass of air, $ R $ is the universal gas constant, and $ T $ is absolute temperature, the equation becomes:

$$ \frac{dP}{dH} = - \frac{MPg}{RT} $$

For a given atmospheric layer where temperature $ T $ varies linearly with geopotential altitude $ H $, such as $ T = T_b + L(H - H_b) $, where $ L $ is the temperature lapse rate, the equation can be integrated to yield pressure as a function of geopotential altitude.
 For a layer with constant temperature ($ L = 0 $), the solution is:

$$ P = P_b \exp\left(-\frac{Mg(H - H_b)}{RT_b}\right) $$

For a layer with a constant lapse rate ($ L \neq 0 $), the solution is:

$$ P = P_b \left( \frac{T_b}{T_b + L(H - H_b)} \right)^{\frac{Mg}{RL}} $$

In effect, removing $ P_b $ from the left side of the equation above gives the pressure ratio. This, so long as $ P $ and $ P_b $ are in the same units, the equation holds for any pressure units.

These equations are used in the standard atmosphere model, where temperature profiles are defined in discrete layers based on geopotential altitude.
 The standard atmosphere is defined in terms of geopotential altitude to simplify calculations by using a constant $ G $ instead of the variable gravitational acceleration.
 The geopotential altitude $ H $ is related to geometric altitude $ Z $ by the equation:

$$ H = \frac{R_e Z}{R_e + Z} $$

where $ E $ is the Earth's radius.
 This relationship accounts for the variation of gravity with altitude and ensures that a small change in geopotential altitude produces the same change in gravitational potential energy as a change in geometric altitude at sea level.

## Constants used here

These are the specific numerical values for the constants used in the pressure equations above. These are not necessarilty the same as those used in the ICAO Doc 7488/3.

|Symbol|Name|Value|Unit|
|------|----|----|-----|
|$$ M $$|Molar mass of Earth's air|28_964.42|g/mol|
|$$ R $$|Universal gas constant|831.325|J/(kmol·K)|
|$$ g $$|Standard gravity at sea level|9.806_65|m/s²|
|$$ R_e $$|Earth's radius|6_356_766|m|
|$$ P_b $$|Base pressure at layer base|Varies by layer, 101_325 @ SL|Pa|
|$$ L $$|Temperature lapse rate|Varies by layer, -6.5 @ SL|K/km|
|$$ T_b $$|Base temperature at layer base|Varies by layer, 288.15 @ SL|K|
|$$ H_b $$|Base geopotential altitude at layer base|Varies by layer, 0 @ SL|m|