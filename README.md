# Droplet Evaporation 0D model
## How to install

## How to use it

## Technical backgroud

# Droplet Evaporation Simulation Model

The 0-D model used for the simulation of a single droplet evaporation was developed according to a diffusion model based on mass and energy balance well reported in literature [35,36], following the classical evaporation model (CEM) proposed by Spalding [37] and Godsave [38]. This model was further improved following Abramzon-Sirignano approach [39,40], enabling a better account for advective mass and energy transport. A detailed description of the governing equations and limitations of this model was recently reviewed by Pinheiro et al. [41]. Mass and energy balances can be described according to equations (4.1) and (4.2), respectively:

\[\frac{dm_d}{dt} = -\dot{m}_d\]  (4.1)

\[m_d c_{pl} \frac{dT_d}{dt} = -Q_s\]  (4.2)

where \(m_d\) is the droplet mass, \(\dot{m}_d\) is the mass evaporation rate, \(c_{pl}\) is the specific heat capacity of the liquid droplet, \(T_d\) is the droplet temperature, and \(Q_s\) is the net power transferred from the environment contributing to droplet temperature increase.

Equation (4.3) can be rewritten as a function of droplet diameter (\(D_d\)), assuming a spherical and homogeneous droplet:

\[\frac{dD_d}{dt} = -\frac{2\dot{m}_d}{\pi \rho_l D_d^2}\]  (4.3)

According to the CEM model, the evaporation flow rate can be described by equation (4.4):

\[\dot{m}_d = \pi D_d D_{vm} \rho_m Sh_m \ln(1 + B_M)\]  (4.4)

where \(D_{vm}\) is the vapor diffusion coefficient, \(\rho_m\) is the density, \(Sh_m\) is the Sherwood number, and \(B_M\) is the Spalding mass transfer number.

It's worth mentioning that all physical properties were calculated in the gas-vapor film region, considering the presence of both gas and water vapor species. The Spalding mass transfer number can be calculated using equation (4.5):

\[B_M = \frac{Y_{vs} - Y_{vg}}{1 - Y_{vs}}\]  (4.5)

where \(Y_{vs}\) and \(Y_{vg}\) refer to the vapor mass fraction at the droplet-gas interface and in the environment away from the droplet surface, respectively.

...

## Assumptions and Limitations

It should be remembered that the model here proposed is characterized by some assumptions and limitations that, according to the literature, can marginally affect the obtained results. However, since the output of the model will be compared to experimental results, it’s worth highlighting the main hypotheses on which the model is based.

1. Gas temperature during evaporation and droplet transport inside the control volume was assumed constant.
2. Heat transfer resistance inside an evaporating droplet was assumed to be negligible.

...
