# Droplet Evaporation 0D model
## How to install

## How to use it

## Technical backgroud

# Droplet Evaporation Simulation Model

The 0-D model used for the simulation of a single droplet evaporation was developed according to a diffusion model based on mass and energy balance well reported in literature [^35^][^36^], following the classical evaporation model (CEM) proposed by Spalding [^37^] and Godsave [^38^]. This model was further improved following Abramzon-Sirignano approach [^39^][^40^], that enabled a better account for advective mass and energy transport. A detailed description of the governing equations and limitations of this model was recently reviewed by Pinheiro et al. [^41^]. Mass and energy balances can be described according to equations (4.1) and (4.2), respectively:

[![\frac{dm_d}{dt} = -\dot{m}_d \tag{4.1}](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdm_d%7D%7Bdt%7D%20%3D%20-%5Cdot%7Bm%7D_d%20%5Ctag%7B4.1%7D)](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdm_d%7D%7Bdt%7D%20%3D%20-%5Cdot%7Bm%7D_d%20%5Ctag%7B4.1%7D)

[![m_d \cdot c_{pl} \cdot \frac{dT_d}{dt} = -Q_s \tag{4.2}](https://latex.codecogs.com/svg.latex?m_d%20%5Ccdot%20c_%7Bpl%7D%20%5Ccdot%20%5Cfrac%7BdT_d%7D%7Bdt%7D%20%3D%20-Q_s%20%5Ctag%7B4.2%7D)](https://latex.codecogs.com/svg.latex?m_d%20%5Ccdot%20c_%7Bpl%7D%20%5Ccdot%20%5Cfrac%7BdT_d%7D%7Bdt%7D%20%3D%20-Q_s%20%5Ctag%7B4.2%7D)

where \(m_d\) is the droplet mass, \(\dot{m}_d\) is the mass evaporation rate, \(c_{pl}\) is the specific heat capacity of the liquid droplet, \(T_d\) is the droplet temperature, and \(Q_s\) is the net power transferred from the environment, contributing to droplet temperature increase.

Equation (4.1) can be rewritten as a function of droplet diameter (\(D_d\)), under the hypothesis of a spherical and homogeneous droplet:

[![\frac{dD_d}{dt} = -\frac{2 \cdot \dot{m}_d}{\pi \cdot \rho_l \cdot D_d^2} \tag{4.3}](https://latex.codecogs.com/svg.latex?%5Cfrac%7BdD_d%7D%7Bdt%7D%20%3D%20-%5Cfrac%7B2%20%5Ccdot%20%5Cdot%7Bm%7D_d%7D%7B%5Cpi%20%5Ccdot%20%5Crho_l%20%5Ccdot%20D_d%5E2%7D%20%5Ctag%7B4.3%7D)](https://latex.codecogs.com/svg.latex?%5Cfrac%7BdD_d%7D%7Bdt%7D%20%3D%20-%5Cfrac%7B2%20%5Ccdot%20%5Cdot%7Bm%7D_d%7D%7B%5Cpi%20%5Ccdot%20%5Crho_l%20%5Ccdot%20D_d%5E2%7D%20%5Ctag%7B4.3%7D)

According to the CEM model, widely discussed and applied in different fields of research to study droplet evaporation, based on mass and energy balance, the evaporation flow rate can be described according to the following equation:

[![\dot{m}_d = \pi \cdot D_d \cdot D_v \cdot \rho_m \cdot Sh_m \cdot \ln(1 + B_M) \tag{4.4}](https://latex.codecogs.com/svg.latex?%5Cdot%7Bm%7D_d%20%3D%20%5Cpi%20%5Ccdot%20D_d%20%5Ccdot%20D_v%20%5Ccdot%20%5Crho_m%20%5Ccdot%20Sh_m%20%5Ccdot%20%5Cln(1%20%2B%20B_M)%20%5Ctag%7B4.4%7D)](https://latex.codecogs.com/svg.latex?%5Cdot%7Bm%7D_d%20%3D%20%5Cpi%20%5Ccdot%20D_d%20%5Ccdot%20D_v%20%5Ccdot%20%5Crho_m%20%5Ccdot%20Sh_m%20%5Ccdot%20%5Cln(1%20%2B%20B_M)%20%5Ctag%7B4.4%7D)

where \(D_v\) is the vapor diffusion coefficient, \(\rho_m\) is the density, \(Sh_m\) is the Sherwood number, and \(B_M\) is the Spalding mass transfer number. It’s worth mentioning that all the physical properties were calculated in the region of the gas-vapour film, considering the presence of both gas and water vapour species. Moreover, as described below, the temperature effect on the value of physical properties was included.

The Spalding mass transfer number can be calculated using equation (4.5):

[![B_M = \frac{Y_{vs} - Y_{vg}}{1 - Y_{vs}} \tag{4.5}](https://latex.codecogs.com/svg.latex?B_M%20%3D%20%5Cfrac%7BY_%7Bvs%7D%20-%20Y_%7Bvg%7D%7D%7B1%20-%20Y_%7Bvs%7D%7D%20%5Ctag%7B4.5%7D)](https://latex.codecogs.com/svg.latex?B_M%20%3D%20%5Cfrac%7BY_%7Bvs%7D%20-%20Y_%7Bvg%7D%7D%7B1%20-%20Y_%7Bvs%7D%7D%20%5Ctag%7B4.5%7D)

where \(Y_{vs}\) and \(Y_{vg}\) refer to the vapour mass fraction at the droplet-gas interface and in the environment away from the droplet surface, respectively. While \(Y_{vg}\) can be calculated using the known vapour content in the gas phase far away from the droplet surface, \(Y_{vs}\) can be determined applying the classical ideal gas law (Equation (4.6)) under the assumption of liquid-vapour thermodynamic equilibrium at the droplet interface (Raoult’s and Clausius-Clapeyron’s laws, Equation (4.7)):

[![Y_{vs} = \frac{\chi_{vs} W_v}{\chi_{vs} W_v + \chi_{gs} W_g} \tag{4.6}](https://latex.codecogs.com/svg.latex?Y_%7Bvs%7D%20%3D%20%5Cfrac%7B%5Cchi_%7Bvs%7D%20W_v%7D%7B%5Cchi_%7Bvs%7D%20W_v%20+%20%5Cchi_%7Bgs%7D%20W_g%7D%20%5Ctag%7B4.6%7D)](https://latex.codecogs.com/svg.latex?Y_%7Bvs%7D%20%3D%20%5Cfrac%7B%5Cchi_%7Bvs%7D%20W_v%7D%7B%5Cchi_%7Bvs%7D%20W_v%20+%20%5Cchi_%7Bgs%7D%20W_g%7D%20%5Ctag%7B4.6%7D)

[![\chi_{vs}^{eq} = \frac{p_{sat}}{p_g} = \frac{p_{atm}}{p_g} \exp{\left(\frac{L_v M_v}{R} \left(\frac{1}{T_b} - \frac{1}{T_d}\right)\right)} \tag{4.7}](https://latex.codecogs.com/svg.latex?%5Cchi_%7Bvs%7D%5E%7Beq%7D%20%3D%20%5Cfrac%7Bp_%7Bsat%7D%7D%7Bp_g%7D%20%3D%20%5Cfrac%7Bp_%7Batm%7D%7D%7Bp_g%7D%20%5Cexp%7B%5Cleft(%5Cfrac%7BL_v%20M_v%7D%7BR%7D%20%5Cleft(%5Cfrac%71%7BT_b%7D%20-%20%5Cfrac%71%7BT_d%7D%5Cright)%5Cright)%7D%20%5Ctag%7B4.7%7D)](https://latex.codecogs.com/svg.latex?%5Cchi_%7Bvs%7D%5E%7Beq%7D%20%3D%20%5Cfrac%7Bp_%7Bsat%7D%7D%7Bp_g%7D%20%3D%20%5Cfrac%7Bp_%7Batm%7D%7D%7Bp_g%7D%20%5Cexp%7B%5Cleft(%5Cfrac%7BL_v%20M_v%7D%7BR%7D%20%5Cleft(%5Cfrac%71%7BT_b%7D%20-%20%5Cfrac%71%7BT_d%7D%5Cright)%5Cright)%7D%20%5Ctag%7B4.7%7D)

where \(\chi_{vs}\) and \(\chi_{gs}\) are respectively the surface vapour and gas mass fraction of the liquid gas mixture, \(M_v\) is the vapour molar mass fraction, \(W_g\) is the gas molar mass, \(L_v\) is the enthalpy of vaporization of the liquid droplet, \(p_g\) is the gas pressure, \(R\) is the gas constant, \(T_b\) is the liquid boiling temperature at standard condition.

Sherwood number (\(Sh\)) accounts for the increase of mass transfer due to gas motion and droplet slip velocity through Reynolds (\(Re\)) and Schmidt (\(Sc\)) numbers. According to Ranz-Marshall empirical correlation [^42^], Sh number can be evaluated using equation (4.8):

[![Sh_m = 2 + \left(Re_d^{1/2} Sc_m^{1/3}\right) \tag{4.8}](https://latex.codecogs.com/svg.latex?Sh_m%20%3D%202%20%2B%20%5Cleft(Re_d%5E%7B1/2%7D%20Sc_m%5E%7B1/3%7D%5Cright)%20%5Ctag%7B4.8%7D)](https://latex.codecogs.com/svg.latex?Sh_m%20%3D%202%20%2B%20%5Cleft(Re_d%5E%7B1/2%7D%20Sc_m%5E%7B1/3%7D%5Cright)%20%5Ctag%7B4.8%7D)

where

[![Re_d = \frac{\rho_g D_d |u_g - u_d|}{\mu_m} \tag{4.9}](https://latex.codecogs.com/svg.latex?Re_d%20%3D%20%5Cfrac%7B%5Crho_g%20D_d%20%7Cu_g%20-%20u_d%7C%7D%7B%5Cmu_m%7D%20%5Ctag%7B4.9%7D)](https://latex.codecogs.com/svg.latex?Re_d%20%3D%20%5Cfrac%7B%5Crho_g%20D_d%20%7Cu_g%20-%20u_d%7C%7D%

