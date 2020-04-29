## Sensor placement for the analysis of seismic surface waves


This repository stores the code used in the publication 
[Sensor placement for the analysis of seismic surface waves: sources of error, design criterion and array design algorithms](https://doi.org/10.1093/gji/ggt489) published on the Geophysical Journal International. [How to cite us](#bibtex-entry).

Find extensive documentation and examples here: [stefanomarano.github.io/SensorPlacement](https://stefanomarano.github.io/SensorPlacement)


 * [Abstract](#abstract)
 * [Example](#example)
 * [BibTex Entry](#bibtex-entry)

## Abstract

Seismic surface waves can be measured by deploying an array of seismometers on the surface of the earth. The goal of such measurement surveys is, usually, to estimate the velocity of propagation and the direction of arrival of the seismic waves. In this paper, we address the issue of sensor placement for the analysis of seismic surface waves from ambient vibration wavefields. First, we explain in detail how the array geometry affects the mean-squared estimation error of parameters of interest, such as the velocity and direction of propagation, both at low and high signal-to-noise ratios (SNRs). Secondly, we propose a cost function suitable for the design of the array geometry with particular focus on the estimation of the wavenumber of both Love and Rayleigh waves. Thirdly, we present and compare several computational approaches to minimize the proposed cost function. Numerical experiments verify the effectiveness of our cost function and resulting array geometry designs, leading to greatly improved estimation performance in comparison to arbitrary array geometries, both at low and high SNR levels.

## Example

The following picture depicts optimized array layouts for 6, 12, and 15 sensors. On the bottom the corresponding array responses are shown.

![Examplary array layouts][ArrayLayout]

[ArrayLayout]: img/ArrayLayout.png "Examplary array layouts"
 
## BibTex Entry

```
@Article{2014:MarFaeLu,
  Title                    = {Sensor Placement for the Analysis of Seismic Surface Waves: Sources of Error, Design Criterion and Array Design Algorithms},
  Author                   = {Stefano Maran\`o and Donat F\"ah and Yue M. Lu},
  Journal                  = {Geophys. J. Int.},
  Year                     = {2014},
  Month                    = {June},
  Number                   = {3},
  Pages                    = {1566-1581},
  Volume                   = {197},
  Doi                      = {10.1093/gji/ggt489},
}
```
