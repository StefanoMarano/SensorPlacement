=========================
Scientist's Documentation
=========================

This pages describes the motivation and the approach behind the sensor placement algorithm. More details are found in the original publication:

* Stefano Maranò, Donat Fäh, and Yue M. Lu, **"Sensor Placement for the Analysis of Seismic Surface Waves: Sources of Error, Design Criterion and Array Design Algorithms"**,  Geophys. J. Int. (2014) 197 (3): 1566--1581. `Free access online <http://academic.oup.com//gji/article/197/3/1566/652229/Sensor-placement-for-the-analysis-of-seismic?guestAccessKey=fbe9d23c-d168-4b3a-8946-118c73b482ad>`__, `doi:10.1093/gji/ggt489 <https://doi.org/10.1093/gji/ggt489>`__.

The proposed sensor placement algorithm was developed to design planar arrays in seismology. Other possible domains of application include acoustic and radar. The most important assumption behind this work is that the wave should cross the array with a pklane wavefront.

.. contents::
   :local:
   :backlinks: top

Motivation
##########

The presence of outliers can severely downgrade the estimation performance. Outliers are related to the presence of sidelobes in the array response.

In the image below, we evaluate the performance of a maximum likelihood estimator at three different signal-to-noise ratio (SNR) levels. The wavevector of a wave with plane wavefront is estimated using the `WaveDec <http://mercalli.ethz.ch/~marra/WaveDec/>`__ software. The red cross pinpoint the location of the true wavevector. The black crosses pinpoint the estimated wavevector for different noise realizations.

At low SNR, the wavenumber estimates appears to be distributed randomly. The signal is dominated by the strong noise and the estimates carry no information about the true wavevector.

At intermediate SNR, many estimates cluster near the true wavevector. Many other estimates are clustered away from the true wavevector. They are clustered near the sidelobe of the expected loglikelihood function (shown in the background). This large errors are called gross errors or outliers.

At high SNR, all the estimates are clustered near the true wavevector. Their variance is well described by the Cramér--Rao bound.

+---------------------------------------------------------------------------+
|                                                                           |
| .. figure:: images/ThreeRegions.png                                       |
|   :alt: Array layout                                                      |
|   :width: 100%                                                            |
|   :figwidth: 90%                                                          |
|   :figclass: align-center                                                 |
|                                                                           |
|   Wavenumber estimate at low signal-to-noise ratio (SNR), at intermediate |
|   SNR and at high SNR.                                                    |
+---------------------------------------------------------------------------+


Important quantities
####################

Given :math:`N_s` sensor positions :math:`\mathbf{p}_n\in\mathbb{R}^2` we define the sampilng pattern as the sum of :math:`N_s` Dirac delta centered at the sensor positions

.. math::
   :label: eqSamplingPattern
   :nowrap:

   \begin{eqnarray*}
      h(\mathbf{p})=\sum_{n=1}^{N_s}\delta(\mathbf{p}-\mathbf{p}_{n})\,.
   \end{eqnarray*}




Then we consider the two-dimensional Fourier transform of the sampling pattern. Let :math:`\boldsymbol{\kappa}=\left(\kappa_x,\kappa_y\right)^T` denote the wavevector (spatial frequency). The Fourier transform is

.. math::
   :label: eqArrayResponse
   :nowrap:

   \begin{eqnarray*}
      H(\boldsymbol{\kappa}) & = & \mathcal{F}\left\{ h(\mathbf{p}))\right\} \\
                                              & = & \int_{\mathbb{R}^2} \sum_{n=1}^{N_s}\delta(\mathbf{p}-\mathbf{p}_{n}) \textrm{exp}\left(-i \boldsymbol{\kappa}^T \mathbf{p} \right) \,\textrm{d}\mathbf{p} \\
                                              & = & \sum_{n=1}^{N_s} \textrm{exp}\left( -i \boldsymbol{\kappa}^T \mathbf{p}_n \right)\,.
   \end{eqnarray*}
	
The squared complex modulus of the Fourier transform, :math:`\left|H(\boldsymbol{\kappa})\right|^{2}`, is known as array response.

+---------------------------------------------------------------+---------------------------------------------------------------+
|                                                               |                                                               |
| .. figure:: images/SamplingPattern.png                        | .. figure:: images/ArrayResponse.png                          |
|   :alt: Sampling pattern                                      |   :alt: Array response                                        |
|   :width: 90%                                                 |   :width: 90%                                                 |
|   :figwidth: 60%                                              |   :figwidth: 60%                                              |
|   :figclass: align-center                                     |   :figclass: align-center                                     |
|                                                               |                                                               |
|   Sampling pattern :math:`h(\mathbf{p})` of a 5 sensor array. |   The array response                                          |
|                                                               |   :math:`\left|H(\boldsymbol{\kappa})\right|^{2}`.            |
|                                                               |                                                               |
+---------------------------------------------------------------+---------------------------------------------------------------+
 

Relationship with likelihood function
#####################################

Love and Rayleigh likelihood functions are related to the array response. It can be shown that the likelihood function of a Love wave :math:`\ln(p_{\mathbf{Y}}(\mathbf{y}|\boldsymbol{\kappa}))` is realted to the array response as

.. math::
   :label: LoveWave
   :nowrap:

   \begin{eqnarray*}
       \mathop{\mathbb{E}} \left\{ \ln\left(p_{\mathbf{Y}}\left(\mathbf{y}|\boldsymbol{\kappa}\right)\right)\right\} \propto f_{\textrm{L}}(\psi,\breve{\psi})\left|H(\boldsymbol{\kappa}-\breve{\boldsymbol{\kappa}})\right|^{2}\,,
   \end{eqnarray*}

where :math:`\mathop{\mathbb{E}}` denotes the expectation operator, :math:`p_{\mathbf{Y}}\left(\mathbf{y}|\boldsymbol{\kappa}\right)` is the likelihood function. True wave parameters are detoted with the superscript :math:`\breve{}`.

Similarly, for Rayleigh waves, it can be shown that

.. math::
   :label: RayleighWave
   :nowrap:

   \begin{eqnarray*}
   \mathop{\mathbb{E}} \left\{ \ln\left(p_{\mathbf{Y}}\left(\mathbf{y} |\boldsymbol{\kappa},\xi\right)\right)\right\} \propto f_{\textrm{R}}(\psi,\breve{\psi},\xi,\breve{\xi})\left|H(\boldsymbol{\kappa}-\breve{\boldsymbol{\kappa}})\right|^{2}\,.
   \end{eqnarray*}

The functions :math:`f_{\textrm{L}}` and :math:`f_{\textrm{L}}` multiply the array response changing its shape. They reflect the contribution to the likelihood function of the three-component sensors. It is important to stress that these two functions do not depend on the sensor positions. For an acoustic wave measured at a microphone array (scalar sensor), these factor would disappear.

The following pictures show graphically the relationship between likelihood function and array response. Observe the shift corresponding to the true wavevector :math:`\breve{\boldsymbol{\kappa}}` and compare with the array response depicted in the previous section.

+---------------------------------------------------------------+---------------------------------------------------------------+
|                                                               |                                                               |
| .. figure:: images/LL_Love.png                                | .. figure:: images/LL_Rayleigh.png                            |
|   :alt: Love wave likelihood                                  |   :alt: Rayleigh wave likelihood                              |
|   :width: 90%                                                 |   :width: 90%                                                 |
|   :figwidth: 60%                                              |   :figwidth: 60%                                              |
|   :figclass: align-center                                     |   :figclass: align-center                                     |
|                                                               |                                                               |
|   Likelihood function for a Love wave. The true wavenumber    |   Likelihood function of a Rayleigh wave.                     |
|   is :math:`\breve{\boldsymbol{\kappa}}`.                     |                                                               |
+---------------------------------------------------------------+---------------------------------------------------------------+

More details concerning the relationship between likelihood function and array response are given in [Maranò_et_al_2014b]_.

Proposed cost function
######################

.. figure:: images/MinimizationRegion.png
  :alt: Array response
  :align: right
  :width: 80%
  :figwidth: 30%

  The region :math:`\mathcal{K}` on the wavenumber (spatial frequency) plane is defined by :math:`\kappa_{\textrm{min}}` and :math:`\kappa_{\textrm{max}}`.

Our aim is to reduce the sidelobes of the array response in a certain spatial badwidth of interest. The region of interest is the `annulus <https://en.wikipedia.org/wiki/Annulus_(mathematics)>`__ :math:`\mathcal{K}` defined by a minimum and maximum wavenumber, :math:`\kappa_{\textrm{min}}` and :math:`\kappa_{\textrm{max}}`, respectively.

We formulate the following optimization problem, minimizing the largest sidelobe in the region :math:`\mathcal{K}`:

.. math::
   :label: eqContinuosMinimization
   :nowrap:

   \begin{eqnarray*}
       \min_{\mathbf{p}_{1},\mathbf{p}_{2},\ldots,\mathbf{p}_{N_s}}	\max_{\boldsymbol{\kappa}\in\mathcal{K}}\left|H(\boldsymbol{\kappa},\mathbf{p}_{1},\mathbf{p}_{2},\ldots,\mathbf{p}_{N_s})\right| \\
	\mathcal{K}=\{\boldsymbol{\kappa}:\kappa_{\textrm{min}}\leq\left\Vert \boldsymbol{\kappa}\right\Vert _{2}\leq2\kappa_{\textrm{max}}\}\,.
    \end{eqnarray*}

This problem is very difficult to optimize. In fact, the minimization variables :math:`\mathbf{p}_{1},\mathbf{p}_{2},\ldots,\mathbf{p}_{N_s}` appear in the argument of complex exponentials, cf. Eq. :eq:`eqArrayResponse`.

.. WARNING::
   The values :math:`\kappa_{\textrm{min}}` and :math:`\kappa_{\textrm{max}}` define region :math:`\mathcal{K}` where the sidelobes are minimized. They should not be confused with the smallest and largest resolvable wavenumber by the optimized array (i.e., array resolution limits).
   
   The array resolution limits are clearly related with the extent of the region :math:`\mathcal{K}`. The smallest resolvable wavenumber is typically slightly smaller than :math:`\kappa_{\textrm{min}}`. The largest resolvable wavenumber is typically :math:`\kappa_{\textrm{max}}/2`.

Discretization and relaxation
#############################

Instead of dealing with the optimization problem of Eq. :eq:`eqContinuosMinimization` directly, we restrict the possible sensor position to arbitrary discrete locations. We introduce the vector :math:`\mathbf{x}\in\{0,1\}^N` to represent the presence or absence of a sensor at discrete locations.

The discretized problem is

.. math::
   :label: eqDiscreteMinimization
   :nowrap:

   \begin{eqnarray*}
	\min_{\mathbf{x}}  \left\Vert \mathbf{F}\mathbf{x}\right\Vert _{\infty} \\
	\mathbf{A}\mathbf{x}  =  \mathbf{b}  \\
	\mathbf{x}  \in \{0,1\}^{N}
    \end{eqnarray*}
 
where :math:`\mathbf{F}:\mathbb{R}^{N}\to\mathbb{C}^{M}` is a linear operator computing the array response :math:`H` at :math:`M` spatial-frequency points. The `infinity norm <https://en.wikipedia.org/wiki/Uniform_norm>`__ returns the largest complex modulus of the array response, :math:`\left\Vert \mathbf{x}\right\Vert _{\infty}=\max\{|x_{1}|,|x_{2}|,\ldots\}`.

Let :math:`\boldsymbol{\kappa}_m` be the :math:`m`-th spatial frequency and let :math:`\mathbf{p}_n` be the position of the :math:`n`-th possible sensor location. The element :math:`m,n` of :math:`\mathbf{F}` is

.. math::
   :label: eqFourierOpreator
   :nowrap:

   \begin{eqnarray*}
	\left[\mathbf{F}\right]_{m,n} & = & \textrm{exp} \left( -i \boldsymbol{\kappa}^T_m \mathbf{p}_n\right)\,.
    \end{eqnarray*}

A linear constraint specifying the number of sensors :math:`\sum_{n=1}^{N} x_n = N_s` is enforced within :math:`\mathbf{A}\mathbf{x}  =  \mathbf{b}`.

The ojective function in :eq:`eqDiscreteMinimization` is convex, a major improvement from :eq:`eqContinuosMinimization`! The problem is still very hard because of the binary constraint on the vector :math:`\mathbf{x}`.
 
As a last step, we relax the problem. Instead of considering the largest sidelobe in terms of complex modulus (:math:`\left| \mathbf{F}\mathbf{x}\right|`), we consider the absolute value of the real and imaginary parts (:math:`\left|\textrm{Re}\left(\mathbf{F}\mathbf{x}\right)\right|` and :math:`\left|\textrm{Im}\left(\mathbf{F}\mathbf{x}\right)\right|`). With this relaxation, the objective function becomes linear.

The relaxed problem, after introducing the dummy variable :math:`y\in\mathbb{R}` is

.. math::
   :label: eqRelaxedMinimization
   :nowrap:

   \begin{eqnarray*}
	\min_{y}\,	y \\
	\left(\begin{array}{c}
	 \textrm{Re}\left(\mathbf{F}\right)\\
	 \textrm{Im}\left(\mathbf{F}\right)\\
	-\textrm{Re}\left(\mathbf{F}\right)\\
	-\textrm{Im}\left(\mathbf{F}\right)
	\end{array}\right) \mathbf{x} \preceq y\mathbf{1} \\
	\mathbf{A}\mathbf{x}	=\mathbf{b} \\
	y\in\mathbb{R}	,\,\mathbf{x}\in\{0,1\}^{N} \\
    \end{eqnarray*}

where :math:`\mathbf{1}` is a vector of :math:`1` s of length :math:`4M`.

The optimization problem of Eq. :eq:`eqRelaxedMinimization` is addressed numerically as a `mixed integer program (MIP) <https://en.wikipedia.org/wiki/Linear_programming>`__.

How to choose the parameters?
#############################

.. note::
   **TODO** Here we explain how to choose the various parameters. How they affect the results and the computational complexity.

Stretching
##########

.. note::
   **TODO** Explain application of scaling property to stretch a normalized array

Bibliography
############

.. [Maranò_et_al_2014b] Stefano Maranò, Donat Fäh, and Yue M. Lu, **"Sensor Placement for the Analysis of Seismic Surface Waves: Sources of Error, Design Criterion and Array Design Algorithms"**,  Geophys. J. Int. (2014) 197 (3): 1566--1581. `Free access online <http://academic.oup.com//gji/article/197/3/1566/652229/Sensor-placement-for-the-analysis-of-seismic?guestAccessKey=fbe9d23c-d168-4b3a-8946-118c73b482ad>`__.
