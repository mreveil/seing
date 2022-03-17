================
Fingerprints
================


Below we discuss fingerprints currently implemented in SEING as well as others in the pipeline
for implementation.



Behler-Parinello
-----------------

Behler-Parinello (BP) also called "Gaussian" fingerprints are local fingerprints based on symmetry
functions. Two of such symmetry functions are given by the radial and angular componenets :math:`G^{rad}` and :math:`G^{ang}` below where summations run over all neighbors :math:`j` and :math:`k` separated by distances :math:`R_{ij}` and :math:`R_{ik}` with respect to atom :math:`i` within a cutoff distance :math:`R_c` around :math:`i`. :math:`\theta_{ijk}` is the angle between atoms :math:`i,j and k`. :math:`\eta`, :math:`R_s`, :math:`\lambda` and :math:`\zeta` are parameters whose values are chosen by the user. :math:`f_c` is a cutoff function used to ensure a smooth transition to zero at the :math:`R_c`. For more information, see: [BP]_

.. math::

   G^{rad}_i = \sum_j e^{-\eta(R_{ij}-R_s)^2}f_c(R_{ij})

.. math::

   G^{ang}_i = 2^{1-\zeta}\sum_{j,k\neq i} (1+\lambda \cos \theta_{ijk})^\zeta e^{-\eta(R_{ij}^2+R_{ik}^2+R_{jk}^2)^2}f_c(R_{ij})f_c(R_{ij})f_c(R_{ij})

.. math::

   f_c =
   \begin{cases}
   & 0.5[\cos(\frac{\pi R_{ij}}{R_c})+1]~~\text{for}~~~R_{ij}\leq R_c
   \\
   & 0 ~~~~~~~~~~~~~~~~~~~~~~~~~~      \text{for}~~~ R_{ij} > R_c \\
   \end{cases}


The SEING implementation of the BP fingerprint only requires the parameters :math:`\eta`, :math:`\lambda` and :math:`\zeta` as :math:`R_s` is automatically set to zero.


AGNI
--------

The [AGNI]_ method was developed as a framework for machine learning force field development in which the forces are calculated directly without going through energy predictions. The associated fingerprint is given by :math:`V_{i,\alpha,k}` below where :math:`\alpha` denotes the direction (x,y or z) of the force between atoms :math:`i` and :math:`j` separated by distance :math:`r_{ij}`. The parameter :math:`w` corresponds to the width of Gaussians placed at positions :math:`a_k` within a cutoff distance :math:`R_c`. Similarly to the BP fingerprint, :math:`f_c` is the cutofff function ensuring a smooth transition to zero at :math:`Rc`.

.. math::

   V_{i,\alpha,k} = \sum_{j\neq i} \frac{r_{ij}^\alpha}{r_{ij}} \frac{1}{\sqrt{2\pi w}}e^{-0.5(\frac{r_{ij}-a_k}{w})^2}f_c(r_{ij})

For the SEING implementation of the AGNI fingerprint, gaussian centers are uniformly chosen between 0 and the cutoff disance :math:`R_c`. The only parameter necessary is the dimensionality of the fingerprint with determines the number of such Gaussian centers to generate. This fingerprint *does not* support derivative calculations.


Bispectrum
------------

Bispectrum fingerprints for representaion of chemical environments were proposed by Bartok et al. and are based on teh decomposition of a local atomic density function with respect to 4D spherical harmonics. The bispecturm representation is then build based on the coefficients :math:`c_{m'm}^j` of the decomposition given below. For more information, please consult the original paper on the development of [Bispectrum]_ fingerprints.

.. math::

   B_{j_1,j_2,j} = \sum_{m'_1,m_1=-j_1}^{j_1} \sum_{m'_2,m_2=-j_2}^{j_1} \sum_{m',m=-j}^{j}  c_{m'm}^jC_{j_1m_1j_2m_2}^{jm}C_{j_1m'_1j_2m'_2}^{jm'}c_{m'_1m_1}^{j_1}c_{m'_2m_2}^{j_2}

Within SEING, only the parameter :math:`j_{max}` is needed (suggested value: 5) to generate bispectrum fingerprints. Please note that this type of fingerprints is relatively slow compared to other fingerprints currently implemented. Derivatives are supported.


Zernike
---------

Zernike fingerprints are similar to Bispectrum fingerprints in the sense that they are based on decomposition of a local atomic density wrt basis sets. However, Zernike fingerprints are based on decomposition wrt zernike polynomials and 3D spherical harmonics. The general formula is given below. Please consult this paper which describes the [Zernike]_ method.

.. math::

   \rho(\tilde{r},\theta,\phi) = \sum_{n=0}^{\inf} \sum_l \sum_{m=-l}^l c_{nl}^mZ_{nl}^m(\tilde{r},\theta, \phi) ~~~~ \text{for} ~~ n-l \geq 0

To use the Zernike fingerprint in SEING, only the parameter :math:`n_max` is needed (suggested value: 5). Derivatives are supported.


PRDF (Coming Soon)
-------------------



Contact Matrix (Coming Soon)
-----------------------------


SPRINT (Coming Soon)
---------------------



Citations
------------

.. [BP] J. Behler and M. Parrinello, Phys. Rev. Lett., 2007, 98, 146401.
.. [AGNI] T. D. Huan, R. Batra, J. Chapman, S. Krishnan, L. Chen and R. Ramprasad, npj Comput. Mater., 2017, 3, 89–109.
.. [Bispectrum] A. P. Bartók, M. C. Payne, R. Kondor and G. Csányi, Phys. Rev. Lett., 2010, 104, 136403.
.. [Zernike] A. Khorshidi and A. A. Peterson, Comput. Phys. Commun., 2016, 207, 310–324.


