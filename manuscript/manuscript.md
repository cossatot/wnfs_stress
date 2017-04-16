---
title: "Topographic modulation of fault kinematics in the Himalaya and Tibet"
author:
- name: Richard Styron
  affiliation: Earth Analysis, Seattle, WA

abstract: "In many locations in the Himalaya and Tibet, extensional stepovers on
strike-slip faults occur beneath pre-existing topographic highs. An influential
physical model of orogens, explaining contemporaneous high-elevation normal
faulting and low-elevation reverse faulting, holds that horizontal tectonic
compression is invariant across the orogen while  vertical stress varies with
topography, changing the balance of stresses. This model is two-dimensional and
requires topography to be supported isostatically, and therefore cannot fully
describe strike-slip to normal fault transitions beneath mountains to small
mountain ranges, as this is a three-dimensional deformation field and
topography of this wavelength is supported isostatically. I introduce a 3D
elastic model describing the modulation of fault kinematics by
shorter-wavelength topographic stress, and show how the model can place
constraints on the tectonic stress field. I then calculate the topographic
stress field on the Western Nepal Fault System, and use topographic stresses
and observed fault kinematics to invert for the tectonic stress field. The
results yield a maximum tectonic compression of 0–0.2 rho gz and minimum
tectonic compression of -0.1–0.1 rho gz, and reproduce kinematics from normal,
strike-slip and thrust faults and earthquakes in and around western Nepal,
including the 2015 Gorkha earthquake. This demonstrates that where vertical and
a horizontal principal stress are near equal, 1-10 km scale variations in
topography can change fault kinematics, and that pre-existing topography can
influence the location of subsequent faults and stepovers."

---


# Introduction

In the Himalaya and Tibet, it has been observed for decades that thrust
earthquakes take place at the low-elevation rangefronts surrounding the
plateau, while normal and strike-slip earthquakes occur in the elevated
interior (e.g., *Molnar and Tapponnier, 1978; Elliott et al., 2010*) (Figure
1a). *Molnar and Lyon-Caen (1988)* offered an influential physical explanation
for this simultaneous low-elevation reverse faulting and high-elevation normal
faulting: Horizontal tectonic compression $\sigma_H$ (integrated over the
crustal column) is essentially spatially invariant, while vertical compression
$\sigma_V$ at depth varies with the height of the overlying terrain. At low
elevations, $\sigma_H > \sigma_V$, leading to crustal thickening, while at
higher elevations, $\sigma_V > \sigma_H$, causing crustal extension.  

\begin{figure*}[bt] 
\centering
\includegraphics{./figures/wnfs_tib_map_.pdf}
\label{fig_tib_map}
\caption{Active faulting of the Tibetan plateau. \textbf{1a}: View of the
         orogen. \textbf{1b}: Map of the southeastern KF and WNFS. \textbf{1c}:
         Map of the Tibrikot-Dogari segment of the WNFS. KF=Karakoram Fault.
         ATF=Altyn Tagh Fault. GM=Gurla Mandhata Detachment. G=Gar Basin.
         H=Humla Fault. TF=Talphi Fault. TDF=Tibrikot-Dogari Faults. BGF=Bari
         Gad Fault. Focal mechanisms from the GCMT catalog 
         (\emph{Ekström et al. 2012}). Faults from HimaTibetMap v. 1.2
         (\emph{Styron et al., 2010})}
\end{figure*}

These observations and hypotheses are focused on orogens in scope and scale. 
The lower limit to which they are applicable is a natural question, but not one 
that has received much attention. The physical model of *Molnar and Lyon-Caen 
(1988)* rests on the assumption of isostatically-supported topography, which 
is valid over $10^2$–$10^3$ km scales (i.e., orogens), but may not be over 1-10
km (mountain to mountain range) scales, where topography is mostly supported
elastically (e.g., *Bollinger et al. 2004*). Additionally, topographic slopes at
these scales may impart locally substantial stresses in the upper crust, which
may not be a factor at larger scales. The motivating observations are typically
earthquake focal mechanisms, which are too sparse spatially and too narrow
temporally to describe deformation at 1-10 km resolution.

However, neotectonic mapping provides an alternate description of the 
deformation field averaged over longer timescales ($10^3$-$10^5$ years) and 
offers better spatial coverage and resolution as well. Recent neotectonic maps 
of the Himalaya and Tibet reveal a rich interaction between fault type and
elevation over 1-100 km scales. These observations document changes in the
kinematics of large fault systems where the faults cross relatively smaller
topographic features (ridges to mountain ranges) (Figure 1); as the topographic
relief is often known to be older than the active faults, it is likely that
topography modifies fault kinematics at these scales as well, though the
mechanisms may be different. Documenting this phenomenon at is the first
objective of this research.

If topography is capable of modifying fault kinematics, then the total
(tectonic plus topographic) stress field must resolve on the faults with shear
stress directions that are consistent with the directions of fault slip.
Therefore, if the topographic stresses can be calculated, these stresses and
the fault kinematics may be used to constrain the tectonic stress field.
Precise estimation of tectonic stress is a major advancement for understanding
earth processes. Stress is generally unknown at the order of magnitude level,
despite being a fundamental physical property of the earth, and the primary
control on the earth's deformation. Accurate stress estimates are critical for
understanding the distribution of seismicity in time and space, as well as the
physical properties and evolution of faults and orogens. This quantification is
the second objective of this research. 

# Fault kinematic transitions and topography in Tibet and the Himalaya

Many fault zones within the Tibetan plateau and vicinity show transitions in
fault kinematics with changes in topography (Figure 1, S1–S4). These are best
displayed in transtensional fault zones in the elevated interior of the
plateau, consistent with the hypothesis that Tibet is essentially at the
maximum elevation that can be sustained by horizontal tectonic compression
(e.g., *Molnar and Lyon-Caen, 1988*). For example, the sinistral Longmu
Co–Gozha Co fault system has a major extensional stepover where the fault
crosses the Kunlun Shan (Figure S1), where the 2008 $M_w$ 7.1 Yutian earthquake
(*Xu et al., 2013*) occurred. Left-lateral faulting continues to the northeast
of the high mountains as the faults merge with the Altyn Tagh Fault. Similarly,
the transtensional Yibug Caka and Mugu Purou rifts in central Tibet (Figure S2)
show local extensional stepovers where topography is elevated (*Taylor et al.,
2003*; *Ratschbacher et al., 2011*), and the conjugate strike-slip systems to
their south link to rifts where regional elevation steps higher (Figure S3).

Additionally, most isolated topographic highs on the plateau outside of
transtensional zones are cut by active normal faults that do not extend far
into the lower surroundings; the Gangdese (Figure 1b) and Tanggula
Ranges (Figure S4) are prime examples. In these locations, the horizontal
differential stress may not be sufficient for faulting, but may be overcome by
$\sigma_V$ to cause localized extension.

Perhaps the most clear example of topographic modulation of fault kinematics is
on the southeastern Karakoram–Western Nepal Fault System, which undergoes three
distinct extensional stepovers, one for each instance in which the fault system
intersects pre-existing topographic highs (Figure 1b). 

The Karakoram Fault (KF) is a major dextral strike-slip fault on the boundary
between the northwestern Himalaya and southwestern Tibet. The KF is dominantly
strike-slip through most of its length but has a transtensional zone where it
cuts through the Gangdese Range called the Gar Basin (*Sanchez et al., 2010*),
and terminates at the Pulan Graben (dominated by the Gurla Mandhata Detachment,
GMD) essentially where the KF hits the northern Himalaya (*Murphy et al.,
2002*). Some or all of the KF slip is transferred to the GMD and south into
the Himalaya along the Humla Fault (*Murphy and Burgess, 2006*). 

Though the nature of fault connectivity remains unclear, it is likely that 
dextral slip continues through the Himalayan wedge along the Western Nepal 
Fault System (WNFS). Dextral and normal slip has been observed on the Tibrikot 
and Dogari segments of the WNFS (*Murphy et al., 2014; Silver et al., 2015*).
Additional right-lateral offsets and fault traces have been observed through
reconnaisance mapping of the Talphi and Bari Gad Faults to the northwest and
southeast of the Tibrikot fault (*Nakata, 1989*).

A striking feature of the overall geometry of the KF-GMD-WNFS faults is that
extensional stepovers occur wherever the strike-slip faults encounter locally
high terrain (Figure 1b). The higher terrain in all cases predates the
strike-slip faulting along the KF-WNFS system: The Gangdese Range was a
regional topographic high and sediment source by the Oligocene (*Leary et al.,
2016*), predating the post-middle Miocene sediments in Gar Basin (*Sanchez et
al., 2010*), while the Himalaya was uplifted to near modern elevation by the
Miocene (*Garzione et al., 2000*), before extension. Extension along the GMD
has resulted in great uplift of the footwall (>7700 m), but normal faulting
must cause a net decrease in regional elevation despite local footwall uplift.
Additionally, faults associated with the current dextral-normal slip regime cut
north-dipping brittle and ductile fault fabrics associated with the uplift of
the Himalaya (*Silver et al., 2015*). It may be that the WNFS becomes
transpressive in its southern extent, as the Bari Gad fault nears the active
Himalayan frontal folds and thrusts; this is observed on the conjugate Altyn
Tagh fault system on Tibet's northern margin (*Cowgill et al., 2004*).

# Three-dimensional topographic and tectonic stresses

The regional-scale topography and stress relationship described here share the 
central concept developed in the orogen-scale models, that changes in fault 
style result from variable, primarily vertical topographic stresses superposed 
on relatively invariant, primarily horizontal tectonic stresses (Figure 2).
However, the regional model has some key differences, as well.

\begin{figure}[b!]
\includegraphics{./figures/block_diagram.pdf}
\label{fig_block_diagram}
\caption{Block diagram indicating the stress and topography relationships.
         Horizontal stresses $\sigma_H$ and $\sigma_h$ are essentially constant,
         though vertical stress $\sigma_V$ varies greatly with short-wavelength
         topography, locally modulating fault kinematics.}
\end{figure}

First, this problem is inherently three-dimensional. Faults of all kinematic
types in the Himalaya and Tibet accommodate ~N-S shortening, ~E-W extension, or
both (Figure 1). This allows us to expand the model. Whereas the
two-dimensional model has $\sigma_H > \sigma_V$ (reverse faulting) in the
lowlands and $\sigma_V > \sigma_H$ (normal faulting) in the highlands on
parallel faults, the three-dimensional model has $\sigma_V> \sigma_H> \sigma_h$
in the lowlands, $\sigma_H>\sigma_V>\sigma_h$ (strike-slip faulting) in areas
of moderate elevations, and $\sigma_H> \sigma_h>\sigma_V$ in the highlands
(Figure 2). Adding another dimension to the model allows us to consider faults
of all styles and orientations (not only those striking perpendicularly to the
2D model's cross-section), and to calculate the full 3D stress tensor field.
Futhermore, if $\sigma_h$ and $\sigma_H$ are not near equal, there may be a
large elevation gap between reverse and normal faulting, leading to large
uncertainties in stress estimations in the 2D model (*Richardson and Coblentz,
1994*). By considering all three fault styles and principal stresses, the
uncertainties are much reduced.

The second major difference is that shorter-wavelength topography is supported 
elastically rather than isostatically (e.g., *Bollinger et al., 2004*). This
means that topographic stresses may vary dramatically over short horizontal and
vertical distances, slopes may impart locally strong horizontal stress, and the
perturbation to the stress field produced by topography extends outward and
downward rather than simply being a simply vertical sum of the weight of the
overlying rocks. As a result, short-wavelength topographic stresses may be
spatially variable and resolve very differently along strike and down-dip on a
through-going fault.

Finally, the regional model does not require invariance of tectonic stress over
$10^2-10^3$ km, only over $10^0-10^2$ km. Tectonic stress may change over
longer distances due to changes in boundary conditions (e.g., plate driving
forces) and lithospheric rheology. While *Molnar and Lyon-Caen* (*1988*) do
provide compelling arguments for stress invariance across orogens, such
invariance is not a requirement of this model. Consistent application of stress
inversions using the regional model (as below) in many locations throughout an
orogen may serve as a test of orogen-scale tectonic stress invariance.

## Topographic stress

To test the hypothesis that the topography-fault kinematics relationship is
based on a varying topographic stress superposed on a laterally-invariant
tectonic stress field, I seek to reproduce the observed fault kinematics by
calculating the topographic and tectonic stress fields, and resolving them on
3D fault models in the region, following methods outlined in *Styron and
Hetland* (*2015*). The topographic stress calculations are deterministic, as
there is little uncertainty in the topography and allowable variation in the
Earth's elastic parameters does not meaningfully modify the results. The
tectonic stresses are solved for using a Bayesian inversion scheme. To avoid
overfitting and assess the veracity of results beyond the study region, the
inversion is performed on two relatively well-studied faults, the Gurla
Mandhata and Tibrikot-Dogari faults, and then validated on additional
deformation data in the region: a coseismic slip model from the 2015 $M_w$ 7.8
Gorkha, Nepal earthquake (*Galetzka et al., 2015*) and pre-Gorkha focal
mechanisms throughout the region. The two faults modeled here were selected
because the fault geometry and kinematics are well known through field (*Murphy
et al., 2002; Murphy et al., 2014*) and thermochronological (*McCallister et
al., 2014*) studies; the other faults in the WNFS have not received sufficient
study to model confidently.

The topographic stresses are calculated through elastic halfspace methods 
following *Liu and Zoback (1992)* and *Styron and Hetland (2015)*. These 
methods involve the convolution of vertical and horizontal loading functions
describing the distribution of topographic loading on the halfspace surface
with Green's functions describing the propagation of stresses in the halfspace
from vertical and horizontal point load on the surface. This results in a 3D
array with the 3x3 topographic stress tensor calculated at every point (500 m
horizontal resolution, 1 km vertical resolution) in a ~840x700 km region. The
halfspace surface is set to sea level, though calculations above 1500 m below
sea level are discarded due to concerns of overestimating shallow topographic
stress.

The Gurla Mandhata and Tibrikot-Dogari fault traces are extended to depth based
on contraints from structural data and thermal modeling. The fault surfaces are
made into a triangular mesh, and the stress tensors are then interpolated onto
them using barycentric interpolation from the three vertices and centroid of
each triangular fault patch. The rake of the maximum shear stress on each fault
patch is calculated based on the strike and dip of that fault patch. 

## Tectonic stress

I then solve for the allowable tectonic stresses through a Bayesian inversion, 
seeking to minimize the misfit between the rake of the resolved total stress
tensor (topographic plus tectonic) and the observed slip rake. The tectonic
stress field $T$ is assumed to increase linearly with depth below the halfspace
surface (*Townend and Zoback, 2004*), and so is scaled to be a fraction of
lithostatic pressure below the halfspace surface (i.e., $\rho g z$, where
$\rho$ = 2700 kg m$^{-3}$). $T$ is horizontal, and has three components:
$T_{\max}$, $T_{\min}$ and $T_{az}$ (the azimuth of $T_{\max}$). $T_{\max}$ has
a uniform prior from [0—1) $\rho g z$, $T_{\min}$ has a uniform prior of [-1—1)
$T_{\max}$ ($T_{\min}$ is by definition smaller than $T_{\max}$ and therefore
cannot be independently defined in terms of $\rho g z$), and $T_{az}$ has a
uniform prior of 0°—359°. Each sample of $T$ is then rotated to $T_{N-S}$,
$T_{E-W}$, and $T_{N-E}$ (the horizontal shear stress) and added to the
topographic stress tensor at each point.

For each of 1 million samples, the mean rake misfit $\bar{\lambda^m}$ is
calculated as the mean of the absolute value of the rake differences between
the observed slip rake and modeled maximum shear stress rake. Then, the
relative likelihood of each sample set is calculated as
$p(D|T) = \frac{\exp(\kappa \cos \bar{\lambda^m})}{\exp(\kappa \cos
\bar{\lambda^m_{\max}})}$
where $\kappa$ is a scale term, reflecting the uncertainty in the rake data.

The posteriors $p(T|D)$ are then sampled proportionally to the relative
likelihood, following Bayes' rule: $p(T|D) \propto p(T) \, p(D|T)$.

## Results

Topographic stresses tend to be in the direction of fault slip, particularly
for the dip-slip faults, including the Gorkha rupture plane, which is loaded
in a thrust sense by slope-induced subhorizontal compression. Tectonic stresses
are not working against topography in the Himalaya.

\begin{figure}[tb]
\includegraphics{./figures/stress_hists.pdf}
\label{fig_stress_results} \caption{Results of the stress inversion}
\end{figure}

The results of the tectonic stress inversion are shown in Figure 3. The
maximum posterior values for the joint posterior distribution (i.e., the
location of the highest posterior probability density in the 3-variable space)
are $T_{\max} = 0.1 \, \rho g z$, (~2.7 MPa km$^{-1}$ depth), $T_{\min} = -0.1
\, \rho g z$, and $T_{az}= 20°$. The mean absolute misfit between the observed
and modeled fault rakes for the maximum posterior model is 26°. The 1-D
marginals are somewhat similar; $T_{\max}$ has a mode at < 0.05 $\rho g z$,
$T_{\min}$ has a mode near 0, with a tensile skew, and $T_{az}$ has a mode at
20°, parallel to the direction of the Indo-Asian convergence (*e.g., Gan et
al., 2007*).

These results agree well with regional earthquake data not used in the
inversion: The maximum-likelihood tectonic stresses were scaled to depth and
added to the topographic stress tensor at each point in a coseismic slip model
from the 2015 Gorkha, Nepal earthquake (*Galetzka et al., 2015*) and pre-Gorkha
focal mechanisms from the central Himalaya and southern Tibet (*ISC catalog*).
The total stress tensors were resolved on each fault plane and the predicted
shear stress rake was compared to the observed slip rake. The rakes matched
very well (<30° misfit) for nearly all data points, regardless of whether the
data were from thrust, normal or strike-slip earthquakes (Figure S5). This
confirms the predictive power of this simple model where spatially-varying
topographic stress coupled with depth-scaled tectonic stress control a
complicated deformation field.

To test the effects of topographic stress (versus simply fault geometry) on
replicating the shear stresses on the faults, the stress inversion procedure
was repeated without topographic stresses on the fault planes, while holding
all else constant. The results yield a most-likely model with $T_{\max}=0.05 \,
\rho g z$, $T_{\min}=-1.15 \, \rho g z$, and $T_{az}=18°$; $\bar{\lambda^m}$ is
about 10° higher. The strong tension for $T_{\min}$ is required to induce
normal-sense shear on the extensional stepovers in the absence of strong
vertical compression underneath topography. Though the misfit is acceptable, it
is unclear how orogen-parallel tension greater than $\rho g z$ could be
generated in the Himalaya; block divergence due to variably-oblique convergence
along the curved Himalayan front should induce some tension (*McCaffrey and
Nábělek, 1998*), although $\rho g z$ is quite high. Additionally, unlike
topographic stress, tectonic stress does not predict the location of
extensional stepovers, it simply is able to match the slip rake on the existing
stepovers to some degree.

# Discussion

In Tibet and the Himalaya, topography likely modulates fault kinematics over
~10 km scales by locally changing the relative magnitudes of $\sigma_V$ to
$\sigma_H$ and $\sigma_h$. Pre-existing topographic highs produce high
$\sigma_V$ in the crust beneath, causing extensional stepovers in younger
strike-slip faults cutting through the topography. This phenomenon is only
possible where the larger-scale balance of stresses is such that $\sigma_V >
\sigma_H$ under topographic highs but $\sigma_H > \sigma_V$ in adjacent lower
locations. By computing topographic stress, the orientions and magnitudes of
$\sigma_H$ and $\sigma_h$ can be tightly constrained.

In the study areas, the topographic relief (not necessarily the modern
elevation) predates the current tectonic regime and associated faults, and is
therefore capable of controlling the location of releasing bends in strike-slip
faults, as well as isolated grabens (e.g., in the Gangdese range). This may be
common in orogens with a polyphase or protracted history (yielding enough
paleorelief) which finally reach a broad equivalence between $\sigma_V$ and
$\sigma_H$. However, evidence of this process may be erased with erosion, and
the process may even reverse as stepover-produced topography builds: *Cowgill et
al.* (*2004*) suggest that the kinematics of some stepovers on the Altyn Tagh
fault change due to increasing topography and $\sigma_V$.

The tectonic stresses estimated are, perhaps, low for creating the world's
current highest mountain range. They are significantly lower than those
estimated at $T_{\max} = 0.5 - 1 \, \rho g z$ in eastern Tibet with the same
methods  (*Styron and Hetland, 2015*), or $T_{\max}= 1.0 \, \rho g z$ measured
from the upper 2.5 km of the SAFOD pilot hole on the San Andreas (*Hickman and
Zoback, 1994*). However, they are more in line with other estimates from
orogens those undergoing active thrusting at their base and normal
faulting at elevation: *Bollinger et al.* (*2004*) estimate maximum deviatoric
stress to be 35 MPa at ~10 km depth (~0.13 $\rho g z$) in the Himalaya, based
on similar observations to those presented here; *Richardson and Coblentz*
(*1994*) found $T_{\max}= ~0.1-0.2 \, \rho g z$ to be plausible in the upper
crust of the Andes based on finite element modeling; and *Copley et al.*
(*2009*) found an average $T_{\max}$ of 20 MPa in the upper 15 km of Albania,
resulting in a gradient of 0.1 $\rho g z$ if $T_{\max}$ increases linearly with
depth.

# References
Bollinger, L., Avouac, J. P., Cattin, R., & Pandey, M. R. (2004). Stress
buildup in the Himalaya. Journal of Geophysical Research: Solid Earth,
109(B11).

Copley, A., Boait, F., Hollingsworth, J., Jackson, J., & McKenzie, D. (2009).
Subparallel thrust and normal faulting in Albania and the roles of
gravitational potential energy and rheology contrasts in mountain belts.
Journal of Geophysical Research: Solid Earth, 114(B5).

Elliott, J. R., Walters, R. J., England, P. C., Jackson, J. A., Li, Z., &
Parsons, B. (2010). Extension on the Tibetan plateau: recent normal faulting
measured by InSAR and body wave seismology. Geophysical Journal International,
183(2), 503-535.

Galetzka, J., Melgar, D., Genrich, J.F., Geng, J., Owen, S., Lindsey, E.O., Xu,
X., Bock, Y., Avouac, J.P., Adhikari, L.B. and Upreti, B.N., 2015. Slip pulse
and resonance of the Kathmandu basin during the 2015 Gorkha earthquake, Nepal.
Science, 349(6252), pp.1091-1095.

Garzione, C.N., Dettman, D.L., Quade, J., DeCelles, P.G. and Butler, R.F.,
2000. High times on the Tibetan Plateau: Paleoelevation of the Thakkhola
graben, Nepal. Geology, 28(4), pp.339-342.

International Seismological Centre, On-line Bulletin, http://www.isc.ac.uk,
Internatl. Seismol. Cent., Thatcham, United Kingdom, 2013.

Leary, R., Orme, D.A., Laskowski, A.K., DeCelles, P.G., Kapp, P., Carrapa, B.
and Dettinger, M., 2016. Along-strike diachroneity in deposition of the Kailas
Formation in central southern Tibet: Implications for Indian slab dynamics.
Geosphere, 12(4), pp.1198-1223.

Liu, L. and Zoback, M.D., 1992. The Effect of Topography on the State of Stress
in the Crust: Application to the Site of the Cajon Pass Scientific Drilling
Project. Journal of Geophysical Research, 97(B4), pp.5095-5108.

McCaffrey, R. and Nabelek, J., 1998. Role of oblique convergence in the active
deformation of the Himalayas and southern Tibet plateau. Geology, 26(8),
pp.691-694.

McCallister, A.T., Taylor, M.H., Murphy, M.A., Styron, R.H. and Stockli, D.F.,
2014. Thermochronologic constraints on the late Cenozoic exhumation history of
the Gurla Mandhata metamorphic core complex, Southwestern Tibet. Tectonics,
33(2), pp.27-52.

Molnar, P. and Lyon-Caen, H., 1988. Some simple physical aspects of the
support, structure, and evolution of mountain belts. Geological Society of
America Special Papers, 218, pp.179-208.

Molnar, P. and Tapponnier, P., 1979. Active tectonics of Tibet. Journal of
Geophysical Research: Solid Earth, 83(B11), pp.5361-5375.

Murphy, M.A. and Burgess, W.P., 2006. Geometry, kinematics, and landscape
characteristics of an active transtension zone, Karakoram fault system,
southwest Tibet. Journal of Structural Geology, 28(2), pp.268-283.

Murphy, M.A., Taylor, M.H., Gosse, J., Silver, C.R.P., Whipp, D.M. and
Beaumont, C., 2014. Limit of strain partitioning in the Himalaya marked by
large earthquakes in western Nepal. Nature Geoscience, 7(1), pp.38-42.

Murphy, M.A., Yin, A., Kapp, P., Harrison, T.M., Manning, C.E., Ryerson, F.J.,
Lin, D. and Jinghui, G., 2002. Structural evolution of the Gurla Mandhata
detachment system, southwest Tibet: Implications for the eastward extent of the
Karakoram fault system. Geological Society of America Bulletin, 114(4),
pp.428-447.

Nakata, T., 1989. Active faults of the Himalaya of India and Nepal. Geological
Society of America Special Papers, 232, pp.243-264.

Richardson, R.M. and Coblentz, D.D., 1994. Stress modeling in the Andes:
Constraints on the South American intraplate stress magnitudes. Journal of
Geophysical Research: Solid Earth, 99(B11), pp.22015-22025.

Sanchez, V.I., Murphy, M.A., Dupré, W.R., Ding, L. and Zhang, R., 2010.
Structural evolution of the Neogene Gar Basin, western Tibet: Implications for
releasing bend development and drainage patterns. Geological Society of America
Bulletin, 122(5-6), pp.926-945.

Silver, C.R., Murphy, M.A., Taylor, M.H., Gosse, J. and Baltz, T., 2015.
Neotectonics of the Western Nepal Fault System: Implications for Himalayan
strain partitioning. Tectonics, 34(12), pp.2494-2513.

Styron, R.H. and Hetland, E.A., 2015. The weight of the mountains: Constraints
on tectonic stress, friction, and fluid pressure in the 2008 Wenchuan
earthquake from estimates of topographic loading. Journal of Geophysical
Research: Solid Earth, 120(4), pp.2697-2716.

Styron, R., Taylor, M. and Okoronkwo, K., 2010. Database of active structures
from the Indo-Asian collision. Eos, Transactions American Geophysical Union,
91(20), pp.181-182.

Taylor, M., Yin, A., Ryerson, F.J., Kapp, P. and Ding, L., 2003. Conjugate
strike-slip faulting along the Bangong-Nujiang suture zone accommodates coeval
east-west extension and north-south shortening in the interior of the Tibetan
Plateau. Tectonics, 22(4).

Townend, J. and Zoback, M.D., 2000. How faulting keeps the crust strong.
Geology, 28(5), pp.399-402.

