# A Temporal–Radial Visualization of Past-Light-Cone Observation in FLRW Cosmology

**Author:** Valdemar Rolfsen  
**Status:** Working draft v0.4 — literature-positioned revision for author review  
**Date:** 8 September 2026  
**Working note:** Revised from v0.3 after a targeted literature review. This revision retains the FLRW visualization as the established core, narrows novelty claims, and adds a research program connecting the representation to anisotropic cosmology and local cosmic-web vorticity. The companion change log separates established results, literature comparisons, and proposed extensions.

---

## Abstract

Astronomical observations connect emission events in the past to an observer through null trajectories, rather than through spatial curves on the observer's present constant-time hypersurface. This paper presents a temporal–radial visualization of that distinction in Friedmann–Lemaître–Robertson–Walker (FLRW) cosmology. On a monotonically expanding interval, the scale factor is mapped to drawing radius and comoving radial separation to drawing angle. The resulting image of an incoming light ray moves inward when traced toward earlier emission times. For simply connected, positively curved spatial slices, a hyperspherical construction also represents the intrinsic geometry of each slice. For other spatial curvatures, the same polar drawing is a schematic map of radial spacetime, not an isometric embedding of the full spatial geometry.

We specify consistent scale-factor and distance conventions, derive the physical radial metric in the plotting coordinates, and compare the representation with a conventional conformal-time diagram. Power-law examples show that the appearance of winding depends on both expansion history and plotting scale. A targeted literature review shows that the underlying null-geodesic physics, polar or hyperspherical projections, and even geodesic spiralling in anisotropic cosmologies have substantial precedents. The present construction should therefore not be interpreted as a new cosmological theory or as evidence for universal rotation. Its narrower contribution is an explicit FLRW baseline that may be useful for comparing diagrammatic winding with physically meaningful anisotropic or vortical effects. We outline such an extension using the kinematics of the cosmological congruence and Bianchi models, and distinguish global vorticity from the observed rotation of cosmic-web filaments. No new dynamical law or observational fit is claimed.

## 1. Introduction

Astronomy is often described as looking outward. We point a telescope toward a distant galaxy and speak of looking across space. What reaches the telescope, however, is light from an earlier emission event, not a simultaneous view of a remote point on the present spatial hypersurface. In the homogeneous FLRW background, increasing line-of-sight comoving distance along our past light cone corresponds to earlier emission time. This statement must not be extended indiscriminately to every quantity called a cosmological distance; angular-diameter distance, for example, need not increase monotonically with redshift [1].

The distinction does not originate in cosmic expansion. Finite light-travel time also makes us see the past in a static spacetime. Expansion changes the relationship among emission time, redshift, and distance. A visualization of cosmological observation should make both facts clear.

This paper develops a temporal–radial picture: cosmic history is drawn as a sequence of spatial hypersurfaces, with earlier slices placed at smaller drawing radii on an expanding branch. A past-directed observation curve crosses those slices. The word *inward* describes this mapping, not an invariant direction in spacetime or motion toward a center of the universe.

The central distinction is therefore:

> A spatial curve at fixed cosmic time describes simultaneous separation. A light ray connecting an emission event to an observer follows a null trajectory through spacetime.

This is standard relativistic geometry, not a new physical principle [5,6]. The contribution proposed here is to define the temporal–radial representation precisely, show how standard null propagation appears in it, and identify where it is useful or misleading. No priority claim for the use of polar or hyperspherical spacetime diagrams is made. Related constructions exist in both standard cosmological diagrammatics and less formal visualization work, while genuine geodesic spiralling in anisotropic homogeneous cosmologies has been studied for decades [8,13–15]. Section 12 therefore separates what is established from what remains a proposed research direction.

## 2. Spacetime as a Foliation of Spatial Hypersurfaces

Let \(\Sigma_t\) be a constant-cosmic-time spatial hypersurface. On the cosmological interval considered here, the spacetime has the product structure

\[
\mathcal M\simeq I\times\Sigma,\qquad
\mathcal M=\bigcup_{t\in I}\Sigma_t.
\tag{1}
\]

Cosmic time \(t\) is proper time along the fundamental comoving observers. We use a dimensionless scale factor normalized by

\[
a(t_0)=1,\qquad H(t)=\frac{\dot a(t)}{a(t)},
\tag{2}
\]

where \(t_0\) is the observation time. A dot denotes differentiation with respect to cosmic time. A comoving radial coordinate \(r\) will carry units of length. In the positively curved case, a separate dimensionless coordinate \(\chi=r/R_c\) is useful, with \(R_c\) the present curvature radius. These conventions avoid using one symbol as both a dimensionless scale factor and a length.

The temporal–radial construction assumes \(a(t)>0\) and \(H(t)>0\) on the interval being drawn. Then scale factor is a one-to-one label of cosmic time, although it is not proportional to time in general. Earlier events have smaller drawing radius. In a contracting interval this ordering reverses; across a bounce or recollapse the same scale factor can occur at different times. At a turning point, using scale factor as a differentiable time coordinate can fail. Such histories require separate branches or an explicit time label.

We work with background FLRW geometry, comoving sources and observer, and radial propagation at fixed sky direction. Perturbations, peculiar velocities, lensing, and radiative transfer are not modeled except where their omission is stated. The simple spherical light-cone intersections below are restricted to regions without caustics or multiple global light paths. The discussion uses the simply connected spatial geometries; nontrivial global identifications require additional treatment [5,6,9].

## 3. Relation to FLRW Geometry

With these conventions, a convenient form of the FLRW metric is

\[
ds^2=-c^2dt^2+a^2(t)\left[dr^2+\mathcal S_K^2(r)d\Omega^2\right],
\qquad d\Omega^2=d\theta^2+\sin^2\theta\,d\phi^2.
\tag{3}
\]

Here \(K\) is the constant spatial-curvature parameter in the normalized comoving metric, with units of inverse length squared. Define

\[
\mathcal S_K(r)=
\begin{cases}
R_c\sin(r/R_c), & K>0,\\
r, & K=0,\\
R_c\sinh(r/R_c), & K<0,
\end{cases}
\qquad R_c=|K|^{-1/2}\quad(K\ne0).
\tag{4}
\]

For \(K\ne0\), the physical curvature-radius magnitude at time \(t\) is \(a(t)R_c\). For \(K=0\), there is no finite physical curvature radius. This is the usual FLRW geometry expressed with a dimensional radial coordinate; it is equivalent to the frequently used dimensionless-coordinate convention [5,6,11].

In the simply connected closed case, put

\[
\chi=\frac{r}{R_c},\qquad A(t)=R_c a(t).
\tag{5}
\]

Equation (3) becomes

\[
ds^2=-c^2dt^2+A^2(t)
\left[d\chi^2+\sin^2\chi\,d\Omega^2\right],
\tag{6}
\]

and \(\Sigma_t\) is a round three-sphere of radius \(A(t)\). The observer-centered radial chart has \(0\le\chi\le\pi\), with the usual coordinate degeneracies at its poles; expressions for a single direct radial path below are used before the antipodal cut locus.

The closed form is closest to the hyperspherical picture, but that does not establish positive curvature observationally. As a dated example, the 2020 analysis of Efstathiou and Gratton reported \(\Omega_{k0}=0.0004\pm0.0018\) for its combined data analysis, consistent with spatial flatness [3]. This paper neither updates those constraints nor fits a curvature model. The geometrical distinction between a literal closed-slice embedding and a general visualization is essential regardless of the preferred observational parameters.

## 4. The Embedding Diagram Is Not Extra Physical Space

For the closed case, each individual spatial slice can be embedded in auxiliary Euclidean four-space as

\[
\mathbf X(t,\chi,\theta,\phi)=A(t)\mathbf n(\chi,\theta,\phi),
\tag{7}
\]

where

\[
\mathbf n=
\left(\cos\chi,\,
\sin\chi\cos\theta,\,
\sin\chi\sin\theta\cos\phi,\,
\sin\chi\sin\theta\sin\phi\right).
\tag{8}
\]

At fixed \(t\), the induced Euclidean metric is \(A^2d\Omega_3^2\), exactly the spatial part of Eq. (6). Thus this is an isometric embedding of each closed spatial slice. It is not an isometric embedding of the four-dimensional Lorentzian spacetime. Allowing \(A\) to vary gives the auxiliary Euclidean line element

\[
d\ell_E^2=dA^2+A^2d\Omega_3^2,
\tag{9}
\]

whereas the physical time term is \(-c^2dt^2\). Drawing time as radius does not change that sign or turn time into an accessible external spatial direction.

For flat or negatively curved slices, concentric round hyperspheres cannot preserve the full intrinsic spatial metric in Eq. (3). In those cases, a temporal–radial drawing represents only the chosen mapping of expansion history and radial separation. Its circles are time contours, not claims that complete flat or hyperbolic spaces are round spheres. Section 7 gives a mapping that works for the radial sector of any of these geometries while retaining this distinction.

Two statements should therefore be kept separate: the surrounding drawing space is auxiliary, and the full hyperspherical slice embedding is literal only in the specified closed case. Neither statement implies that the universe expands into a physical higher-dimensional room.

## 5. Expansion as Hypersurface Scaling

Consider a comoving source at fixed coordinate \(r_s\). Its radial proper separation from the observer, measured along a constant-time slice on the selected direct path, is

\[
D(t)=a(t)r_s.
\tag{10}
\]

Differentiation at fixed \(r_s\) gives

\[
\dot D(t)=\dot a(t)r_s=H(t)D(t).
\tag{11}
\]

In the closed notation this is \(D=A\chi_s\) and \(H=\dot A/A\). The relation concerns the Hubble-flow separation of comoving objects, not bound systems or the changing coordinate position of a propagating photon. It follows from homogeneous spatial scaling and does not identify a spatial center [5,6].

The scale factor is not determined by this differentiation. Predicting \(H_0\), acceleration, or an expansion history requires dynamical assumptions and observational constraints. For context, the Planck 2018 base-\(\Lambda\)CDM analysis, published in 2020, inferred

\[
H_0=(67.4\pm0.5)\,\mathrm{km\,s^{-1}\,Mpc^{-1}},
\tag{12}
\]

while the 2022 SH0ES analysis reported the baseline local value

\[
H_0=(73.04\pm1.04)\,\mathrm{km\,s^{-1}\,Mpc^{-1}}.
\tag{13}
\]

These are the results of those specific analyses [2,4], retained as historical examples of the Hubble tension, not as a current measurement survey or evidence for this visualization. Numerically, a chosen value of \(70\,\mathrm{km\,s^{-1}\,Mpc^{-1}}\) is approximately \(0.0716\,\mathrm{Gyr}^{-1}\). This is an instantaneous fractional rate of distance growth, not a fixed percentage growth law throughout cosmic history.

For the simply connected closed three-sphere,

\[
V_3(t)=2\pi^2A^3(t),\qquad
\frac{\dot V_3}{V_3}=3H(t).
\tag{14}
\]

More generally, the volume of a fixed comoving region scales as \(a^3\), whether or not the complete spatial slice has finite volume. The hyperspherical picture illustrates these scaling relations but supplies no additional cause for them.

## 6. The View Curve

The *view curve* is the diagrammatic image of an incoming radial null geodesic. The photon travels futureward from emission to reception; reconstructing its path from the observer means tracing the same trajectory backward. Cosmic time is a convenient parameter here, not in general an affine parameter.

For fixed \(\theta,\phi\), Eq. (3) gives

\[
0=-c^2dt^2+a^2(t)dr^2,\qquad
\frac{dr}{dt}=-\frac{c}{a(t)}.
\tag{15}
\]

The negative sign is the incoming branch when \(t\) increases toward reception. With \(r(t_0)=0\),

\[
r(t)=\int_t^{t_0}\frac{c\,dt'}{a(t')}.
\tag{16}
\]

For an emission event at \(t_e\), the source's line-of-sight comoving distance is \(D_C(t_e)=r(t_e)\). In the closed notation,

\[
\chi(t)=\int_t^{t_0}\frac{c\,dt'}{A(t')},\qquad
D_C(t_e)=R_c\chi(t_e).
\tag{17}
\]

These are the defining relations of the temporal–radial construction, with their units, orientation, and normalization explicit. A nonzero null tangent cannot lie in a spacelike constant-time hypersurface, whose induced metric is positive definite. Crossing time slices is the physical fact; appearing to move inward is the graphical choice.

For comoving emission and observation in the background geometry, \(1+z=a^{-1}(t_e)\). Consequently [1],

\[
D_C(z)=c\int_0^z\frac{dz'}{H(z')}.
\tag{18}
\]

For a prescribed mixture of noninteracting pressureless matter, radiation, spatial curvature, and a cosmological constant, the background expression is

\[
\begin{aligned}
H(z)&=H_0 E(z),\\
E^2(z)&=\Omega_{r0}(1+z)^4+\Omega_{m0}(1+z)^3\\
&\quad+\Omega_{k0}(1+z)^2+\Omega_{\Lambda0}.
\end{aligned}
\tag{19}
\]

All density parameters are evaluated at \(t_0\), sum to unity for this component inventory, and \(\Omega_{k0}=-Kc^2/H_0^2\). The radiation term assumes relativistic species with \(a^{-4}\) scaling, matter has \(a^{-3}\) scaling, and \(\Lambda\) is constant. Precision treatments of species transitioning between relativistic and nonrelativistic behavior require their actual density evolution rather than assigning them one fixed scaling across all epochs [7,10,11]. Equation (19) is an input expansion model, not an output of the diagram.

## 7. The Temporal–Radial Appearance of the View Curve

### 7.1 A plotting map for any spatial curvature

Choose a positive constant reference length \(L\) and define

\[
R(t)=L a(t),\qquad \psi=\frac{r}{L},\qquad
x=R\cos\psi,\qquad y=R\sin\psi.
\tag{20}
\]

Here \(R\) is drawing radius and \(\psi\) is drawing angle. Along the received ray,

\[
\psi(t)=\frac{1}{L}\int_t^{t_0}\frac{c\,dt'}{a(t')}.
\tag{21}
\]

For the closed case with \(L=R_c\), \(R=A\) and \(\psi=\chi\), recovering the original hyperspherical cross-section. In a flat universe, \(L\) is entirely a plotting reference, not a curvature radius. For any curvature, a different \(L\) defines a radial visualization but does not, in general, preserve the full hyperspherical slice embedding.

The physical sky angles \(\theta,\phi\) are suppressed in this two-dimensional representation. The angle \(\psi\) encodes radial separation along a selected direction, not direction on the sky. A comoving worldline is drawn at fixed \(\psi\); a constant-time radial spatial curve is drawn at fixed \(R\). The arriving observer is at \((R,\psi)=(L,0)\), not at the center of the drawing.

On the expanding branch, tracing the light ray backward decreases \(R\) and increases \(\psi\). This may look like inward winding, but neither a complete turn nor a spiral of any specific kind is guaranteed. Changing \(L\) changes the angular sweep without changing \(a(t)\), \(D_C\), or any observable. Figure 3 illustrates this freedom for a flat example.

The polar drawing also discards multiples of \(2\pi\). It is therefore only a local coordinate picture unless the angular range is restricted; longer trajectories need an unwrapped \(\psi\) label or separate branches. Artificial overlap in the drawing is not evidence of closed spatial topology. In an actual closed model, antipodal focusing and global light paths require their own treatment.

### 7.2 The physical metric in the plotting coordinates

On the assumed interval \(H>0\), the relation \(R=La(t)\) can be inverted locally to give \(t=t(R)\). Since \(dR=H R\,dt\) and \(dr=L\,d\psi\), the radial part of Eq. (3) becomes

\[
\boxed{
 ds_{\mathrm{rad}}^2=
 -\frac{c^2}{H^2(t(R))R^2}\,dR^2+R^2\,d\psi^2.
}
\tag{22}
\]

This equation is the coordinate-transformed FLRW radial metric, derived directly from Eqs. (3) and (20). It is not the Euclidean polar metric \(dR^2+R^2d\psi^2\). Radius is temporal in Eq. (22) despite being drawn as an ordinary radius on the page. For completeness, the suppressed angular sector is

\[
a^2(t)\mathcal S_K^2(L\psi)\,d\Omega^2.
\tag{23}
\]

Only in the closed case with \(L=R_c\) does this equal \(R^2\sin^2\psi\,d\Omega^2\). Thus the radial metric transformation does not silently turn flat or hyperbolic spatial slices into spheres.

The incoming null branch in the plotting coordinates satisfies

\[
\frac{d\psi}{dR}=-\frac{c}{H(t(R))R^2}.
\tag{24}
\]

Equations (22)–(24) supply the physical interpretation of the drawn curve. They are derived coordinate relations, not additional dynamical laws.

### 7.3 Analytic examples

To examine shape without fitting cosmological data, prescribe a power-law history

\[
a(t)=\left(\frac{t}{t_0}\right)^p,\qquad
u=\frac{t}{t_0},\qquad 0<p<1.
\tag{25}
\]

Here \(p\) is a power-law exponent, not pressure. Direct integration of Eq. (21) yields

\[
\psi(t)=\frac{ct_0}{L(1-p)}\left[1-u^{1-p}\right],
\qquad R(t)=L u^p.
\tag{26}
\]

The boundary condition \(\psi(t_0)=0\) and derivative \(\dot\psi=-c/(La)\) verify the result. The angular sweep remains finite as \(t\to0^+\), with \(\psi_{\max}=ct_0/[L(1-p)]\). In the same formal history, the line-of-sight particle-horizon integral is \(ct_0/(1-p)\). The finite angular range does not generally describe a logarithmic spiral.

For the distinct coasting case \(p=1\),

\[
\psi(t)=\frac{ct_0}{L}\ln\frac{t_0}{t},\qquad
R(\psi)=L\exp\left[-\frac{L}{ct_0}\psi\right].
\tag{27}
\]

This does produce a logarithmic spiral in the polar drawing. Its formal past-distance integral diverges as \(t\to0^+\), despite the finite interval in cosmic time. These statements follow from the prescribed functions; they do not assert that either history describes the full observed universe or, for arbitrary curvature, solves a specified matter model.

### 7.4 Matched comparison figures

Figures 1 and 2 use the same flat-background illustrative history \(a=u^{2/3}\), with \(L=ct_0\). Define a dimensionless conformal coordinate by \(d\eta=c\,dt/(La)\), choosing its origin so that \(\eta=3u^{1/3}\). Then \(\eta_0=3\) and the incoming ray is \(r/L=3-\eta\). Two emission events on that ray are labeled \(E_A\) and \(E_B\):

| Event | \(t/t_0\) | \(a(t)\) | \(r/L=\psi\) | \(\eta\) |
|:--|--:|--:|--:|--:|
| \(p_0\), reception | \(1\) | \(1\) | \(0\) | \(3\) |
| \(E_A\), later emission | \(27/64\) | \(9/16\) | \(3/4\) | \(9/4\) |
| \(E_B\), earlier emission | \(1/8\) | \(1/4\) | \(3/2\) | \(3/2\) |

The two emissions can be understood as photons received from two different comoving sources on the same line of sight; they need not be interactions of one photon. The source worldlines and a present-slice spatial separation are drawn in both figures. All plotted curves are calculated from the displayed analytic relations. There is no fit to observations.

![Conformal-time representation of the illustrative history. The same incoming radial null trajectory contains the labeled emission events and ends at the observer today. Dotted horizontal contours are constant-time slices, dashed vertical lines are comoving worldlines, and the dash-dot segment lies entirely on the present slice. The arrow traces the null trajectory toward the past. Both axes are dimensionless; this is a one-direction radial diagram, not an all-sky light cone.](figures/figure-1-conformal.svg){width=92%}

![Temporal–radial representation of exactly the same events and curves as Figure 1, using \(R/L=a\) and \(\psi=r/L\). Constant-time contours become circular arcs and comoving worldlines become radial segments; the present-slice separation follows the outer arc rather than the null curve. Only the angular sector needed for this radial patch is drawn. The origin is a formal \(a\to0\) plotting limit, not the observer or a spatial center. The drawing angle is not a sky angle, and the flat spatial geometry is not being embedded as a round sphere.](figures/figure-2-temporal-radial.svg){width=92%}

![Plotting-scale dependence for the same flat history \(a=u^{2/3}\). Each curve uses \(\psi=(ct_0/L)\,3(1-u^{1/3})\), plotted in its own normalized coordinates \((x/L,y/L)\). Only the reference length changes: the physical expansion history and comoving distance at each time are identical. Consequently, different apparent winding does not represent different photon rotation or different cosmologies. Curves stop at \(u=10^{-6}\); their limiting angular sweeps are \(6\), \(3\), and \(1.5\) radians.](figures/figure-3-plotting-scale.svg){width=88%}

## 8. The Observable Horizon

A particle horizon is a limit on causal contact for an assumed cosmological history, not a wall at the edge of space. With the normalization in Eq. (2), its present line-of-sight comoving radius is

\[
D_{\mathrm{ph}}(t_0)=c\int_{t_i}^{t_0}\frac{dt}{a(t)},
\tag{28}
\]

provided this integral converges. The lower limit \(t_i\) specifies the past boundary or beginning of the history being modeled. In a big-bang history with \(a\to0\) and monotonic expansion, this may also be written

\[
D_{\mathrm{ph}}(t_0)=c\int_0^\infty\frac{dz}{H(z)}.
\tag{29}
\]

The second expression is not a substitute for specifying the early-time history. A finite cosmic age alone is insufficient for convergence, as the coasting example in Eq. (27) demonstrates. If the integral is finite and \(a(t)<1\) throughout the earlier interval, Eq. (28) exceeds \(c(t_0-t_i)\), directly because its integrand contains \(1/a(t)>1\).

Several boundaries or scales must be distinguished. The **particle horizon** concerns signals that could have reached us since the adopted past limit. The **last-scattering emission surface** concerns the origin of the primordial CMB photons received now, within a thin-surface approximation. The **Hubble radius**, \(c/H(t)\), is the instantaneous proper separation at which the background recession rate \(HD\) equals \(c\); it is not generally a causal horizon. A future **event horizon**, when one exists, is defined instead by the convergent integral over the remaining future history [8].

In particular, last scattering is not the lower endpoint of Eq. (28). Photon opacity makes the causal reach of a model different from direct electromagnetic access to its earliest epochs. A diagram of the particle horizon does not remove that observational limitation [7,10].

The often-quoted observable-universe radius of roughly \(46\) billion light-years is a present-day particle-horizon distance in standard cosmological examples, not the radius of a physical boundary or of an embedding sphere [8]. Identifying it with \(R_c\) in a closed model would impose \(\Omega_{k0}=-c^2/(H_0^2R_c^2)\), a separate curvature assumption. Identifying it with the arbitrary plotting length \(L\) would merely choose a drawing scale.

In compact models, a radial null-distance integral can extend beyond the range of a unique shortest spatial path. The simple horizon-as-sphere picture then needs global causal analysis rather than repeated use of a single observer-centered radial chart. The diagrams here remain in the direct, caustic-free regime.

## 9. Why the Early Universe Appears in Every Direction

Let \(p_0\) be the observer's reception event. In standard notation, \(J^-(p_0)\) is the **causal past**, containing events connected to \(p_0\) by timelike as well as null curves. It is not only the past null cone [9]. Let \(\mathcal N^-(p_0)\) denote the null cone traced backward from \(p_0\), restricted here to the direct portion before caustics or cut points. On this portion it coincides with \(\partial J^-(p_0)\).

In the idealization that last scattering takes place on a single background hypersurface \(\Sigma_{t_*}\), the observed last-scattering emission surface is

\[
\boxed{
\mathcal S_{\mathrm{LSS}}(p_0)
=\mathcal N^-(p_0)\cap\Sigma_{t_*}
=\partial J^-(p_0)\cap\Sigma_{t_*}.
}
\tag{30}
\]

Under these assumptions, its coordinate radius is

\[
r_*=\int_{t_*}^{t_0}\frac{c\,dt}{a(t)},
\tag{31}
\]

while \(\theta,\phi\) cover the sky. This is a two-dimensional emission locus on the earlier three-dimensional slice. By contrast, \(J^-(p_0)\cap\Sigma_{t_*}\) includes the interior region, and would not isolate the observed thin emission surface. We do not equate the CMB itself with a set of events: the CMB is the received radiation, and Eq. (30) identifies an idealized origin of its primary photons.

Last scattering actually has finite duration. A visibility function, describing the probability distribution of last-scattering time for received photons, replaces the infinitely thin surface in a more complete treatment. Subsequent scattering and gravitational effects also modify the observed radiation [10,12]. Equation (30) is therefore a background, thin-last-scattering approximation, not a full CMB transfer calculation.

Every sky direction supplies a past-directed null trajectory, and in this approximation each reaches the common last-scattering epoch. The resulting emission locus surrounds the observer in angular coordinates without placing that observer at a preferred center of space. It is not a material shell existing around us today. The associated comoving locations can nevertheless be projected to a shell on a present slice as an inferred spatial construction.

This geometry explains the all-sky location of the emission surface. It does not, by itself, explain the CMB's blackbody spectrum, near-uniform temperature, or anisotropy pattern, and it does not solve the horizon problem. Those questions concern thermal history, initial conditions, and perturbation dynamics [7,10]. Nor should the drawing angle in Figure 2 be mistaken for the physical sky angles used here.

## 10. Light-Cone Observation and Cosmological Distance Measures

Observations arrive along null trajectories, but distances inferred from them may describe spatial separations on specified constant-time slices. There is no contradiction between those statements. For a comoving source observed at \(t_e\), Eqs. (10) and (16) give

\[
D_{\mathrm{proper}}(t;t_e)=a(t)D_C(t_e),\qquad
D_{\mathrm{proper}}(t_0;t_e)=D_C(t_e).
\tag{32}
\]

The last equality uses \(a(t_0)=1\) and the selected direct radial path. The present proper distance is inferred from a spacetime model; the received photon did not travel along that present spatial line. This is the distinction illustrated by the present-slice curves and the null curves in Figures 1 and 2.

Other useful measures have different meanings [1]. The light-travel-time distance is \(c(t_0-t_e)\). The transverse comoving distance is \(D_M=\mathcal S_K(D_C)\). On the direct branch, the background angular-diameter and luminosity distances obey

\[
D_A=\frac{D_M}{1+z},\qquad
D_L=(1+z)D_M=(1+z)^2D_A,
\tag{33}
\]

with the luminosity relation assuming photon conservation and transparent propagation in the metric description. These expressions are not generally equal to the line-of-sight comoving distance. In particular, the curvature mapping affects transverse rather than radial distance, and angular-diameter distance need not be monotonic in redshift [1].

The Euclidean arclength of the temporal–radial drawing is another, purely auxiliary quantity. Using Eq. (20),

\[
L_E=\int_{t_e}^{t_0}
\sqrt{\dot R^{\,2}+R^2\dot\psi^{\,2}}\,dt
=\int_{t_e}^{t_0}\sqrt{\dot R^{\,2}+c^2}\,dt,
\tag{34}
\]

because the incoming ray satisfies \(R\dot\psi=-c\), or equivalently \(R|\dot\psi|=c\). The sign does not affect the squared expression. This length depends on the auxiliary plotting scale and is neither \(D_C\) nor the null spacetime interval. Equation (22), not the page's Euclidean geometry, determines which drawn curves are physically null.

## 11. What the Visualization Does Not Claim

The temporal–radial construction adds no independent equation for \(a(t)\). In general relativity, for matter and radiation modeled as a perfect fluid and with a separate cosmological constant, the background equations include [5,6,11]

\[
H^2=\frac{8\pi G}{3}\rho-\frac{Kc^2}{a^2}
+\frac{\Lambda c^2}{3},
\tag{35}
\]

and

\[
\frac{\ddot a}{a}
=-\frac{4\pi G}{3}\left(\rho+\frac{3P}{c^2}\right)
+\frac{\Lambda c^2}{3}.
\tag{36}
\]

Here \(\rho\) is mass-equivalent density (energy density divided by \(c^2\)) and \(P\) is pressure. Both exclude the separately written cosmological-constant contribution. The capital \(P\) also distinguishes pressure from the power-law exponent \(p\) in Section 7. Matter evolution and an equation of state must accompany these equations to specify a history.

The diagram can display the consequences of an assumed accelerating or decelerating history, but it does not explain dark energy or resolve the Hubble tension. Within exact FLRW, the fundamental comoving congruence has vanishing shear and vorticity, so apparent winding in the temporal–radial plot is not evidence for physical rotation. It implies neither an external expansion arena nor a privileged spatial center. Any extension to physical rotation or anisotropy requires a spacetime model beyond FLRW and observables defined independently of the plotting coordinates.

The construction also has explicit representational limits: it uses scale factor as a time label only on an appropriate branch; the polar drawing can identify different unwrapped angles; the full hyperspherical spatial interpretation applies only to the closed case; and the simple observed-surface picture omits finite-width and perturbed light propagation. These are conditions on the use of the visualization, not proposed modifications of FLRW physics.

## 12. Discussion

### 12.1 Comparison with conformal-time diagrams

For the same arbitrary length \(L\) used in Eq. (20), define dimensionless conformal time by

\[
d\eta=\frac{c\,dt}{La(t)}.
\tag{37}
\]

Then the radial metric and incoming ray become

\[
ds_{\mathrm{rad}}^2=L^2a^2(t)\left[-d\eta^2+d\psi^2\right],
\qquad \psi=\eta_0-\eta.
\tag{38}
\]

With equally scaled axes, radial null trajectories are straight at \(45\) degrees in \((\psi,\eta)\). This standard representation makes causal relations particularly transparent [8]. The temporal–radial map instead uses the scale factor directly as radius, displaying simultaneous spatial scaling through the sizes of its constant-time arcs.

Figures 1 and 2 demonstrate both the equivalence and the tradeoff. Their event labels, source worldlines, emission times, and null trajectory correspond exactly. Figure 1 displays null propagation as a straight line, but scale factor must be read from the slice labels or calculated from \(a(\eta)\). Figure 2 displays scale-factor ratios directly as radius ratios, but the Euclidean slope and arclength no longer express the causal metric, and the angular sweep depends on the plotting convention.

Thus the representation is a complement to a conformal diagram, not a correction to it. Nor does the temporal–radial picture establish a new observational equivalence between time and space. It reorganizes information already supplied by the FLRW metric and expansion history.

### 12.2 Relation to the balloon analogy and proposed use

The balloon analogy motivates the closed-slice construction by showing how distances between comoving points can increase without a preferred point on the surface. The temporal–radial picture changes the emphasis: it draws a sequence of slices and an observation trajectory through them, making simultaneity and emission time different graphical objects. For nonclosed geometry, that visual inheritance must not be mistaken for a literal intrinsic embedding.

A plausible use is introductory explanation of the difference between a source's earlier emission event and its inferred present separation. The matched figures provide an example, not a measured improvement in comprehension. A pedagogical evaluation could ask readers to identify the observer, distinguish null from constant-time curves, interpret the CMB emission locus, and judge which aspects change when \(L\) is varied. Comparison with conventional diagrams would then test whether the representation reduces confusion or introduces new misconceptions about spatial closure and centers.

No such evaluation is reported here, and a comprehensive priority review remains outside the present revision. The established distance, horizon, and causal-geometry literature supplies the physical basis [1,8,9]. The intended contribution is the explicit mapping, its illustrative implementation, and a clear account of its limitations.

### 12.3 Related representations and limits on novelty

The physical content of the temporal–radial construction overlaps strongly with established representations. Conformal-time diagrams make radial null propagation especially simple by mapping the FLRW radial metric to a conformally Minkowskian form, as in Eq. (38). Proper-distance diagrams, including the treatment of cosmological horizons by Davis and Lineweaver, explicitly distinguish photon trajectories from simultaneous spatial separations [8]. Closed FLRW slices have long been represented through hyperspherical embeddings, and polar projections of cosmological spacetime have also appeared in pedagogical and independent visualization work.

A particularly close non-peer-reviewed example is Sherard's *Visualizing space-time and null geodesics with co-moving S-1 coordinates in an expanding universe*, which uses a compact circular spatial representation to display null geodesics, Hubble flow, horizons, and observable-universe structure [13]. Its assumptions and transformation are not identical to Eqs. (20)–(24), but its existence narrows any claim that a circular or polar visualization of expanding cosmology is itself novel. A broader priority claim would require a systematic literature review beyond the targeted search undertaken for this revision.

The present representation should consequently be judged by a narrower standard: whether the explicit map

\[
(R,\psi)=(La,\,r/L)
\]

and the corresponding metric in Eq. (22) provide a useful common language for comparing expansion, causal propagation, and departures from FLRW. The strongest potential contribution is not that a photon can be drawn as a winding curve, but that the FLRW curve can serve as a mathematically controlled **zero-vorticity reference representation**.

### 12.4 Diagrammatic winding versus physical geodesic spiralling

Physical spiralling of null geodesics is not unique to this visualization and is known in anisotropic cosmology. Barrow, Juszkiewicz, and Sonoda analyzed geodesic spiralling in nearly isotropic homogeneous cosmological models and related it to observable CMB structure [14]. Bianchi VII\(_h\) models in particular can generate characteristic spiral temperature patterns on the sky. Later Bayesian analyses and Planck studies found no evidence for a physically coupled Bianchi VII\(_h\) component and placed stringent bounds on cosmic vorticity; the Planck 2013 analysis reported approximately

\[
\left(\frac{\omega}{H}\right)_0<8\times10^{-10}
\qquad (95\%\ \mathrm{CL})
\tag{39}
\]

within the class of models tested [15]. These results mean that a visible spiral in a temporal–radial FLRW plot cannot be interpreted as evidence for universal rotation.

The relevant distinction is kinematical. Let \(u^\mu\) be the four-velocity of a cosmological matter congruence and

\[
h_{\mu\nu}=g_{\mu\nu}+\frac{u_\mu u_\nu}{c^2}
\tag{40}
\]

its spatial projector. The vorticity tensor is

\[
\omega_{\mu\nu}
= h_\mu{}^{\alpha}h_\nu{}^{\beta}
\nabla_{[\beta}u_{\alpha]}.
\tag{41}
\]

Exact FLRW has \(\omega_{\mu\nu}=0\) and vanishing shear for the fundamental comoving congruence. An anisotropic or vortical cosmology can instead alter both the matter congruence and photon propagation. A scientifically meaningful extension of the temporal–radial picture would therefore compare the FLRW baseline with null geodesics in a specified non-FLRW metric and ask whether deviations in the drawing can be tied to invariant or observer-tetrad-defined quantities such as shear, vorticity, image rotation, or direction drift.

A schematic decomposition such as

\[
W_{\mathrm{total}}=W_{\mathrm{FLRW}}+\Delta W
\tag{42}
\]

is useful only as notation. Here \(W\) must not be defined as raw coordinate winding, because that quantity can be changed by the plotting scale \(L\) or by coordinate transformations. A future paper would need to define \(\Delta W\) through a physically fixed tetrad or directly measurable angular transport. Establishing such a definition is an open problem, not a result of the present article.

### 12.5 Local cosmic-web rotation is not universal rotation

Large structures can carry angular momentum even in a globally almost-FLRW universe. Simulations have found substantial coherent rotation of matter around cosmic filaments [16], and a 2021 observational analysis reported a statistical redshift signature consistent with filament spin [17]. More recently, Tudorache et al. reported a thin chain of 14 H I-rich galaxies embedded in a larger \(\sim15\) Mpc filament at \(z=0.032\). Galaxies on opposite sides of the filament spine show opposite line-of-sight motions, from which the authors infer a characteristic rotational velocity of approximately \(110\,\mathrm{km\,s^{-1}}\) for the dense central region [18].

This observation concerns **local large-scale-structure dynamics**, not rotation of the universe as a whole. In the weak-field, sub-horizon limit one may decompose the matter velocity schematically as

\[
\mathbf v = H\mathbf r + \mathbf v_{\mathrm{pec}},
\tag{43}
\]

with local peculiar-flow vorticity

\[
\boldsymbol\omega_{\mathrm{pec}}
=\frac{1}{2}\nabla\times\mathbf v_{\mathrm{pec}}.
\tag{44}
\]

A rotating filament can have nonzero local vorticity while the FLRW background has zero global vorticity. Conflating the two would be a category error. Nevertheless, this separation suggests a possible visualization problem: whether a generalized temporal–radial framework can display an FLRW expansion baseline, a localized vortical flow, and the null paths through that region without confusing coordinate winding with physical rotation.

### 12.6 A concrete research program

The literature comparison changes the scientific target of the project. The next step is not to argue that the temporal–radial spiral reveals hidden cosmic rotation. It is to test whether the representation provides a useful **comparative geometry** across increasingly general cosmological models.

A minimal program would proceed in four stages:

1. **FLRW baseline.** Keep Eqs. (20)–(24) as the zero-shear, zero-vorticity reference and quantify explicitly which curve properties depend on the arbitrary plotting scale \(L\).
2. **Anisotropic but nonrotating case.** Apply an analogous construction to Bianchi I, where directional scale factors introduce shear without background vorticity. This tests whether the representation can separate anisotropic expansion from apparent polar winding.
3. **Spiralling homogeneous case.** Treat a Bianchi VII\(_h\) model using a physically defined observer tetrad, reproduce known geodesic-spiralling behavior, and determine whether a coordinate-independent residual relative to the FLRW baseline can be defined [14,15].
4. **Perturbed/local case.** Introduce a localized vortical peculiar-velocity field representative of a filament and ray-trace through a perturbed metric or suitable weak-field approximation. The aim would be to visualize the coexistence of global expansion and local angular momentum, not to reinterpret filament rotation as universal vorticity [16–18].

Success would require more than an attractive plot. The generalized representation would need either to simplify an existing calculation, expose a useful invariant relation, produce a new observable diagnostic, or measurably improve inference or understanding relative to conventional diagrams. Failure to meet one of those criteria would leave the construction primarily pedagogical. This makes the extension falsifiable as a *methodological contribution*: its value can be tested against established alternatives rather than assumed from visual intuition.

## 13. Conclusion

A cosmological observation connects a past emission event to a present observer along a null trajectory. A present-slice spatial distance is a different construction, even when it is inferred from that observation. The temporal–radial visualization expresses the distinction by mapping scale factor to radius and radial comoving separation to angle on a monotonically expanding interval.

Its defining map can be summarized as

\[
\boxed{
R(t)=La(t),\qquad
\psi(t)=\frac{1}{L}\int_t^{t_0}\frac{c\,dt'}{a(t')}.
}
\tag{45}
\]

Earlier points on the received ray then lie inward in the drawing. For a closed model and \(L=R_c\), this specializes to the closed-slice hyperspherical notation \(R=A\), \(\psi=\chi\). For other curvature choices it is a radial visualization, not a round embedding of the complete spatial slices. The physical interpretation comes from the Lorentzian metric in Eq. (22), not from Euclidean distances on the page.

The construction clarifies the intended distinction between light-cone observation and spatial separation, while keeping the particle horizon, last scattering, and the Hubble radius separate. It neither predicts an expansion history nor explains the CMB's thermal properties. Its detailed curve shape is representation-dependent, as the analytic examples and scale comparison show.

The intuitive statement can therefore be retained with a precise qualification:

> In this temporal–radial representation of an expanding interval, looking outward along the past light cone is drawn as looking inward through cosmic history.

This is a visual description of standard spacetime geometry. The literature reviewed in this revision also shows that polar cosmological visualizations and genuine geodesic spiralling both have precedents, so the present FLRW construction should not be presented as evidence for new dynamics or universal rotation. Its more promising role is as a controlled zero-vorticity baseline for comparing expansion-induced diagrammatic winding with physically defined anisotropic or vortical effects. Whether that extension yields a useful invariant diagnostic, computational simplification, or measurable pedagogical advantage remains an open research question.

## References

[1] D. W. Hogg, “Distance measures in cosmology,” arXiv:astro-ph/9905116 (1999; revised 2000). DOI: 10.48550/arXiv.astro-ph/9905116. Sections 4–7 and 10 provide the distance conventions used here.

[2] Planck Collaboration, N. Aghanim et al., “Planck 2018 results. VI. Cosmological parameters,” *Astronomy & Astrophysics* **641**, A6 (2020). DOI: 10.1051/0004-6361/201833910. arXiv:1807.06209.

[3] G. Efstathiou and S. Gratton, “The evidence for a spatially flat Universe,” *Monthly Notices of the Royal Astronomical Society: Letters* **496**, L91–L95 (2020). DOI: 10.1093/mnrasl/slaa093. arXiv:2002.06892.

[4] A. G. Riess et al., “A Comprehensive Measurement of the Local Value of the Hubble Constant with 1 km/s/Mpc Uncertainty from the Hubble Space Telescope and the SH0ES Team,” *The Astrophysical Journal Letters* **934**, L7 (2022). DOI: 10.3847/2041-8213/ac5c5b. arXiv:2112.04510.

[5] S. Weinberg, *Gravitation and Cosmology: Principles and Applications of the General Theory of Relativity*. Wiley (1972).

[6] S. M. Carroll, *Spacetime and Geometry: An Introduction to General Relativity*. Addison-Wesley (2004).

[7] E. W. Kolb and M. S. Turner, *The Early Universe*. Addison-Wesley (1990).

[8] T. M. Davis and C. H. Lineweaver, “Expanding Confusion: Common Misconceptions of Cosmological Horizons and the Superluminal Expansion of the Universe,” *Publications of the Astronomical Society of Australia* **21**, 97–109 (2004). DOI: 10.1071/AS03040. arXiv:astro-ph/0310808.

[9] E. Minguzzi, “Lorentzian causality theory,” *Living Reviews in Relativity* **22**, article 3 (2019). DOI: 10.1007/s41114-019-0019-x. Section 1.11 defines the causal and chronological relations.

[10] W. Hu and S. Dodelson, “Cosmic Microwave Background Anisotropies,” *Annual Review of Astronomy and Astrophysics* **40**, 171–216 (2002). DOI: 10.1146/annurev.astro.40.060401.093926. arXiv:astro-ph/0110414.

[11] D. Tong, *Lectures on General Relativity*, University of Cambridge lecture notes, Chapter 4, Section 4.6, “A Taste of Cosmology.” Used for the FLRW and Friedmann-equation conventions; factors of \(c\) are restored in this paper.


[12] U. Seljak and M. Zaldarriaga, “A Line of Sight Approach to Cosmic Microwave Background Anisotropies,” *The Astrophysical Journal* **469**, 437–444 (1996). DOI: 10.1086/177793. arXiv:astro-ph/9603033. Sections 2 and 3.4 describe the visibility function and the finite time interval of last scattering.

[13] P. K. Sherard, “Visualizing space-time and null geodesics with co-moving S-1 coordinates in an expanding universe,” revised manuscript (3 March 2024). DOI: 10.13140/RG.2.2.35586.79049. Non-peer-reviewed related visualization; cited for comparison rather than as established cosmological theory.

[14] J. D. Barrow, R. Juszkiewicz, and D. H. Sonoda, “Universal rotation: how large can it be?”, *Monthly Notices of the Royal Astronomical Society* **213**, 917–943 (1985). DOI: 10.1093/mnras/213.4.917. Includes analysis of geodesic spiralling in homogeneous anisotropic cosmologies.

[15] Planck Collaboration, P. A. R. Ade et al., “Planck 2013 results. XXVI. Background geometry and topology of the Universe,” *Astronomy & Astrophysics* **571**, A26 (2014). DOI: 10.1051/0004-6361/201321546. Includes a physically coupled Bianchi VII\(_h\) search and constraints on cosmological vorticity.

[16] Q. Xia, M. C. Neyrinck, Y.-C. Cai, and M. A. Aragón-Calvo, “Intergalactic filaments spin,” *Monthly Notices of the Royal Astronomical Society* **506**, 1059–1072 (2021). DOI: 10.1093/mnras/stab1713. Cosmological simulations of coherent filament rotation.

[17] P. Wang, N. I. Libeskind, E. Tempel, X. Kang, and Q. Guo, “Possible observational evidence for cosmic filament spin,” *Nature Astronomy* **5**, 839–845 (2021). DOI: 10.1038/s41550-021-01380-6.

[18] M. N. Tudorache et al., “A 15 Mpc rotating galaxy filament at redshift z = 0.032,” *Monthly Notices of the Royal Astronomical Society* **544**, 4306–4316 (2025). DOI: 10.1093/mnras/staf2005. Reports a rotating filamentary system with an inferred characteristic rotation velocity of approximately 110 km s\(^{-1}\).
