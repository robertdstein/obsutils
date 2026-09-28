# Validation of `obsutils.spherex`

Comparison against published SPHEREx spectrophotometry, and against independent
catalogue photometry. Run 2026-09-27 against QR2/QR3 as served by IRSA.

## 1. Published comparison

Reference: *SPHEREx as a Frontier for Infrared Transients: Classification of New Galactic
FU Ori Outbursts and Classical Novae*, ApJL (2026),
[arXiv:2602.07129](https://arxiv.org/abs/2602.07129). Per-exposure forced-PSF
spectrophotometry for eight Galactic transients is published at
[Zenodo record 18762027](https://zenodo.org/records/18762027), measured with IRSA's
Spectrophotometry Tool (Tractor forced PSF photometry on Level 2 images).

Three targets were re-measured with `obsutils.spherex` and matched exposure by exposure.
Their timestamps are mid-exposure, ours is `MJD-OBS`, a constant 57.56 s offset.

| target | published points | matched | wavelength offset | wavelength rms |
|---|---|---|---|---|
| WTP16aakapb | 147 | 147 (100%) | +1.99 nm | 1.94 nm (0.080% of lambda) |
| WTP15aaavof | 94 | 94 (100%) | +0.39 nm | 0.53 nm (0.021%) |
| WTP24abacpo | 27 | 27 (100%) | +0.98 nm | 1.04 nm (0.054%) |

**Reproduced exactly:** which exposures cover a target (268/268), and the wavelength each
exposure samples, to under a tenth of a spectral channel (channel width ~35 nm).

**Not reproduced:** the fluxes. Median ratio to the published values is 1.11, 1.80 and -3.46
for the three targets, and the disagreement tracks how much the aperture is background
dominated:

| aperture noise / source flux | n | median ratio | 16-84% |
|---|---|---|---|
| 2-10% | 34 | 1.51 | 1.45 - 1.59 |
| 10-50% | 94 | 1.10 | 0.84 - 1.55 |
| > 50% | 99 | 2.64 | -0.91 - 10.80 |

All eight published targets are Galactic, in star-forming regions with annulus background rms
of 333-971 uJy/pixel, against 21 uJy/pixel for the AT2025abcr field. An 18.5 arcsec aperture
there measures nebulosity and neighbours, not the target. Our formal errors come out 20-45
times larger than the PSF-fit errors, which is the correct error for what an aperture measures
in such a field.

**Conclusion:** these data validate the query and wavelength machinery, and cannot validate
aperture fluxes even in principle. This motivated switching the default to forced PSF
photometry, which is compared against the same data in section 3.

## 2. Independent absolute check, in a clean field

The AT2025abcr field is uncrowded enough to test the flux scale against catalogue photometry.

Unit conversion, checked directly on one image: the solid angle pixel map gives 37.86 arcsec^2
per pixel against a WCS pixel scale of 6.1541 arcsec, and our MJy/sr to uJy factor is 894.3
against 888.9 expected for a 6.15 arcsec pixel.

AllWISE finds four sources inside the r = 18.5 arcsec aperture. Summing them, and using the
host's extended-aperture magnitude rather than its PSF-fit one:

| band | sum of catalogued sources in aperture | this pipeline | ratio |
|---|---|---|---|
| W1 (3.35 um) | 7.82 mJy | 7.62 mJy | 0.97 |
| W2 (4.60 um) | 4.55 mJy | 4.95 mJy | 1.09 |

**The absolute flux scale is right to about 10%.** The apparent 1.5x excess over the host's own
catalogue entry (4.84 mJy at W1) is not an error: a second source 16 arcsec away contributes
2.7 mJy, roughly 35% of the aperture flux. The measured SED is a blend, so any transient signal
is diluted, and a variable neighbour would contaminate the difference spectrum. Check for
neighbours before interpreting a difference.

## 3. Forced PSF photometry against the published values

`method="psf"` renders the ePSF calibration product for the right detector zone at the source's
sub-pixel position and fits its amplitude plus a flat sky. Re-running the three targets above
and matching exposure by exposure, restricted to published SNR > 10 (below that their SNR > 3
cut selects upward fluctuations and biases the comparison):

| n | flux ratio to published | per-exposure scatter | error ratio |
|---|---|---|---|
| 213 | **1.047** | **12.3%** | **1.03** |

For comparison, the fixed-aperture photometry this replaced gave a flux ratio of 1.30 with 89%
scatter and errors 27 times larger at r = 3 pixels, and worse at larger radii. That method has
since been removed from the package.

Per target, PSF fit: WTP16aakapb 0.983 (12% scatter), WTP15aaavof 1.126 (8%).

The selection bias at the faint end is visible as expected, the ratio falling towards zero as
the published SNR approaches their cut:

| published SNR | 3-5 | 5-10 | 10-30 | > 30 |
|---|---|---|---|---|
| median flux ratio | -0.20 | 0.47 | 0.80 | 1.06 |

**The error bars now match IRSA's to 3%**, which is the strongest single check: it says the PSF
model, the variance propagation and the unit conversion are all right. A 5% median flux offset
and 12% per-exposure scatter remain. Part of that is applying the QR3 ePSF calibration to QR2
images, which carry a superseded PSF of their own, and part is that Tractor models neighbours
in these crowded fields while the fit here uses a single flat sky term. For publication-grade
absolute fluxes in a confused field, use IRSA's Spectrophotometry Tool.

## 4. Effect on AT2025abcr

| | aperture r = 3 px | PSF fit |
|---|---|---|
| usable measurements | 320 / 374 | 374 / 374 |
| median per-epoch error | 191 uJy | 54 uJy |
| SED peak | 18.8 mJy at 1.76 um | 1.65 mJy at 1.87 um |
| error inflation factor | 2.0 | 4.2 |
| rms(sigma) of post - pre | 1.09 | 1.14 |
| mean post - pre offset | +119 +- 48 uJy (2.5 sigma) | **-41 +- 25 uJy (-1.6 sigma)** |

The PSF fit measures the nuclear point source; the aperture was measuring the extended host
plus the neighbour, which is why the SED normalisation differs by a factor of ten. Precision per
epoch improves by 3.5 times, and the marginal +120 uJy excess seen with aperture photometry does
not survive: there is no significant transient flux.

## 5. Off-nuclear transients and blending

AT2025abcr is off-nuclear: the host nucleus is 9.4 arcsec from the transient position, which at
6.15 arcsec/pixel is 1.5 pixels, against a PSF FWHM of 5-10 arcsec. A single PSF fitted at the
transient position therefore absorbs nuclear light.

Fitting the nucleus and other catalogued neighbours simultaneously (`--deblend`, or `deblend=True`), which also
tilts the fitted background:

| | single PSF | joint fit + tilted sky |
|---|---|---|
| pre-discovery flux at transient position | 1099 uJy | 772 uJy |
| recovered nucleus flux | - | 4.20 mJy (AllWISE extended: 4.84 mJy) |
| correlation between the two amplitudes | - | -0.13 |

The joint fit removes about 330 uJy of leakage and recovers the nucleus at the right
brightness, and the two amplitudes are only weakly correlated, so the fit really does separate
them.

**But the remaining 772 uJy is not transient flux.** There was no point source at that position
before discovery, so the pre-discovery amplitude is host galaxy light. A point source plus a
smooth sky cannot describe an extended galaxy at this pixel scale: raising the background to a
quadratic changed the difference by under 2 uJy, which is why the background order is no longer
an option. Only differences between epochs are interpretable.

### The difference has a floor too

Applying the same stacking to the companion amplitudes, which come from the very same fits:

| fitted component | post - pre |
|---|---|
| transient position (off-nuclear) | -60 +- 19 uJy |
| host nucleus, 9.4 arcsec away | +41 +- 34 uJy |
| neighbour, 16.2 arcsec away | +13 +- 31 uJy |

The transient position and the nucleus move by similar amounts in opposite directions while the
more distant neighbour does not move, and their sum is nearly unchanged. That is flux being
partitioned differently between two blended components between epochs, not either of them
varying. Roughly 50 uJy is the systematic floor for this target, and there is no TDE detection
above it.

## 6. End-to-end test on a bright transient

The checks above validate the pieces. This one asks whether the pipeline, run exactly as a user
would, recovers an obvious transient.

**V1935 Cen (Nova Centauri 2025)**, discovered 2025-09-22 at V = 6.2, MJD 60940. A nova is the
cleanest possible test: there is no host galaxy, the progenitor is invisible beforehand, and the
amplitude is enormous. It is in the Galactic plane (b = +1.3), so the field is crowded, which
also exercises the companion fitting.

```
spherex-phot --ra 219.340708 --dec -58.794444 --reference-mjd 60940 --deblend --plot
```

416 images cover the position, 129 before the eruption and 287 after. A worked version is in
[spherex_nova.ipynb](spherex_nova.ipynb).

| epoch | t - t_eruption | PSF flux |
|---|---|---|
| MJD 60879 | -61 d | 0.90 mJy |
| MJD 60887 | -53 d | 1.56 mJy |
| MJD 60911 | -29 d | 1.12 mJy |
| MJD 61077 | +137 d | **6.52 mJy** |
| MJD 61085 | +145 d | 5.04 mJy |
| MJD 61111 | +171 d | 3.10 mJy |
| MJD 61248 | +308 d | 3.16 mJy |

A clean rise and decline. Stacked per wavelength, the difference is positive in all 21
populated bins, from +4.9 mJy at 0.8 um to +1.1 mJy at 4.9 um, a smooth red transient SED. The
mean offset is **+2.93 +- 0.06 mJy, 52 sigma**, with 20 of 21 bins individually above 3 sigma.

The pre-eruption stack is the other half of the test. In the mid-infrared it should be zero,
because a quiescent nova progenitor is invisible there, and with the neighbours fitted it is:

| | pre-eruption, blue (< 1.5 um) | pre-eruption, red (> 3.5 um) |
|---|---|---|
| single PSF | +2.93 mJy | -0.26 mJy |
| joint fit + companions | +3.07 mJy | **+0.05 mJy** |

Fitting the companions brings the mid-infrared baseline to zero within the errors. The residual
blue flux is the progenitor system plus unresolved crowding in the plane, which AllWISE does not
fully catalogue.

For contrast, AT2025abcr run the same way gives -0.06 +- 0.02 mJy with **no** bin above 3 sigma.
The pipeline finds a bright transient where there is one and not where there is not.

### A flaw this exposed

The error inflation factor assumes the source does not vary. Measured over all of the nova's
epochs it came out at 17, because the nova's own fading was being counted as noise, and that
suppressed a genuine 52 sigma detection to 3 sigma.

`stack_by_phase` now measures the factor on the pre-reference epochs, where the transient is
absent by construction, falling back to every epoch when that phase covers too few wavelength
bins to be stable. The factor also carries a warning above 6.

| target | calibrated on all epochs | calibrated on pre-reference epochs |
|---|---|---|
| V1935 Cen (varies) | 17.0 | **4.1** |
| AT2025abcr (does not) | 3.0 | **2.8** |

The variable source changes by a factor of four, the non-variable one barely at all, which is
the expected signature. AT2025abcr's conclusion is unchanged either way.

## 7. Archival overlay on the SED plots

`archival_photometry` estimates what a forced fit at a position would have measured before
SPHEREx launched, by summing catalogued sources nearby. It is **off by default**, since it is a
rough reference level rather than a measurement and fails badly in crowded fields; pass
`archival=True`, or `--archival`, to compute it and overlay it on the SED. Each source is weighted by the fraction
of its flux a fit centred on the target picks up, computed from the real ePSF: for a fit that
solves only for amplitude at a fixed position, a neighbour offset by s contributes the
normalised overlap of the two PSFs.

Using the ePSF matters. A Gaussian of the same FWHM badly underestimates the wings:

| separation | Gaussian weight | ePSF weight |
|---|---|---|
| 1.7 arcsec | 0.81 | 0.97 |
| 5.0 arcsec | 0.15 | 0.71 |
| 9.4 arcsec | 0.0014 | 0.21 |
| 16 arcsec | 0.00 | 0.025 |

Bands used, all inside the SPHEREx range: Pan-STARRS i, z, y; 2MASS J, H, Ks; WISE W1, W2.
Pan-STARRS covers declinations above about -30, and is skipped silently for southern targets.

### Using the right magnitudes

The first version summed point-source magnitudes and under-predicted the near infrared by a
factor of two to three. That was a bug, not a limitation: a galaxy is in the point-source
catalogues, but only its core is. The 2MASS point-source entry for the AT2025abcr host is
4.1 mJy at Ks against 10.5 mJy for the whole galaxy.

The estimate now queries the 2MASS extended source catalogue alongside the point-source one and
prefers its magnitudes where a source is in both, and uses the WISE elliptical aperture
magnitudes where WISE flags a source as extended:

| band | point-source magnitudes | extended magnitudes | SPHEREx pre-discovery |
|---|---|---|---|
| 2MASS J (1.24 um) | 0.71 mJy | **1.75 mJy** | 1.38 mJy |
| 2MASS H (1.66 um) | 0.66 mJy | **1.78 mJy** | 1.64 mJy |
| 2MASS Ks (2.16 um) | 0.56 mJy | **1.46 mJy** | 1.13 mJy |
| WISE W1 (3.35 um) | 0.89 mJy | **1.12 mJy** | 0.92 mJy |
| WISE W2 (4.60 um) | 0.54 mJy | **0.68 mJy** | 0.67 mJy |

Agreement is now 2-29% across 1.2 to 4.6 micron, against factors of two to three before.

Pan-STARRS is off by default. Neither of its magnitudes works for a galaxy: the PSF ones
measure only the core, and summing Kron ones double counts, because Pan-STARRS deblends an
extended source into fragments that each carry a near-total Kron magnitude. Switching from PSF
to Kron moved the i band from a ratio of 0.44 to 3.00, with i, z and y disagreeing among
themselves by a factor of two, so both options were unusable.

### Compare it with the single-PSF measurement

The estimate is the total catalogued flux falling in the PSF, which is what one PSF fitted at
the position measures. A joint fit that also solves for companions has deliberately removed
some of that, so the two are no longer the same quantity:

| band | archival, all sources | SPHEREx single PSF | archival, companions excluded | SPHEREx joint fit |
|---|---|---|---|---|
| 2MASS J | 1.75 | 1.38 | 0.00 | 1.09 |
| 2MASS Ks | 1.46 | 1.13 | 0.00 | 0.84 |
| WISE W2 | 0.68 | 0.67 | 0.00 | 0.35 |

Excluding the companions does not fix it, it makes it worse: their catalogue entries carry the
extended light that physically remains at the target position, so the estimate goes to zero
while the joint fit still measures about a milliJansky of host.

### Where it still fails

V1935 Cen, in the Galactic plane, over-predicts by a factor of two in the near infrared and by
ten or more at W1 and W2, where the fitted sky term absorbs the stellar confusion background
that the catalogue sum counts as source flux. The overlay is a reference level for an
uncrowded field, not a calibration, which is what the "estimated from archival photometry"
label is for.
