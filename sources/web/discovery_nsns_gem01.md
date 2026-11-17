# Stefan Ziegenbalg / Northern Sky Narrowband Survey (NSNS): "OIII arc in Gemini"

Primary URL: https://simg.de/nebulae3/gem-01.html (accessed 2026-11-16; raw HTML data/raw/discovery_web/simg_gem01.html)
Survey index: https://simg.de/nebulae3/index.html ; instruments: https://simg.de/nebulae3/instruments.html ; news: https://simg.de/news.html ; RSS: https://simg.de/rss.xml
Catalogue repository: https://github.com/stzie/simg.de-discoveries (GitHub API, accessed 2026-11-16)
Persistent snapshot: DOI 10.5281/zenodo.17215566 ("OIII arc in Gemini", issued 2025-09-27 per DataCite; zenodo.org returned HTTP 403 to our client). Journal: Ziegenbalg 2025, RNAAS 9, 227, doi:10.3847/2515-5172/adfec7.

## RSS item (first announcement on simg.de)

```
<item>
      <title>Detection of an OIII arc in Gemini</title>
      <link>http://www.simg.de/nebulae3/gem-01.html</link>
      <guid isPermaLink="false">20240628::simg.de/nebulae3/gem-01</guid>
      <pubDate>Fri, 28 Jun 2024 00:00:00 +0200</pubDate>
      <category>Northern Sky Narrowband Survey</category>
      <description><![CDATA[
An arc-shaped emission nebula in the constellation Gemini with dominant OIII lines and an apparent length of about 14&deg; is presented.  The object has already been recorded in far UV surveys. What is new, however, is that this object can be seen in visible light and that it is an OIII emitter. 
<p>
Click on the images for more information and high resolution pictures. 
<p>
<a href="http://www.simg.de/nebulae3/gem-01.html"><img src="http://www.simg.de/nebulae3/gem-01-ohs-400.jpg" alt="Detection of an OIII arc in Gemini" width="400" height="660" style="vertical-align:top; margin:5px 2px"></a>
 <a href="http://www.simg.de/nebulae3/gem-01.html"><img src="http://www.simg.de/nebulae3/gem-01b-ohs-400.jpg" alt="Detection of an OIII arc in Gemini" width="400" height="468" style="vertical-align:top; margin:5px 2px"></a>
]]></description>
    </item>
```

## HTML-commented (not displayed) caption on the instruments page

> Selected results
> 
> Below are selected result. Click on the images for detailed descriptions and other views.
> 
>  HII regions and reflection nebulae in Orion constellation. 
>  In this false color composite Hα is mapped to red, blue continuum is mapped to green, and red continuum is mapped to blue. Reflection nebulae appear green to blue, while HII regions are red.
>  This region of the sky is full of famous nebulae, such as the Orion Nebula (M42), Flame Nebula, Horsehead Nebula, Barnard's Loop, and the Witch Head Nebula.
> 
>  Molecular clouds In constellations Taurus and Perseus in front of HII regions.
>  In this false color composite Hα is mapped to red, blue continuum is mapped to green, and red continuum is mapped to blue. Reflection nebulae appear green to blue, while HII regions are red.
>  The bright reddish structure near the center of the image if the California nebula. The greenish nebula at the bottom are the Pleiades.
> 
>  True color image of the Milky Way form the constellation Vulpecula to Aquila, dominated by massive dark nebulae in front of star clouds.
> 
>  An arc-shaped emission nebula in the constellation Gemini with dominant OIII lines.
>  In this false color composite OIII is mapped to red, Hα is mapped to green and SII is mapped to blue.
>  The object has already been recorded in far UV surveys. What is new, however, is that this object can be seen in visible light and that it is an [OIII] emitter. 
> 
>  Interestingly, the arc has been detected in the preview images of the first clear night after the camera array reached its full configuration.
>  Furthermore, it is the first image that contains data from all three emission lines. So, this detection can be seen as the first light of the full array configuration.

## Catalogue entries as first committed (commit 2e1a2ac0, 2024-06-28T14:45:08Z, file new/k.md)

```
	# Name        |         RA |       DEC |     GLON |     GLAT | diam. | diam. | Angle | Comment
	# ------------+------------+-----------+----------+----------+-------+-------+-------+----------
	K1a           |  7 43 56.7 | +18 38 21 | 201.5168 | +19.7097 | 857.3 | 146.4 |   175 | Large (14° × 2.5°) arc with dominant [OIII] emission lines. In a small region (about 3°, denoted as “K1b”), H-alpha and [SII] are dominant. Probably a supernova remnant.
	K1b           |  7 42 39.7 | +21 04 02 | 199.0072 | +20.3643 | 183.8 | 120.0 |   173 | Part of K1a with dominant H-alpha and [SII] emission lines.
	K1c           |  7 57 38.2 | +17 27 01 | 204.0453 | +22.2458 | 441.8 | 133.1 |     1 | Faint filaments that probably belong to K1a/b and emit light at all three wavelengths.
```

## Rename/resize (commit f8273e1d, 2025-06-11T09:40:37Z, "Renamed: O201.5 -> O200.0, O204.0 -> O205.0")

```
-    O201.5                   |  7 43 56.7 | +18 38 21 | 201.5168 | +19.7097 | 857.3 | 146.4 |   175 | K1a    | 2024-06-28 | Large (14° x 2.5°) arc with dominant [OIII] emissions. In a small region (about 3°, denoted as 'K1b'), H-alpha and [SII] are dominant. Probably a supernova remnant.
-    (O201.5a)                |  7 42 39.7 | +21 04 02 | 199.0072 | +20.3643 | 183.8 | 120.0 |   173 | K1b    | 2024-06-28 | Part of K1a with dominant H-alpha and [SII] emissions.
-    (O201.5b)                |  7 50 53.2 | +12 31 38 | 208.1539 | +18.7685 | 181.4 | 44.21 |   175 |        | 2025-05-22 | Another OIII filament that probably belongs to O201.5.
-    O204.0                   |  7 57 38.2 | +17 27 01 | 204.0453 | +22.2458 | 441.8 | 133.1 |     1 | K1c    | 2024-06-28 | Faint filaments that probably belong to K1a/b and emit light at all three wavelengths.
+    O200.0                   |  7 41 48.4 | +19 56 34 | 200.0324 | +19.7501 |  1015 | 181.5 |   173 | K1a, O201.5  | 2024-06-28 | Large (17° x 2.5°) arc with dominant [OIII] emissions
+    (O200.0a)                |  7 42 39.7 | +21 04 02 | 199.0072 | +20.3643 | 183.8 | 120.0 |   173 | K1b, O201.5a | 2024-06-28 | Part of O200.0 only visible in H-alpha and [SII]
+    (O200.0b)                |  7 50 53.0 | +12 31 23 | 208.1575 | +18.7661 | 181.1 | 43.90 |   175 | O201.5b      | 2025-05-22 | Very faint OIII filament that appear to belong to O200.0
+    O205.0                   |  7 55 18.1 | +16 14 54 | 204.9999 | +21.2606 | 309.9 | 99.86 |   167 | K1c, O204.0  | 2024-06-28 | Filaments that might belong to Monogem Ring and emit at all three wavelengths
```

## Repository commit log (latest 24)

- 2025-09-25T12:31:11Z 58c39561 Fixed object O026.2
- 2025-07-01T13:26:38Z 9de9652c Removed: B069, B075, B121, O117.3, O248.0
- 2025-06-11T09:40:37Z f8273e1d Renamed: O201.5 -> O200.0, O204.0 -> O205.0
- 2025-05-28T14:16:31Z 14f24869 Updated decriptions
- 2025-05-28T12:20:50Z 0f2eb2a7 updated: H120.0, O109.4, O100.8
- 2025-05-25T20:55:35Z 4c0d8fb9 Added B111+2. Updated comments for B111+3, B111+6.
- 2025-05-22T19:31:30Z 79396a27 * Added numerous OIII objects
- 2024-11-12T13:00:30Z 32bd57c9 Updated object H106.4 due to avialibility of new data
- 2024-11-01T20:22:05Z 7a8e6837 Added new object and possible ionization source.
- 2024-10-27T19:25:02Z ec14175d * Added new objects
- 2024-10-26T20:58:39Z 3f04251d Added objects
- 2024-10-25T19:51:05Z 6e9e73da Added new objects
- 2024-10-24T22:01:00Z 5f08354f Fixed typos.
- 2024-10-24T21:49:11Z 0dabd101 Added new objects
- 2024-10-24T12:18:44Z 8e776a2e Added new objects.
- 2024-10-21T18:46:32Z 25c1d592 Revised list of uncataloged objects
- 2024-10-17T10:57:01Z acb4d992 updated links
- 2024-10-17T10:51:21Z 13b85226 fixed typo, updated FA1, added GA6
- 2024-06-28T14:47:08Z 7d0a21d7 Fixed syntax errors.
- 2024-06-28T14:45:08Z 2e1a2ac0 Added FA4-FA6 and K1-K5
- 2023-12-04T13:38:09Z 45835786 Added objects J1-J16 and JA1. Updated description for H2 and H6.
- 2023-11-23T01:53:39Z c99f818d Added E6 (confirmed at http://www.simg.de/nebulae1/aql-00.html)
- 2023-02-05T22:04:57Z 2fb57e22 Added objects I*
- 2022-12-01T13:46:43Z eb776146 Removed related  from the tables.

## RNAAS 9, 227 (2025) Section 2.1 (verbatim excerpt via iopscience, Tavily extract)

> "A 17°×3°-large O iii arc was detected in the constellation Gemini. Its center position is located at a R.A. of αc = 7h41m48s and a decl. of δc = 19°54′. ... A central region, about 2.5° in size, is only visible in Hα and [S ii]. The northern region (right in the image) can be seen in all observed wavelengths, while the southern part was detected only in [O iii]."
> "The arc partially overlaps with the Monogem Ring supernova remnant (SNR) ... The faint filaments east of the arc ... probably belong to the Monogem Ring. Since their structure differs from that of the arc, and the northern part of the arc clearly lies outside of the Monogem Ring, it is unlikely that the arc is part of the Monogem Ring SNR."
RNAAS dates (Crossref): received 2025-08-01, accepted 2025-08-06, published online 2025-08-28.

## Current page text (verbatim, HTML stripped)

An arc-shaped emission nebula in the constellation Gemini, with dominant [OIII] lines and an apparent length of about 17°, 
was detected in data from the Northern Sky Narrowband Survey .
The object could be a fragment of a supernova remnant from the same progenitor as the pulsar Geminga, and it may interact with the Monogem Ring.
This possibility is discussed below , after a presentation of the images and object data in the next section.

This image is a false color composite where [OIII] (505.7 nm) is mapped to red, Hα is mapped to green, and [SII] (671.7 nm and 673.0 nm) is mapped to blue. All three color channels contain some continuum light, which appears magenta in the image due to the choice of the color factors. 
Stars are partially subtracted (with different factors for each color so that the average star color is white).

The southern (lower) part of the large arc-shaped object is only visible in [OIII] light (red). A smaller part, about 3° in diameter, is bright only in Hα and [SII] light (cyan).
Some fainter filaments east of the arc (left in the image) and the northern part of the arc are visible in all three wavelengths.

The bright object at the lower right is SH2-274.

Click on the image for a full scale version or follow 
this link for an annotated presentation of the mentioned objects.

 Object data

The center ot the OIII arc is located at a right ascension of 7:41:48 h and a declination of 19:56°. 
Its size is about 17° in the north/south direction and 2.5° in the east/west direction.

[OIII]

Hα

[SII]

Luminance measurement regions

The monochrome images show the OIII arc separately in the light of the three emission lines. Stars are partially subtracted.
The color image marks the regions for which intensity measurements were performed (see below).

The following table contains the results of intensity measurements for the regions marked in the upper image.
The average intensity in Rayleighs is given in each case. The region labeled bg was used as a background reference. 

 Region Hα [OIII] [SII]

 K1a 0.02 2.04 0.06

 K1a' 1.68 1.53 1.14

 K1b 2.20 0.63 2.27

 K1c 0.24 0.72 1.54

For the determination of the intensity, some parameters had to be estimated.
The transmission ratio of the optics and filters was assumed to be 75% for [OIII] and Hα and 70% for the [SII] doublet.
The atmospheric extinction under best conditions was assumed to be 20%.
The errors resulting from these estimations only distort the absolute measured values by a factor,
but have no influence on the intensity ratios between the wavelengths. More relevant 
are the errors caused by the background calibration, as they distort the intensity 
ratios between the individual emission lines. 
The accuracy of the background reference (the two regions marked bg in the image above) is limited by the fact that the entire region is full of faint background nebulae.

 Classification and discussion of the superstructure

The arc seems to be a part of a larger structure. To better assess the situation, a wide field view has been created from data of DR0.2 of the Northern Sky Narrowband Survey .

32° × 32° view with data from DR0.2 of the Northern Sky Narrowband Survey . 
Click on the image for a 80° × 60° view with a resolution of up 10″, or 

this link for an annotated presentation of the arc parts.

As above, in this false-color composite 
[OIII] is mapped to red, Hα is mapped to green, and [SII] is mapped to blue. 
The brightest stars are added to visualize the regions contaminated by strong continuum light. 
Unlike the image above, this version contains no other continuum light (apart from the mentioned stars).

The image shows several supernova remnants (SNRs) between the arc and the galactic plane. The large structure is discussed below. 
The three small SNRs in that region are 
 G206.7+5.9 ,
 G209.8+2.3 , and
 G205.5+0.5 .

 Visual light

A closer look at the faint filaments east of the arc (left in the images) reveals that they seem to cross the arc, as sketched in panel (a) of the next figure.
Furthermore, unlike the arc, these filaments have a structure similar to the other SNR (red and blue) near the galactic plane.
Connecting them seems natural. 
The yellow dotted ellipse in panel (b) of the following figure shows the outline of the projection of a hypothetical bubble containing the SNR filaments that appear to belong together.

(a)

(b)

 Panel (a): The eastern (left) filaments appear to cross the arc.

 Panel (b): The SNR outline constructed from visible filaments (yellow) matches well with a single Monogem Ring according Knies et al. (white)
and the double SNR scenario (magenta) from the same authors. In all cases, the northern part of the arc lies outside these SNRs and also has a different structure.
Thus, it is very unlikely that the arc is a filament of these SNRs.

Click on the images for JavaScript presentations. (Plots can be toggled using the '2' key to improve visibility.)

 X-rays

In X-rays, a SNR denoted as the Monogem Ring can be observed. Its progenitor is most likely the same as that of the Pulsar PSR B0656+1.
This ring is plotted in white in panel (b) of the figure above. Size and position are taken from Knies et al. (2024) , who analyzed the region 
using the latest X-ray data from the eROSITA telescope array. The same authors also propose the existence of a second SNR in a similar distance and with a similar size.
This two SNR scenario is plotted in magenta in the image above, where the eastern (left) SNR is designated G205.6+12.4 while the western one is still denoted as Monogem Ring .

Both scenarios are a good match to the visual structures, including the faint filaments east of the arc.
However, the arc does not appear to be part of the Monogem Ring or G205.6+12.4.

 Far UV

The arc can also be seen in Galex far-UV (134 nm to 179 nm) images as a thin filament.
However, since the data is incomplete and the region is full of similar structures, the arc cannot be recognized as a separate object in Galex images.

Far UV spectroscopy was performed by Kim et al. in the range of 90 nm to 175 nm. 
In the low resolution CIV (154.8 nm and 155.0 nm) emission lines image (Fig. 1 in their work), only the region corresponding to the southern part of the arc (at a galactic longitude from 201° to 208°) becomes visible.
In the paper, this region is denoted as R1 , and it approximately coincides with the region where the arc apparently overlaps with the Monogem Ring or G205.6+12.4, respectively.
Due to the low resolution, it is uncertain whether these emissions originate from the SNR or are the result of some interaction. 
It is very unlikely that the emissions come directly from the arc, because the northern part (at a galactic longitude from 193° to 201°) is not visible.
Furthermore, the authors detected only CIII and CIV emissions, but no OIII. In visible light, however, [OIII] is dominant.
The authors interpret this structure as “the blast wave [of the Monogem SNR] with an isolated cloud”.

 Interaction between Monogem Ring and OIII arc?

There is a noticeable similarity between the emission line intensities of the arc and the position of the Monogem Ring or G205.6+12.4:

 In the southern part, where the arc and the Monogem Ring apparently overlap, the arc is visible only in [OIII].
 Near the apparent boundary of the Monogem Ring, there is a region where only Hα and [SII] emissions are detected. In this region, the arc also appears distorted.
 In the northern part, which lies outside the Monogem Ring, all three observed emission lines are visible.

If the ionization of the arc is caused by the cooling (and subsequent recombination) of hot plasma, this correlation could be explained by differences in temperature.
It should be mentioned that interaction can only take place within the Monogem Ring or on its shell, because the Monogem Ring expands with hypersonic velocity.

 3D structure of the arc

The OIII arc is probably (a patch of a) thin shell. What we see is the two-dimensional projection of this structure.
Even if the shell is assumed to be spherical, it is difficult to reconstruct it from only a small fragment.
However, due to the high aspect ratio, the projected arc must lie close to the perimeter of the projected sphere.

A rough estimate of the region where the center of a spherical shell could be located is indicated by the yellow dotted ellipse in the figure below.
The yellow square marks the center that would be expected if the arc lies exactly on the apparent boundary of a spherical shell.

Possible superstructures of the OIII arc. Click on the image for an annotated JavaScript presentation.

 Yellow: Assuming the arc is part of a thin spherical shell, the dotted ellipse marks the region of likely sphere centers.
Because of the high aspect ratio, the arc is expected to lie close to the apparent boundary.
The yellow square indicates the center if the arc lies almost exactly on that boundary.

 Magenta: The arrow shows the motion of the pulsar Geminga, starting at its birthplace and ending at its current position. 
A possible Geminga SNR is marked by the square (center) and dotted circle (outline).

 White: The dotted line marks a filament that cannot be assigned to other structures in the region and therefore may be associated with the possible Geminga SNR.

 Cyan: The dotted circle represents the spherical approximation of the Monogem Ring used for the 3D visualization in the next figure.

 Ghost object?

A search for similar fragments along the approximate position of the spherical shell yielded no results.
One possible explanation is that the object is moving rapidly, so that hypersonic velocities relative to the ambient interstellar medium (ISM) only occur in the forward direction of motion.
A shock front would then form only in that direction. At least a fraction of the ionization energy would come from the kinetic energy inherited from a high-velocity progenitor.

Planetary nebulae with these properties are described by Ogle et al. (2025) and are referred to there as ghost planetary nebulae (GPN),
because they would be invisible without a high velocity relative to the ISM.
A SIMBAD search for nearby (<200pc, corresponds to an arc size of <60pc) white dwarfs or hot subdwarfs with high proper motion (>100 mas/yr, corresponds to 10 km/s at 200 pc) within
a 60°-diameter region around the equatorial coordinates (RA,DEC) = (99.7°, 15.5°) yielded no candidates that could plausibly have formed the arc as a GPN.
 5) & (pmra>83) & (pmdec>55) -->

 Geminga (ghost) SNR?

Another object, however, does move in the right direction: the pulsar Geminga. 
The magenta arrow in the previous figure traces its motion, beginning at its birthplace ( Pellizza et al. (2019) and Faherty et al. (2007) )
and ending at its current position. 

Could the OIII arc be a SNR formed by Geminga’s progenitor? Could its shape be explained by Geminga’s motion?
To address these questions, more data must be gathered from the literature.

 Data

 Faherty et al. (2007) estimate Geminga's current distance at 250 (+120  / -70) pc based on HST measurements.
 Pellizza et al. (2019) reported an age of 342 kyr, a maximum progenitor mass of 15 solar masses, and a birthplace distance between 90 and 240 pc.
Unfortunately, these distance estimates rely on an outdated and inaccurate parallax measurement for Geminga’s current position – 157 (+59 / -34) pc (from Caraveo et al., 1996 ; see Faherty et al., 2007 for a discussion).
I therefore adopt the same distance for the birthplace as for the current position, since this would also explain the possible interaction with the Monogem Ring (see below). 
In addition, a 0° inclination of the bow shock appears to provide the best fit in Figure 3 of Caraveo et al. (2003) .
With these parameters, Geminga's apparent motion of 17° corresponds to a travel distance of about 74 pc and a velocity of 210 km/s.

 Modelling

The upper limit of the kinetic energy inherited from the progenitor is therefore less then 6.0·10 48  erg (velocity: 210 km/s; mass: 15 solar masses of the progenitor minus 1.4 solar masses of the neutron star).
This value is far below the characteristic explosion energy of 10 51 erg. 
Thus, the motion should have no significant impact on SNR evolution until most of the explosion energy is radiated away in the snowplough phase.

It is further assumed that the kinetic energy inherited from the progenitor is preserved (i.e., not radiated away through shock ionization) and can be treated independently of the SNR energetics.
This implies that the SNR as a whole will be slowed down as more mass is swept up.

These assumptions, combined with the analytic SNR evolution approximation of Blondin et al. (1998) , form the basis of a simple simulation written in C .
The underlying mathematics and usage instructions are documented in the comments of the source code (see link).

 Simulation results

The remaining parameters used for the simulation are a characteristic explosion energy of 10 51  erg and a particle density of 1.7·10 -2  cm -3 .
The density was fitted to match the position of the OIII arc, as described later in this section. 
The results are shown in magenta in the figure above. The small square marks the simulated center of the SNR, located 31 pc from its birthplace and halfway between the birthplace and Geminga’s current position.

Because the mass swept up by the motion (not by the expansion) is large compared to the initial mass, the SNR is decelerated to a velocity of about 49 km/s, or 23% of the initial speed.
Although this is probably greater than the sound speed in the ISM (depending on temperature), this does not affect visibility, since the expansion speed remains much higher – about 130 km/s.
Thus, the resulting velocity relative to the ambient ISM is likely supersonic across the entire shell, meaning the SNR cannot be categorized as a ghost object.
 6  K, see Knies et al., 2024 ), so that a shockwave could form.
(Interestingly, the speed of soun in pure atomit oxygen is 4 times lower)
-->

The radius of the SNR according to the evolution model from Blondin et al. (1998) is 116 pc, shown as the dotted magenta circle in the figure above.
The three-dimensional situation is sketched in the next figure, which also illustrates the possible interaction with the Monogem Ring. 
With the assumed distances of 300 pc for the Monogem Ring ( Knies et al., 2024 ) and 250 pc for both the birthplace and current position of Geminga (and thus also the center of the Geminga SNR),
this is a perfect match: The [OIII]-only part lies within the Monogem Ring, the Hα+[SII] part lies on the boundary, and the region emitting all three wavelengths lies outside.
However, this alignment depends strongly on the uncertain distances. For example, adopting the (unlikely) birthplace distances from Pellizza et al. (2019) would place the entire arc outside the Monogem Ring.

The perfect match was also obtained by fitting the ISM density to 1.7·10 -2  cm -3 .
With a density of 5·10 -3  cm -3 (density around the Monogem Ring derived by Knies et al., 2024 ), the SNR radius would be considerably larger – about 142 pc – and the travel distance increase to 40 pc.

The simulation results depend mainly on ISM density and distance (and only weakly on explosion energy, progenitor mass, and velocity).
The critical parameters are highly uncertain. Within reasonable limits, a change in one parameter can usually be compensated by adjusting the others. 
The same applies to the models themselves: a more advanced simulation might yield different results, but as long as the parameters remain uncertain, compensating adjustments should be possible.

If the sketch in the following figure correctly represents the 3D situation, the distance of the OIII arc would be about 260 pc. 
At this distance, the length of the arc would be almost 80 pc.

Sketch of the 3D situation showing the possible Geminga SNR and the Monogem Ring.
The plotting plane is defined by the position of the viewer (the Sun) and the centers of the Monogem Ring and the possible Geminga SNR.
These points are marked by small crosses.

Thin black, green, red, and blue lines mark the shells of the SNRs at heights of 0 pc, 25 pc, 50 pc and 75 pc above the plotting plane.
Thick lines overlaid on them indicate the location of the OIII arc at the same height, assuming the arc is a fragment of a thin shell. 

The positions of the OIII arc were derived from the image, with two possible solutions.
The more distant one is considered likely because it explains the possible interaction with the Monogem Ring: the [OIII]-only part (red thick line) lies within the Monogem Ring, 
the Hα+[SII] part (cyan) lies on the boundary, and the part emitting all three wavelengths (magenta) lies outside.

 Discussion and conclusion

There are two Hα+[SII] filaments that cannot be assigned to any known structures in the region and may therefore belong to the Geminga SNR.
These are marked by a white dotted line in the figure above . 
The longer structure crosses the galactic plane and thus cannot be associated with the 

Orion-Eridanus superbubble , which is located below the galactic plane.

It is also possible that some filaments usually attributed to the Orion-Eridanus superbubble actually belong to the Geminga SNR. 
Depending on its 3D orientation (see Pon, 2016 ), the Orion-Eridanus superbubble may overlap with the Geminga SNR as constructed here.
However, the intersection is too small to explain why only a small part of the Geminga SNR is visible.

According to the simulation in Blondin et al. (1998) , the Geminga SNR (as constructed here) is in 
transition between the Sedov phase and the radiative (or snowplough) phase. This transition is characterized by instabilities. 
I could not found analytical models for the temperature in this phase, but
the 2.7·10 5  K predicted by Padmanabhan (2001) for the Sedov phase can likely be treated as an upper limit.
At this temperature, plasma recombines and the SNR should become visible at optical wavelengths. 
Could the instabilities account for the discrepancy between prediction and observation?
Cloud the velocity inherited from the progenitor make the difference, even if it is small compared to the expansion velocity?

Another possible explanation for the discrepancy between prediction and observation is a non-uniform ISM. 
In denser regions, SNRs evolve faster and remain smaller. 
If the SNR exploded into a low-density region toward the Monogem Ring, the OIII arc may be the only part still visible, while the rest has already cooled and possibly merged with the ISM.

An interaction between the Monogem Ring and the Geminga SNR could explain the differing emission line ratios, but not the OIII arc as a whole, since part of the arc always lies outside the Monogem Ring.

In general, the size (about 80 pc for the OIII arc and >200 pc for the entire SNR if uniform density is assumed) should not be surprising.
A look at other galaxies, e.g. M101 and NGC 2303 , reveals [OIII]-only structures in outer (low-density) regions with diameters of more than 400 pc.

In conclusion, it seems plausible that the OIII arc belongs to the Geminga SNR. 
At present, this is also the only explanation for the object. 
Nevertheless, for a secure classification, many questions remain open – especially why only a small fragment is visible.

 Further observations

 Weinberger et al. discovered a [SII]+Hα filament. which is also visible in the wide field image presented above. 
A few similar objects are located in the same region. Some of them are bright in [SII], while others are bright in [OIII].
They appear to be extensions of the elongated southern filament of the Monogem Ring or G205.6+12.4, respectively.
The mentioned objects can be found in 

this Javascript presentation . 

 Image data

Images where captured with a camera array which is described on the instruments page .

Image data are:

 View #1
 View #2

 Center position:
 RA: 7:44h, DEC: 18:30°
 RA: 6:12h, DEC: 14°

 FOV:
 8.75°×14.31° (RA×DEC, through center)
 60°×60° (RA×DEC, through center)

 Orientation:
 North is up
 North is up

 Scale:
 10 arcsec/pixel (in center at full resolution)
 10 arcsec/pixel (in center at full resolution)

 Projection type:
 Stereographic
 Stereographic

 Exposure times (Sum of exposure times of all single exposures used to calculate the image):

 Hα: 80 h
 OIII: 55 h
 SII: 52 h

Data from DR0.2 of the Northern Sky Narrowband Survey 

 Image processing

Both images were processed differently. The second image (60°×60° FOV) was calculated from DR0.2 of the Northern Sky Narrowband Survey . (Follow the link for a detailed description.)
The first image (8.75°×14.31° FOV) is older and is not continuum-subtracted (in contrast to the newer one). Its processing is described as follows.

All image processing steps are deterministic and none of the algorithms use machine learning (often referred to as “AI”), which tends to generate plausible looking fake details.
The software used can be downloaded here .

The image processing steps were:

 Bias correction, dark current subtraction, flatfield correction, noise estimation
 Alignment and brightness calibration using stars from reference image
 Stacking with outlier rejection, background estimation and optimal weighting based on noise estimation
 Star subtraction where star positions and intensities are extracted from continuum images
 Denoising and deconvolution of both components (stars and residual)
 Dynamic range compression using non-linear high-pass filter
 Color composition and tonal curve correction

 Citing / persistent resources

If referencing a volatile webpage is not appropriate, a persistent snapshot of this page is available via the DOI 10.5281/zenodo.17215566 .

 References

John M. Blondin, Eric B. Wright, Kazimierz J. Borkowski, and Stephen P.
 Reynolds.
 Transition to the Radiative Phase in Supernova Remnants.
 ApJ , 500(1):342–354, June 1998.
[  DOI  ]

Patrizia A. Caraveo, Giovanni F. Bignami, Roberto Mignani, and
 Laurence G. Taff.
 Parallax Observations with the Hubble Space Telescope Yield the
 Distance to Geminga.
 Astrophysical Journal Letters , 461:L91, April 1996.
[  DOI  ]

P. A. Caraveo, G. F. Bignami, A. De Luca, S. Mereghetti,
 A. Pellizzoni, R. Mignani, A. Tur, and W. Becker.
 Geminga's Tails: A Pulsar Bow Shock Probing the Interstellar
 Medium.
 Science , 301(5638):1345–1348, September 2003.
[  DOI  ]

Jacqueline Faherty, Frederick M. Walter, and Jay Anderson.
 The trigonometric parallax of the neutron star Geminga.
 Astrophysics and Space Science , 308(1-4):225–230, April 2007.
[  DOI  | 
 http  ]

I. J. Kim, K. W. Min, K. I. Seon, J. W. Park, W. Han, J. H. Park,
 U. W. Nam, J. Edelstein, R. Sankrit, and E. J. Korpela.
 Far-Ultraviolet Observations of the Monogem Ring.
 The Astrophysical Journal , 665(2):L139–L142, August 2007.
[  DOI  | 
 http  ]

Jonathan R. Knies, Manami Sasaki, Werner Becker, Teng Liu, Gabriele Ponti, and
 Paul P. Plucinsky.
 A new understanding of the gemini-monoceros x-ray enhancement from
 discoveries with erosita, 2024.
[  arXiv  | 
 http  ]

Patrick Ogle, Mark Petersen, Tim Schaeffer, Lewis McCallum, Alberto
 Noriega-Crespo, R. Michael Rich, Biny Sebastian, Carl Bjork, Steeve
 Body, Sendhil Chinnasamy, Marcel Dreschsler, Tarun Kottary, Yann
 Sainty, Patrick Sparkman, and Xavier Strottner.
 SDSO1 is a Ghost Planetary Nebula Bow Shock in Front of M31.
 arXiv e-prints , page arXiv:2507.15834, July 2025.
[  DOI  | 
 arXiv  ]

T. Padmanabhan.
 Theoretical Astrophysics, Vol 2 .
 2001.
[  DOI  ]

L. J. Pellizza, R. P. Mignani, I. A. Grenier, and I. F. Mirabel.
 On the local birth place of Geminga.
 Astronomy and Astrophysics , 435(2):625–630, May 2005.
[  DOI  | 
 arXiv  | 
 http  ]

Andy Pon, Bram B. Ochsendorf, João Alves, John Bally, Shantanu
 Basu, and Alexander G. G. M. Tielens.
 Kompaneets Model Fitting of the Orion-Eridanus Superbubble. II.
 Thinking Outside of Barnards Loop.
 ApJ , 827(1):42, August 2016.
[  DOI  | 
 arXiv  | 
 http  ]

R. Weinberger, S. Temporin, and B. Stecklum.
 Detection of an optical filament in the Monogem Ring.
 Astronomy and Astrophysics , 448(3):1095–1100, March 2006.
[  DOI  | 
 arXiv  | 
 http  ]

This work made use of the 
 SIMBAD database , maintained by CDS, Strasbourg, France, and
the Gaia public archive .
 -->
