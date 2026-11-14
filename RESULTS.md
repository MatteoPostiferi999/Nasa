# Results v0: is the ground under Norilsk subsiding, and why?

NASA Space Apps 2026, *Dancing with the SARs*. Website: `docs/` (GitHub Pages).
Everything below is reproduced by `./run_v0.sh`; numbers come from
`outputs/zone_fits.csv`, `outputs/zone_comparisons.csv`, `data/ts/*_meta.json`.

## One-paragraph answer

**Thaw subsidence is real, reproducible, and controlled by the ground, not
by the industry.** Comparing every 30 m pixel with its own ~2 km
neighbourhood, the map of seasonal motion in summer 2020 (Sentinel-1A)
correlates with the one in 2025 (Sentinel-1C) at r = 0.36, against
0.007 ± 0.023 for chance (2025 map shifted by 6 km), i.e. ~15σ. The
atmosphere is independent between summers, so the repeating part is the ground.
Wetland ground (WorldCover 90, 48 km²) sinks 4.7 mm (2020) and 3.7 mm (2025)
more than its surroundings, following √TDD at every date (r = −0.99 and
−0.96): this is how a deepening thaw front behaves (Stefan's law). Industrial
sites do **not** sink more than their surroundings: the big plants
(Nadezhda smelter, TPP-3, ventilation shafts) sink *less*, in both summers.
At zone scale, industry is not detectably different from tundra (opposite
signs in the two summers). Whether the whole lowland is sinking relative to
bedrock is still not measurable: at 10-50 km scale the atmosphere dominates.

## What we measured

| | |
|---|---|
| Sentinel-1 | path 136 descending, 6 IW bursts, HyP3 ISCE2 multi-burst, 40 m, VV |
| Summer 2020 | 8 dates (15 Jun to 19 Sep), 13 pairs after QC, SBAS residual 2.0 mm |
| Summer 2025 | 7 dates (13 Jun to 5 Sep), 11 pairs after QC, SBAS residual 1.9 mm |
| Troposphere | per pair, phase vs height (+ planar ramp) fitted on coherent pixels at 250-1100 m only, subtracted everywhere |
| Coverage | 55 % (2020) and 46 % (2025) of the 30 m grid solved; incidence ~41° |
| NISAR 2026 | 4 GUNW pairs (12 days): asc 111 21 Jun–3 Jul, 20 Aug–1 Sep; desc 6 26 Jun–8 Jul, 25 Aug–6 Sep |
| Climate | ERA5-Land daily air temperature → cumulative thawing degree-days (TDD) |
| Zones | industrial (OSM, 69 km²), TPP-3 site (r = 1 km), city (WorldCover), tundra NW (64 km²), tundra SE (83 km²), plateau bedrock reference (3 patches, 3.9 km²) |

Model tested per zone and summer: `disp(t) = a + b·√TDD(t)` (Stefan-type thaw).

### Local, reproducible signal (step 14, the strongest result)

Per pixel, the seasonal motion is the √TDD fit over the summer; the local map
is that minus its 2 km Gaussian average.

| test | value |
|---|---|
| r(local 2020, local 2025), 23,095 samples 150 m apart | **0.363** |
| same with the 2025 map shifted 6 km (24 directions) | 0.007 ± 0.023 (max 0.052) |
| r of the raw (not high-passed) seasonal maps | 0.013 (atmosphere) |
| r of local map with local relief (DEM error check) | −0.09 / +0.14; partial r after removing relief 0.38 |
| pixels < −10 mm in both summers | 21,808 (19.6 km²) vs 6,699 expected if unrelated (3.3×) |

By land cover, end of summer, vs ~2 km surroundings (mm, ±1σ):

| class | 2020 | 2025 |
|---|---|---|
| herbaceous wetland (48 km²) | **−4.7 ± 0.25** | **−3.7 ± 0.24** |
| moss and lichen (213 km²) | +0.7 | +0.7 |
| grassland (115 km²) | −0.2 | +0.9 |
| tree cover | +2.7 | +2.4 |
| bare / sparse | +3.2 | +1.5 |
| built-up | +3.5 | +0.1 |

Per OSM site (inside minus a 0.3-1.5 km ring): 87 sites, site-level r between
summers 0.41. 10 sites have the same sign and > 2σ in both summers: 7 sink
*less* than their surroundings (Nadezhda plant +16.5 / +7.1 mm, TPP-3
+14.0 / +5.4, ventilation shafts VS-1,2,4 +12.8 / +9.2, Norilskgeologiya,
Kayerkan coal mine, two unnamed), 3 sink more (an unnamed site at 69.211 N
88.342 E −14.6 / −3.3, Mayak ventilation shaft, a quarry). Landfills and
tailings on average −2.8 mm in both summers. Full list:
`outputs/local_facilities.csv`; hotspot clusters: `outputs/local_hotspots.csv`.

Reading: piles, gravel pads and rock fill keep built ground stable; the
ground that gives way on thaw is the wet, ice-rich natural ground. Caveat:
C-band phase is also sensitive to soil moisture, which could add an apparent
signal of a few mm in wetlands; the monotonic √TDD shape over the whole
summer is what thaw predicts, but the two effects are not separated here.

### Paired zone comparisons (after tropospheric correction)

Slope of the per-date difference series against √TDD, in mm per √(°C·day).
Negative = the first zone sinks more per unit of thaw.

| summer | comparison | Δb | σ | z | end-of-summer Δ (mm) |
|---|---|---|---|---|---|
| 2020 | industrial − tundra mean | +0.53 | 0.23 | +2.3 | +13.7 |
| 2020 | TPP-3 site − tundra mean | +0.81 | 0.28 | +2.8 | +24.4 |
| 2020 | tundra NW − tundra SE | +1.38 | 0.76 | +1.8 | +40.3 |
| 2025 | industrial − tundra mean | −0.30 | 0.21 | −1.4 | −10.6 |
| 2025 | TPP-3 site − tundra mean | −0.19 | 0.31 | −0.6 | −6.3 |
| 2025 | tundra NW − tundra SE | −0.33 | 0.51 | −0.6 | −6.2 |

Industry sinks significantly *less* than tundra in 2020 and not
significantly more in 2025; the two tundra controls differ by up to 40 mm.
At 20-30 km separation, zone comparisons still carry turbulent atmosphere;
the local test above does not.

### Reference stability

Each plateau patch compared with the mean of the other two, maximum
deviation over the summer: 2020: 7.2 / 6.6 / 4.8 mm; 2025: 2.2 / 6.9 / 5.4 mm.
No patch moves systematically relative to the others, so all three are kept.
In late September the plateau (~650 m) decorrelates (early snow): pairs with
less than 50 % usable reference pixels were dropped (2 pairs in 2020, 4 in 2025).

### NISAR cross-check

After ionosphere and troposphere screens from the product, out-of-swath
masking and a planar deramp, the median |east-west| motion on flat ground is
7.2 mm (early) and 9.0 mm (late) per 12 days, comparable to the median
|vertical| of 8.9 and 7.5 mm. With single pairs this is at the noise level, so
**NISAR cannot yet confirm the vertical-motion assumption** used for
Sentinel-1. Track 6 descending covers only the south-west half of the AOI.

## What it supports

| Hypothesis | v0 status | Why |
|---|---|---|
| A. Climatic thaw (natural ground subsides with √TDD) | **local evidence** | wetlands sink vs surroundings along √TDD in both summers; absolute motion vs bedrock still atmosphere-dominated (σ 0.25-0.6) |
| B. Anthropogenic thaw (industry > tundra) | **not supported** | zone pairs: opposite signs; site by site, the large plants sink less than their surroundings in both summers |
| C. Stable ground | **not testable yet** | same reason as A: "not significant" here means noisy, not stable |

## What it cannot say

1. **No multi-year trend.** Snow decorrelates C-band every winter, and
   Sentinel-1 acquired nothing over Norilsk in 2022–2024. Summers are compared,
   not linked.
2. **The atmosphere dominates absolute motion.** All lowland zones move
   together by tens of mm from one pass to the next (inter-zone correlation
   0.44–1.0, mostly above 0.85); a full reversal on 19 Sep 2020 shows it is not ground. The
   reference is ~550 m above the lowland zones, so stratified water vapour maps into
   apparent lowland motion. Paired tests cancel it; absolute rates do not.
3. **Error bars are statistical only.** They come from the scatter of 7–8
   dates. Zone-local atmosphere (km scale) is not in them; small zones (TPP-3,
   3 km²) are the most affected. Our synthetic test shows the method
   recovers zone differences to ±0.05 mm/√(°C·day) for large zones and
   ±0.2 for the TPP-3 zone.
4. **Vertical is assumed** (LOS / cos 41°, slopes < 5° only). NISAR does not
   yet confirm this (above).
5. **One tank is one pixel.** At 40 m we see the ground around TPP-3, not
   the foundation of tank no. 5. The two literature coordinates for the tank
   could not be resolved with Sentinel-1 amplitude (snowmelt in spring,
   no clear disappearing target in summer); we use the OSM TPP-3 site.
6. **Two summers.** 2020 was exceptionally warm (MAAT −3.7 °C, TDD 1752),
   2025 colder (−7.4 °C, 1435). Two summers can show a zone difference, not a
   climate trend.

## Method checks

- End-to-end synthetic test (`./run_v0.sh --synthetic`): known per-zone slopes
  and an east motion in the industrial zone are recovered (zone differences
  within 0.05 for large zones; NISAR east 2.8 mm vs 3.0 true).
- Sign conventions checked against product documentation: HyP3 positive phase =
  motion away from the sensor; NISAR LOS unit vectors target-to-sensor (ENU);
  ISCE3 interferogram = reference × conj(secondary).
- Each HyP3 pair carries its own reference point; pairs are re-referenced to
  the plateau patches before inversion (without this the SBAS residual was
  11 mm instead of 3–6 mm; with the empirical troposphere it is 2 mm).

## Next steps (in order of impact)

1. **Weather-model troposphere** (ERA5 via PyAPS, or GACOS) on top of the
   empirical height correction, to make the absolute (zone-scale) signal usable.
   Ground-truth the persistent hotspots with Sentinel-2 / field knowledge.
2. Process the other five summers (2017–2019, 2021, 2026): 495 credits total.
3. Stack all 8 NISAR pairs per track to beat single-pair noise; repeat the
   east/vertical test.
4. Annual late-summer pairs on built-up and rock to link summers.
