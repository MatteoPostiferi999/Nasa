# Results v0: is the ground under Norilsk subsiding, and why?

NASA Space Apps 2026, *Dancing with the SARs*. Website: `docs/` (GitHub Pages).
Everything below is reproduced by `./run_v0.sh`; numbers come from
`outputs/zone_fits.csv`, `outputs/zone_comparisons.csv`, `data/ts/*_meta.json`.

## One-paragraph answer

In the two summers we processed (2020, 2025), **industrial ground in Norilsk
did not subside detectably faster than natural tundra of the same elevation
and slope.** The industrial-minus-tundra difference changes sign between
summers (+9.9 mm in 2020, −16.5 mm in 2025, |z| < 1.7) and is no larger than
the difference between the two tundra controls themselves (34.8 mm in 2020).
This argues against a strong, city-wide *anthropogenic* thaw signal larger
than about 1–1.5 cm per summer (2σ). Whether the whole lowland, tundra
included, is sinking relative to bedrock (*climatic* thaw) **cannot be decided
yet**: the absolute signal is dominated by atmospheric delay, which we have
not removed in v0.

## What we measured

| | |
|---|---|
| Sentinel-1 | path 136 descending, 6 IW bursts, HyP3 ISCE2 multi-burst, 40 m, VV |
| Summer 2020 | 8 dates (15 Jun to 19 Sep), 13 pairs after QC, SBAS residual 3.3 mm |
| Summer 2025 | 7 dates (13 Jun to 5 Sep), 11 pairs after QC, SBAS residual 5.9 mm |
| Coverage | 55 % (2020) and 46 % (2025) of the 30 m grid solved; incidence ~41° |
| NISAR 2026 | 4 GUNW pairs (12 days): asc 111 21 Jun–3 Jul, 20 Aug–1 Sep; desc 6 26 Jun–8 Jul, 25 Aug–6 Sep |
| Climate | ERA5-Land daily air temperature → cumulative thawing degree-days (TDD) |
| Zones | industrial (OSM, 69 km²), TPP-3 site (r = 1 km), city (WorldCover), tundra NW (64 km²), tundra SE (83 km²), plateau bedrock reference (3 patches, 3.9 km²) |

Model tested per zone and summer: `disp(t) = a + b·√TDD(t)` (Stefan-type thaw).

### Paired comparisons (the robust test)

Slope of the per-date difference series against √TDD, in mm per √(°C·day).
Negative = the first zone sinks more per unit of thaw.

| summer | comparison | Δb | σ | z | end-of-summer Δ (mm) |
|---|---|---|---|---|---|
| 2020 | industrial − tundra mean | +0.35 | 0.24 | +1.5 | +9.9 |
| 2020 | TPP-3 site − tundra mean | +0.45 | 0.37 | +1.2 | +23.0 |
| 2020 | tundra NW − tundra SE | +0.23 | 0.95 | +0.2 | +34.8 |
| 2025 | industrial − tundra mean | −0.45 | 0.27 | −1.7 | −16.5 |
| 2025 | TPP-3 site − tundra mean | −0.42 | 0.35 | −1.2 | −16.9 |
| 2025 | tundra NW − tundra SE | +0.10 | 1.08 | +0.1 | −8.8 |

No comparison reaches 2σ in either summer, and the sign flips between
summers. The two controls disagree with each other as much as industry
disagrees with them.

### Reference stability

Each plateau patch compared with the mean of the other two, maximum
deviation over the summer: 2020: 7.0 / 4.4 / 4.3 mm; 2025: 2.9 / 3.8 / 3.0 mm.
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
| A. Climatic thaw (tundra subsides with √TDD) | **not testable yet** | absolute slopes have σ ≈ 0.5–1.1 mm/√(°C·day): atmosphere-dominated |
| B. Anthropogenic thaw (industry > tundra) | **not supported** | paired difference within 2σ in both summers, opposite signs |
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
  11 mm instead of 3–6 mm).

## Next steps (in order of impact)

1. **Tropospheric correction** for Sentinel-1 (ERA5 via PyAPS, or GACOS), then
   re-test A and C. This is the single change most likely to make the absolute
   signal usable.
2. Process the other five summers (2017–2019, 2021, 2026): 495 credits total.
3. Stack all 8 NISAR pairs per track to beat single-pair noise; repeat the
   east/vertical test.
4. Annual late-summer pairs on built-up and rock to link summers.
