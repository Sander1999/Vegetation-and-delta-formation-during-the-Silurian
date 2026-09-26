# Sources and parameter provenance

Accessed 26 September 2026. Parameters selected for convenience or sensitivity remain scenario values even when their governing mechanism is supported by a publication.

## Supplied experimental account

Bortolotti, M. H. and Bert da Costa, S. (2022). *Effects of the Silurian greening on delta structure and delta formation*. TU Delft AESB2230 essay, dated 19 December 2022.

The supplied DOCX is the source for discharge, tilt, experimental durations, vegetation placement and manual interventions. Its measurements are primarily photographs and visual observations. It does not contain the quantitative geometry, flux and vegetation data necessary to calibrate a unique numerical reconstruction. The project uses its experimental question without assuming that every causal interpretation or palaeobotanical statement is established.

## Apparatus and sediment

- [Emriver: EM4 specifications](https://emriver.com/em4/). Current manufacturer dimensions are 3.7 × 1.2 m. These constrain a nominal apparatus-sized domain, not the exact exposed sediment area in 2022.
- [Emriver: modelling media specification sheet](https://emriver.com/wp-content/uploads/2026/02/Media-Data-Sheet-1.pdf). Ground melamine plastic has reported density 1.6 g/cm³. Current Carbondale mean diameters: yellow 1.4 mm, white 1.0 mm, black 0.7 mm and red 0.4 mm. Memphis means differ: 1.02, 0.62, 0.33 and 0.14 mm respectively. Neither exact historical mixture nor colour fractions are supplied in the essay. Do not treat the current specification as a measurement of that experiment.

## Mechanisms and comparison design

- Lauzon, R. and Murray, A. B. (2022). [Discharge Determines Avulsion Regime in Model Experiments With Vegetated and Unvegetated Deltas](https://doi.org/10.1029/2021JF006225). *Journal of Geophysical Research: Earth Surface*, 127, e2021JF006225. This is the published article associated with the vegetation/discharge preprint cited in the essay. It uses DeltaRCM Vegetation, rather than Delft3D FM or an EM4 flume. It supports examining discharge, sediment supply, channel reoccupation and network persistence separately. Its modern emergent-vegetation rules and field-scale parameters are not a calibration for this flume-scale setup.
- Lauzon, R. and Murray, A. B. (2018). [Comparing the Cohesive Effects of Mud and Vegetation on Delta Evolution](https://doi.org/10.1029/2018GL079405). *Geophysical Research Letters*, 45, 10437–10445. Vegetation is represented through reduced lateral sediment transport and increased hydraulic resistance. This supports separating reinforcement from drag and avoiding the claim that a roughness-only treatment resolves root mechanics.
- Nardin, W. and Edmonds, D. A. (2014). [Optimum vegetation height and density for inorganic sedimentation in deltaic marshes](https://doi.org/10.1038/ngeo2233). *Nature Geoscience*, 7, 722–726. The study finds competing hydraulic effects: intermediate vegetation can favour deposition, whereas tall or dense vegetation can keep sand in channels. It motivates resistance sensitivity and measurement of flow diversion, not a prescribed sign for the response.

## Early land plants

- Boyce, C. K. (2008). [How green was Cooksonia? The importance of size in understanding the early evolution of physiology in the vascular plant lineage](https://doi.org/10.1666/0094-8373%282008%29034%5B0179%3AHGWCTI%5D2.0.CO%3B2). *Paleobiology*, 34, 179–194. Fragmentary preservation and variation in physiological capability constrain how confidently a living plant can be reconstructed. This paper does not supply calibrated hydraulic drag or bank reinforcement for the present model. Modern tree-root or marsh properties cannot be assigned to Silurian plants without additional assumptions.

## Numerical framework

- [Deltares: D-Flow Flexible Mesh](https://www.deltares.nl/en/software-and-data/products/delft3d-fm-suite/modules/d-flow-flexible-mesh) and [D-Morphology](https://www.deltares.nl/en/software-and-data/products/delft3d-fm-suite/modules/d-morphology) describe the hydrodynamic and morphological components.
- [Official Delft3D FM manuals](https://content.oss.deltares.nl/delft3dfm2d3d/) are the authority for input keywords, sediment formulations and supported output fields. Verify settings against the installed engine version. General shallow-water, suspended-transport and Exner balances in `EXPERIMENT.md` describe the physical framework; they do not by themselves establish which optional closures a particular engine build executes.

The retained run configuration and verification record should identify the exact engine version, selected transport formula, grain properties, bed porosity, roughness or vegetation implementation, numerical settings and any unsupported options. A source citation supplies a rationale; only inspected configuration and output establish implementation and execution.
