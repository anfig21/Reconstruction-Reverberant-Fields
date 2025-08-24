# Reconstruction of Reverberant Sound Fields
_Source code (MATLAB) for reconstructing room sound fields over **large spatial domains** from **distributed microphone arrays**._

This repository accompanies the method by **Figueroa-Duran & Fernandez-Grande (JASA, 2025)**. The approach models the **direct sound + early reflections** with a wave-based expansion (localising apparent image sources) and reconstructs the **late reverberant field** with kernel methods assuming a **sinc-like spatial correlation** (random-wave-field theory). The result is accurate **interpolation and extrapolation** across metres using compact arrays.

---

## The repo in short
- **Large-aperture reconstruction** from sparse, compact arrays  
- **Early/late split** aligned with room-impulse-response structure  
- **Robust extrapolation** via statistically **synthesised pressure points**  
- Suitable for **6-DoF/navigable audio**, spatial analysis, and room-acoustics studies

---

## Repository structure
- `Main_NBI_Reconstruction.m` – entry point orchestrating the full pipeline  
- `dataAcquisitionNBI.m` – load/organise multichannel data & geometry  
- `dataHandling.m` – utilities for formatting, batching, and I/O  
- `dirDOA_SRP_PHAT.m` – DOA estimation (SRP-PHAT) for direct/early parts  
- `earlyDOA.m`, `earlyRange.m` – localise apparent origins (image-source cues)  
- `reconstructReflection.m` – wave-based reconstruction of direct/early field  
- `kernelReconstructionOverlap.m` – kernel reconstruction of late field + stitching  
- `windowRIR.m`, `peakDetection.m` – split RIR; detect salient arrivals  
- `clusterCoefficients.m` – coefficient handling for clustered sources  
- `plotFreqResponse.m`, `setupPlot.m` – plotting helpers  
- `toolbox/` – additional helpers

> Open each script’s header for parameter notes and expected inputs.

---

## Getting started
1. **Clone**
   ```bash
   git clone https://github.com/anfig21/Reconstruction-Reverberant-Fields

2. **Add to MATLAB path**
   ```bash
   addpath(genpath('Reconstruction-Reverberant-Fields'));

3. **Prepare inputs**
   * Multichannel recordings (RIRs or time-domain signals)
   * Microphone positions (N×3, metres); source positions if known
   * Sampling rate and any array groupings (for distributed arrays)

4. **Configure & run**
   * Edit parameters at the top of `Main_NBI_Reconstruction.m` (grid extent, frequency range, windowing, DOA settings, etc.).
   * Then run:
     ```bash
     Main_NBI_Reconstruction

5. **Outputs**
   * Reconstructed pressure field on the target grid/volume (time/frequency)
   * DOA and early-reflection estimates with optional visualisations

---

## Method overview (very brief)
- Early field (direct + early reflections): wave-based expansion with DOA/range cues to place apparent image sources.
- Late field: kernel reconstruction using a sinc-like spatial correlation; statistically synthesised points stabilise long-range extrapolation.

---

## Requirements
- MATLAB (tested by the authors)
- Standard Signal Processing/Phased Array toolboxes may be useful (see script headers)

---

## Data
The recordings used in the paper are available here:  
**➡︎ [Download the Niels Bohr Institute dataset (Zenodo, v1.2)](https://doi.org/10.5281/zenodo.14212820 "DOI link")**
**➡︎ [Download the DTU Clasroom dataset (data.dtu)](https://doi.org/10.11583/DTU.25867705.v1 "DOI link")**

You can also use your own multichannel recordings; see `dataAcquisitionNBI.m` / `dataHandling.m` for format hints.

---

## Citing
If you use this code or build upon the method, please cite:

A. Figueroa-Durán and E. Fernández-Grande, “Reconstruction of reverberant sound fields over large spatial domains,” Journal of the Acoustical Society of America, 157(1):180–190, 2025. https://doi.org/10.1121/10.0034833

## Contributing
Issues and pull requests are welcome. For questions about the research, please reference the paper above, include a minimal example and refer them to anfig@dtu.dk
