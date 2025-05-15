# TSCF: A statistic-based tool for smoothing temperature sensitivity effects on matric potential dielectric readings

The Temperature Sensitivity Correction Function (TSCF) is introduced as straightforward method to localise, quantify, and smooth the daily-scale periodic noise. This function is based on a dynamic moving average approach using a combination of two centred simple moving average (SMA) with different window size (k). First, a fast SMA (mp_f) with a window size (k1), in this example equal to a hourly period, is used to decrease data noise normally produced in this type of high resolution observations. This magnitude allows for keeping short reactions times produced by fast wetting fronts. The size of k1 can also be set to a different resolution depending on the data resolution and noise intensity (e.g., 3 hours). Secondly, a slow SMA (mp_s) with a window equal to a daily period (k2) is used to smooth the daily oscillation produced by temperature sensitivity effects. The resulting corrected time series is a dynamic combination of fast and slow smoothing routines.

### How does it work?

To produce an algorithm able to recognise and implement a proper combination of these time series, we subtracted the mp_f from the mp_s to obtain a convergence-divergence noise signal (CDNS). This application resembles to the Moving Average Convergence-Divergence (MACD) oscillator proposed by Appel (2005) for economic sciences, with the simplification of using the simple unweighted mean (Savitzky & Golay, 1964) instead of the exponential moving average. The CDNS function allows to dynamically localise and quantify the intensity of the noise signal.

We avoid smoothing important wetting processes in over the dry season. When a wetting front reaches the matric potential sensor under lower water content (more negative potentials), the observed potential will drop towards zero and will change the noise periodical pattern producing an easily-noticeable peak in the CDNS spectrum. The capacity to recognise these peaks is best when using non-linear magnitudes of matric potential (e.g., kPa). This enables us to set the two conditions used for combining both mp_f and mp_s time series into a corrected series of data over an iterative process: i) the matric potential threshold (Tau, τ) at which the temperature sensitivity effects become significant (or visible), and ii) the change in the noise pattern when a wetting process is perceived.

The successful dynamic combination of h_f and h_s ensures that fast wetting fronts are no smoothed and therefore described in high resolution.

### Authors:

[Sebastián Bravo Peña (1,2) ![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0009-0005-3970-6187)\
[Meindert Commelin (1) ![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0001-7460-1915)\
[Ole Wendroth (3)](https://scholar.google.com/citations?user=cf-NkNEAAAAJ&hl=en)

1: Soil Physics and Land Management group; Wageningen University, Wageningen, The Netherlands\
2: Instituto de Ingeniería Agraria y Suelos, Universidad Austral de Chile, Valdivia, Chile\
3: Department of Plant and Soil Sciences, University of Kentucky, Lexington, United States\
Contact: [sebastian.bravopena\@wur.nl](mailto:sebastian.bravopena@wur.nl){.email}

## Performing calculations
