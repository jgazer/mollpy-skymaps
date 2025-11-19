 IBEX ENA Data
===============

Please read this information if you are interested in using the Mollpy-SkyMaps routines to produce ENA all-sky maps from IBEX data.

The Interstellar Boundary Explorer (IBEX) Mission's publicly released data include full-sky maps of Energetic Neutral Atoms (ENA) observations.
As such, these data are not included in this repository; instead, they can be obtained from the official IBEX Data Releases.
For more information about IBEX, including references and proper citation of used data, we refer to the IBEX website.

Instructions:

The IBEX Data Releases are available for download through the IBEX website at:	https://ibex.princeton.edu/DataRelease
The most recent IBEX global ENA maps are in *Data Release 17* (DR17) for IBEX-Lo, and *Data Release 18* (DR18) for IBEX-Hi.
For smooth usage with Mollpy-SkyMaps, download the data release you are interested in, unzip the folder, and move its content to this repository's folder:
	~/Data/IBEX/DR_17/  or  ~/Data/IBEX/DR_18/.
Alternatively, after downloading and unzipping the data, update the path variables in the "config.json" file:
	"path_ibex_lo" and "path_ibex_hi", to the directory in your local drive where the data is saved.

2025-10-31, J. Gasser
