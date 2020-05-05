4/14/2020 - Jacob LaFontaine (jlafonta@usgs.gov)

# Information
Processing of datasets to produce segment-based parameters for the USGS National Hydrologic Model with Geospatial Fabric version 1.1 (GF_v11)

- Parameters computed:
	- `seg_width`
	- `seg_depth`
	- `mann_n`

- Base datasets:
	- Existing PRMS parameters from the NHM paramdb version 1.1:
		- `seg_slope`

- Existing datasets:
	- GF_v11 stream segment shapefile
	- River network shapefiles for North America and Central America from the HydroSHEDS HydroRIVERS website (https://www.hydrosheds.org/page/hydrorivers)

# Requirements
The following external code and datasets are required to complete the workflow below.
- R (https://www.r-project.org)
- River network shapefiles for North America and Central America from the HydroSHEDS HydroRIVERS website (https://www.hydrosheds.org/page/hydrorivers)

# Workflow
Processing:
1. HydroRIVERS shapefiles of North and Central America were merged into one shapefile in ArcGIS.
2. The merged shapefile of HydroRIVERS was converted to a point coverage in ArcGIS.
3. A spatial join was done using the nearest neighbor method to transfer the HydroRIVERS attributes to the GF_v11 stream segments.
4. The new attribute table for GF_v11 stream segments was output to a table.
5. Code was written in R to read the GF_v11 attribute table and output attributes of stream width and depth to .csv files for inclusion in the PRMS parameter database.
6. Manning's n was computed using the `seg_slope` parameter from the PRMS parameter database. The equation for computing Manning's n from segment slope was obtained from Bray (1979). That equation is n = 0.1 * seg_slope ** 0.18
7. The parameter `mann_n` was output to a .csv file for inclusion in the PRMS parameter database.

The R script `create_segment_parameters.R` was used to generate PRMS parameters `seg_depth`, `seg_width`, and `mann_n`.

