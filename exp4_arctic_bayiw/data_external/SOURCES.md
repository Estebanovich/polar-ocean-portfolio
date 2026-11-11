# External data for Exp. 4

Files in this directory are not tracked by git; only this file is. Dates below are the download dates (2026-11-11). `shasum -a 256` of every file is given.

## 1. IBCAO v5.2 bathymetry (`ibcao/`)

Version obtained: **IBCAO v5.2** (released June 2026; grid dated 2026-02-16). Citation: IBCAO Version 5.2 Compilation Group (2026): The International Bathymetric Chart of the Arctic Ocean (IBCAO) Version 5.2. NERC British Oceanographic Data Centre, https://catalogue.ceda.ac.uk/uuid/30ad58b60d8c4a4dab6d64a40640540a (DOI 10.17043/ibcao-5.2; release paper for v5: Jakobsson et al. 2024, https://doi.org/10.1038/s41597-024-04278-w).

The full mosaic and the 16-tile product are 12 GB (> 2 GB), so they were not downloaded. The CEDA GeoTIFFs are uncompressed with one strip per row; `scripts/fetch_ibcao_window.py` downloads only the rows of a window by HTTP range requests and writes it as a GeoTIFF (same values as the originals, georeferencing copied from them). The window is 2050 x 2600 cells of 100 m (full-grid cols 31700-33749, rows 41100-43699; EPSG:3996, x = 261000 to 466000 m, y = -1200000 to -1460000 m), covering 76.5-78.5 N, 12-18 E. Re-running the script reproduces the three rasters bit for bit (checked with `cmp`). Land and the sea floor are in the same raster (positive land elevation, negative depth, metres, vertical datum mean sea level).

| File | Source URL | SHA-256 |
|---|---|---|
| `ibcao_v5_2_2026_depth_100m_svalbard.tif` (float32; depth, no Greenland ice sheet data) | https://dap.ceda.ac.uk/bodc/gebco/ibcao/ibcao_v5.2/no_greenland_ice_sheet_elevation_data/100mx100m_grid_cell_spacing/16_tiles/ibcao_v5_2_2026_depth_part_2_1.tiff (rows 12000-14549) and `..._part_2_0.tiff` (rows 0-49) | `b7d13733cb294e829aab2cfa80b850be35d1b54ca93b9f4f5b7c144c098ddd81` |
| `ibcao_v5_2_2026_sid_100m_svalbard.tif` (uint16; source identifier per cell) | https://dap.ceda.ac.uk/bodc/gebco/ibcao/ibcao_v5.2/100mx100m_data_source_identifier_grid/ibcao_v5_2_2026_sid_100m.tiff | `ce7c82be550424babc56b4d58117ae3c042b783a11cecfb27d6b887cfead0fb4` |
| `ibcao_v5_2_2026_tid_100m_svalbard.tif` (uint8; GEBCO type identifier per cell) | https://dap.ceda.ac.uk/bodc/gebco/ibcao/ibcao_v5.2/100mx100m_data_type_identifier_grid/ibcao_v5_2_2026_tid_100m.tiff | `df825deb3fff9544c9c314d2eeb9697bed32d402138c83b510e909c7dd73b943` |
| `ibcao-2026-source-data.csv` (table of source datasets; the `Number` column is the SID) | `ibcao-2026-source-data.csv` inside https://bolin.su.se/data/s3/upload-next/next-s3-uploads/662435dc-beae-4063-b799-7ba8544492d8/ibcao-5.2.zip (Bolin Centre, https://bolin.su.se/data/ibcao-5.2), extracted by range request from the 50.8 GB zip | `75d0af4dab4eb7825f6129c1878d5047dd2db91e1862b11d92a6e55001c35139` |

Sanity check of the georeferencing (point values): the outer Isfjorden mouth (78.1 N, 13.6 E) is -291 m with TID 11 (multibeam) and Bellsund-side land points are positive.

**Terms and redistribution.** The GEBCO Arctic page says that downloading the data means accepting its disclaimer (not for navigation or safety at sea; provided "as is"). CEDA licenses the data under the UK Open Government Licence v3 ("when using these data you must cite them correctly"); the Bolin Centre record gives the Open Data Commons Attribution License (ODC-By). Both allow redistribution with attribution, so IBCAO may be redistributed with the citation above, but this repository does not commit it (rule: `data_external/*` is ignored).

## 2. NPI coastline, S100 Kartdata (`npi_coast/`)

Norwegian Polar Institute, *Kartdata Svalbard 1:100 000 (S100 Kartdata)*, https://doi.org/10.21334/npolar.2014.645336c7, record https://data.npolar.no/dataset/645336c7-adfe-4d5a-978d-9426fe788ee3 (modified 2024-10-03; the product is updated several times a year, shapefile dated 2023-09/10). No registration needed. ETRS89 / UTM 33N (EPSG:25833).

| File | Source URL | SHA-256 |
|---|---|---|
| `NP_S100_SHP.zip` (290 MB, 92 files) | https://public.data.npolar.no/kartdata/NP_S100_SHP.zip | `61030b67371be1e46ed9e81c744a6132980dbb382bd7eaf871f85c6f7af95260` |
| `extracted/NP_S100_SHP/S100_Land_f.shp` (land polygons) | extracted from the zip | `8a6518be1a9610e8e0d7c06198027053dcf5fe91d4a09903d5f0d01ea327b799` |
| `extracted/NP_S100_SHP/S100_Land_f.dbf` | same | `21931a36457147b4758ccb5a8f1c3ce15f5726b88cb1dc81b2465720281d7ecc` |
| `extracted/NP_S100_SHP/S100_Land_f.shx` | same | `e7e7ac306aa7de7991cdd472500c64c90f7183a342386e6c8d1ab9291f61884f` |
| `extracted/NP_S100_SHP/S100_Land_f.prj` (also `S100_Land_l.prj`) | same | `e65678526a7a86d34803b06c1b2412a74ea6df15b446507b8d2c89548af78c9d` |
| `extracted/NP_S100_SHP/S100_Land_f.cpg` (also `S100_Land_l.cpg`) | same | `aa76969c68d31359f63100513fe53e3c1a3509e8c053804b19041e0fe7caf009` |
| `extracted/NP_S100_SHP/S100_Land_l.shp` (land outline lines) | same | `6b58109446c6cb50e6e4abc5e07e10b569244c6b921b86dce3460d206fcf3bfb` |
| `extracted/NP_S100_SHP/S100_Land_l.dbf` | same | `065c181e2ac213901cc2e72b9e44b386073bd2de8fe7eeffe7673d4fe37f24d5` |
| `extracted/NP_S100_SHP/S100_Land_l.shx` | same | `ffd3cdd03bee98b753047f53e6ac8584f4c08ab2584885f6918b050885944a32` |
| `extracted/NP_S100_SHP/_Lisens_License.txt` | same | `e2c527d86e24ffe95be343285f049020098b9e89701021455db1f379b400a7e4` |

**Licence:** Creative Commons Attribution 4.0 International. The map data must be credited ("Norwegian Polar Institute"), parts are of older origin and not for navigation, and `_Lisens_License.txt` must accompany any copy or derived data.

## 3. UNIS Hydrographic Database (`unis_hd/`)

UNIS HD, University Centre in Svalbard (Skogseth, Ellingsen et al.), https://doi.org/10.21334/unis-hydrography, record https://data.npolar.no/dataset/39d9f0f9-af12-420c-a879-10990df2e22d (modified 2025-11-07); see also Skogseth et al. 2020, https://doi.org/10.1016/j.pocean.2020.102394. Licence: **CC BY 4.0** (the `license` attribute of the files); the contributor citations are in `citationlist.txt`.

The database holds 82 352 single-profile NetCDF files (0-34 E, 75-83 N). Downloaded: the files whose name coordinates are within **77.5-78.0 N and 13.5-17.5 E**, all months and years: **1547 profiles, 1890-2019, 89 MB** (August 639, September 487, July 109; March-April 63). Source: `https://api.npolar.no/dataset/39d9f0f9-af12-420c-a879-10990df2e22d/attachment/<attachment id>/_blob`, one file per profile named `UNIS_HD_<YYYYMMDD>T<HHMMSS>_<lat>_<lon>.nc` (variables PRES, TEMP, PSAL, CNDC, TIME, LATITUDE, LONGITUDE). Every file was checked against the sha256 that the API lists for it (0 failures). The whole-database NetCDF files (`UNIS_HD.nc`, `UNIS_HD_1.nc`) and the decade zips were not downloaded.

| File | SHA-256 |
|---|---|
| `profiles/*.nc` (1547 files) | per-file hashes are those listed by the API; hash of the sorted `shasum -a 256 *.nc` list: `6f10f341312c3f79d85520050d10288bf512a73020d7b4a5a4cc0ddcc75c44f4` |
| `citationlist.txt` | `a041b864f02ed8ec5883a2020cf44da282ce279605787fea29ac28d2c1b1b86e` |

The box is a coarse bounding box and includes parts of Bellsund and the shelf; the Van Mijenfjorden polygon, years and minimum number of profiles are chosen in PD0.

## 4. Satellite images (`sat/`)

Made by `figures/src/fetch_sat.py` (STAC search and the Microsoft Planetary Computer data API, https://planetarycomputer.microsoft.com). Sentinel-2 L2A true colour (B04, B03, B02), mosaic of the granules T33XVG and T33XWG, EPSG:32633, box 475839-560386 E, 8602776-8659607 N (84.5 x 56.8 km), 1800 x 1000 px. Credit: **Contains modified Copernicus Sentinel data (2023, 2024), processed by ESA**; free, full and open access (Copernicus Sentinel data terms).

| File | Scene | Downloaded | SHA-256 |
|---|---|---|---|
| `vmf_s2_2023-07-05.jpg` | 5 Jul 2023, S2A and S2B passes (tile cloud cover 0.04 %), summer without ice | 2026-11-11 | `67454e695e6e4cc94377b2761fc90f8132a7340f83d503e2c09d51099db71ab8` |
| `vmf_s2_2024-03-25.jpg` | 25 Mar 2024, S2A and S2B passes (cloud 5 %), landfast ice in the inner fjord, open water off Bellsund | 2026-11-11 | `4c8f3d3a28b81ecd3a0e0c49ccf1d4be68d4ab50e699902b2494b09c68228c68` |

## Pending for Esteban

Nothing was blocked: every source downloaded without login, key or licence acceptance. Two caveats: (1) only a window of IBCAO was downloaded (the full 100 m grid is 12.6 GB and was deliberately skipped); if the whole mosaic is wanted, download `ibcao_v5_2_2026_depth_100m.zip` from the CEDA URL above. (2) the IBCAO v5.2 release paper is the v5.0 one (no v5.2 paper exists yet).
