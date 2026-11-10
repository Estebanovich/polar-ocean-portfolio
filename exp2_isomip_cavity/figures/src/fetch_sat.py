"""Fetch the satellite images used in Figure 0 of Experiment 2 (all in Antarctic Polar Stereographic, EPSG:3031).

  antarctica_bluemarble.jpg        NASA Blue Marble: Next Generation, shaded relief and bathymetry
                                   (NASA Earth Observatory), whole continent; locator inset of (a).
  ase_modis_2025-12-01.jpg         MODIS Terra corrected reflectance (true colour), 1 Dec 2025, Amundsen Sea
                                   Embayment; panel (a).
  pig_landsat9_2025-12-01.jpg      Landsat 9 OLI-2, Collection 2 Level-2 surface reflectance (bands 4-3-2),
                                   scene LC09_L2SR_233113_20251201_02_T2 (path 233, row 113; 0.04 % cloud), Pine
                                   Island Glacier ice shelf and Pine Island Bay; panel (b). Drawn in a polar
                                   stereographic projection rotated so that the ice flows from left to right (PIG_CRS:
                                   EPSG:3031 with central meridian 112.5 E instead of 0), which also puts the scene
                                   edges along the image axes, so the crop holds no empty (no-data) pixels.

Sources and licences:
  - NASA GIBS / Worldview snapshot API (https://wvs.earthdata.nasa.gov). NASA imagery is free to use;
    credit "NASA Worldview / GIBS" (and NASA Earth Observatory for Blue Marble).
  - Landsat: U.S. Geological Survey, public domain ("Landsat imagery courtesy of the U.S. Geological Survey"),
    rendered from the Cloud-Optimised GeoTIFFs by the Microsoft Planetary Computer data API
    (https://planetarycomputer.microsoft.com), which reprojects to EPSG:3031. The stretch (rescale) is a plain
    linear stretch of the surface-reflectance integers, the same for the three bands.
The scene and the date were chosen as clear-sky summer images (the Landsat scene metadata give 0.04 % cloud;
the MODIS image of the same day was checked by eye). BBOX is xmin,ymin,xmax,ymax in metres (EPSG:3031, or PIG_CRS for the Landsat image).
Usage: fetch_sat.py   (writes into sat/; every file is a JPEG under 0.5 MB)
"""
import io
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parent / "sat"
GIBS = "https://wvs.earthdata.nasa.gov/api/v1/snapshot?REQUEST=GetSnapshot&CRS=EPSG:3031&FORMAT=image/jpeg"
PC = "https://planetarycomputer.microsoft.com/api/data/v1/item/bbox"

BBOX_ANT = (-3400000, -3400000, 3400000, 3400000)
BBOX_ASE = (-2050000, -860000, -1200000, -160000)       # 850 x 700 km
PIG_CRS = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=112.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
BBOX_PIG = (806000, -1459000, 986000, -1319000)          # 180 x 140 km, in PIG_CRS
IMAGES = {
    "antarctica_bluemarble.jpg": GIBS + "&TIME=2020-01-01&BBOX=%d,%d,%d,%d&WIDTH=480&HEIGHT=480"
                                        "&LAYERS=BlueMarble_ShadedRelief_Bathymetry" % BBOX_ANT,
    "ase_modis_2025-12-01.jpg": GIBS + "&TIME=2025-12-01&BBOX=%d,%d,%d,%d&WIDTH=1215&HEIGHT=1000"
                                       "&LAYERS=MODIS_Terra_CorrectedReflectance_TrueColor" % BBOX_ASE,
    "pig_landsat9_2025-12-01.jpg": PC + "/%d,%d,%d,%d/1000x778.png" % BBOX_PIG
                                       + "?collection=landsat-c2-l2&item=LC09_L2SR_233113_20251201_02_T2"
                                         "&assets=red&assets=green&assets=blue&rescale=7000,52000"
                                         "&coord_crs=%s&dst_crs=%s&nodata=0" % ((urllib.parse.quote(PIG_CRS),) * 2),
}

if __name__ == "__main__":
    from PIL import Image
    OUT.mkdir(exist_ok=True)
    for name, url in IMAGES.items():
        with urllib.request.urlopen(url, timeout=120) as r:
            im = Image.open(io.BytesIO(r.read()))
        if im.mode == "RGBA":     # Landsat: refuse a crop that reaches outside the scene (no-data on its border;
            import numpy as np    # a few interior no-data pixels are open water with zero reflectance, drawn black)
            a = np.asarray(im.getchannel("A"))
            assert min(a[0].min(), a[-1].min(), a[:, 0].min(), a[:, -1].min()) == 255, f"{name}: crop leaves the scene"
        im = im.convert("RGB")
        im.save(OUT / name, quality=88, optimize=True)
        print("wrote", OUT / name, im.size)
