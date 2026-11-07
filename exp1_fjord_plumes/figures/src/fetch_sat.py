"""Fetch the two satellite images used in Figure 0 from NASA GIBS (Worldview snapshot API).

Imagery: NASA EOSDIS GIBS / Worldview, free to use with attribution.
  sermilik_modis_2022-08-05.jpg  MODIS Terra corrected reflectance (true colour) + 15 m coastlines,
                                 Sermilik Fjord and Helheim Glacier, SE Greenland
  helheim_hls_2022-07-30.jpg     Harmonized Landsat Sentinel-2 (HLS S30, 30 m), Helheim Glacier terminus.
                                 Contains modified Copernicus Sentinel data 2022.
BBOX is lat_min,lon_min,lat_max,lon_max (EPSG:4326). Usage: fetch_sat.py  (writes into sat/)
"""
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parent / "sat"
API = "https://wvs.earthdata.nasa.gov/api/v1/snapshot?REQUEST=GetSnapshot&CRS=EPSG:4326&FORMAT=image/jpeg"
IMAGES = {
    "sermilik_modis_2022-08-05.jpg": dict(TIME="2022-08-05", BBOX="65.45,-39.0,66.55,-36.9", WIDTH=900, HEIGHT=1100,
                                          LAYERS="MODIS_Terra_CorrectedReflectance_TrueColor,Coastlines_15m"),
    "helheim_hls_2022-07-30.jpg": dict(TIME="2022-07-30", BBOX="66.28,-38.55,66.46,-37.75", WIDTH=1000, HEIGHT=560,
                                       LAYERS="HLS_S30_Nadir_BRDF_Adjusted_Reflectance"),
}

if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    for name, q in IMAGES.items():
        url = API + "".join(f"&{k}={v}" for k, v in q.items())
        urllib.request.urlretrieve(url, OUT / name)
        print("wrote", OUT / name)
