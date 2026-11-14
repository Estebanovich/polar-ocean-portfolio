"""Fetch the satellite images used in Figure 0 of Experiment 4 (writes into sat/; every file is a JPEG under 1 MB).

  vmf_s2_2024-03-25_domain.jpg   Sentinel-2 L2A true colour (B04, B03, B02), 25 Mar 2024, UTM 33N (EPSG:32633),
                                 the model domain with a margin (BBOX_DOMAIN): landfast ice in the inner fjord,
                                 open water in Bellsund; panel (a).
  svalbard_modis_2023-07-05.jpg  MODIS Terra corrected reflectance (true colour), 5 Jul 2023, a clear summer day,
                                 Svalbard (BBOX_SVALBARD, lat/lon, EPSG:4326, the width scaled by cos 78.6 deg so
                                 that the image is close to conformal at Svalbard); locator inset of (a).

The Sentinel-2 image is a mosaic of the granules in SCENES (tiles T33XVG and T33XWG meet inside the fjord).
Scene selection (Task 1b): STAC search (Planetary Computer, collection sentinel-2-l2a, bbox 14-17.5 E,
77.5-78 N, tile cloud cover < 25 %), keeping dates on which the mosaic covers the whole box; chosen by eye
among the cleanest dates with landfast ice. A thin seam between two passes crosses the south-east corner.

Sources and licences:
  - Contains modified Copernicus Sentinel data (2024), processed by ESA. Free, full and open access under the
    Copernicus Sentinel data terms. Rendered by the Microsoft Planetary Computer data API
    (https://planetarycomputer.microsoft.com), reprojected to EPSG:32633 with a plain linear stretch of the
    surface-reflectance integers (the same for the three bands).
  - NASA GIBS / Worldview snapshot API (https://wvs.earthdata.nasa.gov). NASA imagery is free to use; credit
    "NASA Worldview / GIBS".
BBOX_DOMAIN is xmin, ymin, xmax, ymax in metres (EPSG:32633): the model domain (cases.yaml bbox_utm, 455-555 km
E, 8594-8654 km N) with a margin of 4 km. BBOX_SVALBARD is lat_min, lon_min, lat_max, lon_max.
Usage: fetch_sat.py
"""
import io
import math
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parent / "sat"
PC = "https://planetarycomputer.microsoft.com/api/data/v1/item/bbox"
GIBS = "https://wvs.earthdata.nasa.gov/api/v1/snapshot?REQUEST=GetSnapshot&CRS=EPSG:4326&FORMAT=image/jpeg"
BBOX_DOMAIN = (451000, 8590000, 559000, 8658000)        # 108 x 68 km
SIZE_DOMAIN = (1404, 884)                               # 76.9 m per pixel
BBOX_SVALBARD = (76.3, 9.0, 80.9, 30.0)
LAT_SCALE = 78.6
SCENES = {
    "vmf_s2_2024-03-25_domain.jpg": (9000, [
        "S2B_MSIL2A_20240325T122649_R052_T33XWG_20240325T164435",
        "S2B_MSIL2A_20240325T122649_R052_T33XVG_20240325T163818",
        "S2A_MSIL2A_20240325T113721_R123_T33XWG_20240325T160834"]),
}


def render(item, rescale):
    from PIL import Image
    url = (PC + "/%d,%d,%d,%d/%dx%d.png?collection=sentinel-2-l2a&item=%s&assets=B04&assets=B03&assets=B02"
           "&rescale=0,%d&coord_crs=EPSG:32633&dst_crs=EPSG:32633&nodata=0"
           % (BBOX_DOMAIN + SIZE_DOMAIN + (item, rescale)))
    with urllib.request.urlopen(url, timeout=180) as r:
        return Image.open(io.BytesIO(r.read())).convert("RGBA")


def modis_size(width=800):
    la0, lo0, la1, lo1 = BBOX_SVALBARD
    return width, int(width * (la1 - la0) / ((lo1 - lo0) * math.cos(math.radians(LAT_SCALE))))


if __name__ == "__main__":
    import numpy as np
    from PIL import Image
    OUT.mkdir(exist_ok=True)
    for name, (rescale, items) in SCENES.items():
        mosaic = Image.new("RGBA", SIZE_DOMAIN, (0, 0, 0, 0))
        for it in items:
            mosaic = Image.alpha_composite(mosaic, render(it, rescale))
        a = np.asarray(mosaic.getchannel("A"))
        assert (a > 0).mean() > 0.99, f"{name}: mosaic leaves part of the box empty"
        im = mosaic.convert("RGB")
        im.save(OUT / name, quality=84, optimize=True)
        print("wrote", OUT / name, im.size)
    w, h = modis_size()
    url = GIBS + ("&TIME=2023-07-05&BBOX=%g,%g,%g,%g" % BBOX_SVALBARD
                  + f"&WIDTH={w}&HEIGHT={h}&LAYERS=MODIS_Terra_CorrectedReflectance_TrueColor")
    with urllib.request.urlopen(url, timeout=120) as r:
        im = Image.open(io.BytesIO(r.read())).convert("RGB")
    im.save(OUT / "svalbard_modis_2023-07-05.jpg", quality=86, optimize=True)
    print("wrote", OUT / "svalbard_modis_2023-07-05.jpg", im.size)
