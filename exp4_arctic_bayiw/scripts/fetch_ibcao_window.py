"""Extract the Svalbard window of IBCAO v5.2 (depth, SID, TID; EPSG:3996, 100 m) with HTTP range requests.

The CEDA GeoTIFFs are uncompressed, one strip per row, so only the rows of the window are downloaded
(full-grid pixel window cols 31700-33749, rows 41100-43699, covering 76.5-78.5 N, 12-18 E). Depth comes from the
16-tile product (tiles 2_1 and 2_0); SID and TID from the single 100 m grids. Source: CEDA archive
https://dap.ceda.ac.uk/bodc/gebco/ibcao/ibcao_v5.2/ (licence: UK Open Government Licence v3; also ODC-By at Bolin).
Usage: fetch_ibcao_window.py OUTDIR
"""
def parse(d):
    bo='<' if d[:2]==b'II' else '>'
    ver=struct.unpack(bo+'H',d[2:4])[0]
    tags={}
    TS={1:1,2:1,3:2,4:4,5:8,6:1,7:1,8:2,9:4,10:8,11:4,12:8,16:8,17:8,18:8}
    FM={1:'B',2:'c',3:'H',4:'I',5:'II',6:'b',7:'B',8:'h',9:'i',11:'f',12:'d',16:'Q',17:'q',18:'Q'}
    if ver==42:
        off=struct.unpack(bo+'I',d[4:8])[0]; n=struct.unpack(bo+'H',d[off:off+2])[0]; ents=[d[off+2+12*i:off+14+12*i] for i in range(n)]; esz=4; hdr='HHI'
    else:
        off=struct.unpack(bo+'Q',d[8:16])[0]; n=struct.unpack(bo+'Q',d[off:off+8])[0]; ents=[d[off+8+20*i:off+28+20*i] for i in range(n)]; esz=8
    for e in ents:
        if ver==42: tag,typ,cnt=struct.unpack(bo+'HHI',e[:8]); val=e[8:12]
        else: tag,typ,cnt=struct.unpack(bo+'HHQ',e[:12]); val=e[12:20]
        size=TS.get(typ,1)*cnt
        if size<=esz: raw=val[:size]
        else:
            o=struct.unpack(bo+('I' if ver==42 else 'Q'),val[:esz])[0]
            raw=d[o:o+size] if o+size<=len(d) else None
        if raw is None: tags[tag]=('offset-beyond',typ,cnt,val); continue
        if typ==2: tags[tag]=raw
        elif typ in (5,): tags[tag]=struct.unpack(bo+'%dI'%(2*cnt),raw)
        else: tags[tag]=struct.unpack(bo+'%d%s'%(cnt,FM[typ]),raw)
    return tags

import struct, sys, requests, numpy as np
from concurrent.futures import ThreadPoolExecutor
CE="https://dap.ceda.ac.uk/bodc/gebco/ibcao/ibcao_v5.2"
C0,C1,R0,R1=31700,33750,41100,43700   # full-grid pixel window [C0,C1) x [R0,R1)
S=requests.Session()
def get(url,a,b):
    for k in range(5):
        try:
            r=S.get(url,headers={'Range':f'bytes={a}-{b}'},timeout=120)
            if r.status_code==206: return r.content
        except Exception as e: pass
    raise RuntimeError((url,a,b))
def window(url,rows_origin,cols_origin,rows,dtype):
    head=get(url,0,2*1024*1024-1)
    t=parse(head); off=t[273]; assert off[0]!='offset-beyond'
    assert t[259]==(1,) and t[278]==(1,)
    bps=t[258][0]//8; width=t[256][0]
    ca=C0-cols_origin; cb=C1-cols_origin
    def row(r):
        o=off[r]+ca*bps
        return np.frombuffer(get(url,o,o+(cb-ca)*bps-1),dtype=dtype)
    with ThreadPoolExecutor(12) as ex: return np.array(list(ex.map(row,rows))), t
def write_tif(path,arr,x0,y0,geokeys,geoascii,fmt):
    h,w=arr.shape; bps=arr.dtype.itemsize
    data=arr.astype('<'+arr.dtype.char).tobytes()
    ents=[]; extra=b''
    base=8+2+14*12+4  # header + count + entries + next ifd
    def add(tag,typ,cnt,val): ents.append((tag,typ,cnt,val))
    ext=[]
    def ofs(b):
        nonlocal extra
        o=base+len(extra); extra+=b+(b'\0' if len(b)%2 else b''); return o
    ps=struct.pack('<3d',100.,100.,0.); tp=struct.pack('<6d',0,0,0,x0,y0,0); gk=struct.pack('<%dH'%len(geokeys),*geokeys); ga=geoascii
    datapos_placeholder=0
    add(256,4,1,w); add(257,4,1,h); add(258,3,1,bps*8); add(259,3,1,1); add(262,3,1,1)
    add(273,4,1,'DATA'); add(277,3,1,1); add(278,4,1,h); add(279,4,1,len(data)); add(284,3,1,1); add(339,3,1,fmt)
    add(33550,12,3,ofs(ps)); add(33922,12,6,ofs(tp)); add(34735,3,len(geokeys),ofs(gk)); add(34737,2,len(ga),ofs(ga))
    assert len(ents)==15 or True
    n=len(ents); base2=8+2+n*12+4
    assert base2==base+ (n-14)*12, (n,)
    # recompute offsets with right base
    shift=base2-base
    out=bytearray(b'II*\0'+struct.pack('<I',8)+struct.pack('<H',n))
    datapos=base2+len(extra)
    for tag,typ,cnt,val in ents:
        if val=='DATA': val=datapos
        elif tag in (33550,33922,34735,34737): val+=shift
        if typ==3 and cnt==1: out+=struct.pack('<HHIHH',tag,typ,cnt,val,0)
        else: out+=struct.pack('<HHII',tag,typ,cnt,val)
    out+=struct.pack('<I',0)+extra+data
    open(path,'wb').write(out)
if __name__=='__main__':
    out=sys.argv[1]
    X0=-2910000+C0*100.; Y0=2910000-R0*100.
    jobs={}
    for name,url,dt,fmt in [('sid',CE+'/100mx100m_data_source_identifier_grid/ibcao_v5_2_2026_sid_100m.tiff','<u2',1),
                            ('tid',CE+'/100mx100m_data_type_identifier_grid/ibcao_v5_2_2026_tid_100m.tiff','u1',1)]:
        arr,t=window(url,0,0,range(R0,R1),dt)
        print(name,arr.shape,arr.dtype,arr.min(),arr.max())
        ga=t[34737]; gk=t[34735]
        write_tif(f'{out}/ibcao_v5_2_2026_{name}_100m_svalbard.tif',arr,X0,Y0,gk,ga,fmt)
    T=CE+'/no_greenland_ice_sheet_elevation_data/100mx100m_grid_cell_spacing/16_tiles/ibcao_v5_2_2026_depth_part_%d_%d.tiff'
    a,t=window(T%(2,1),29100,29100,range(R0-29100,43650-29100),'<f4')   # tile col 29100, row0 29100
    b,_=window(T%(2,0),43650,29100,range(0,R1-43650),'<f4')
    arr=np.vstack([a,b]); print('depth',arr.shape,np.nanmin(arr),np.nanmax(arr))
    write_tif(f'{out}/ibcao_v5_2_2026_depth_100m_svalbard.tif',arr,X0,Y0,t[34735],t[34737],3)
