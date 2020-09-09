from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, ICRS
from mwa_metadb_utils import getmeta



obsids = [1133775752,
        1137236608,
        1150234552,
        1164110416,
        1182630616,
        1194350120,
        1220886016,
        1220886016,
        1222697776,
        1225462936,
        1255444104,
        1275085816,
        1275092416,
        1275094456,
        1275172216,
        1275177136,
        1275178816,
        1275258616,
        1275431416,
        1275863416,
        1275866536,
        1276725752,
        1278106408]

earth_location = EarthLocation.from_geodetic(lon="116:40:14.93", lat="-26:42:11.95", height=377.8)
for obsid in obsids:
    beam_meta_data = getmeta(service='obs', params={'obs_id':obsid})
    alt = beam_meta_data[u'metadata'][u'elevation_pointing'] #in sexidecimal
    az = beam_meta_data[u'metadata'][u'azimuth_pointing'] 
    
    ra = beam_meta_data[u'metadata'][u'ra_pointing']
    dec = beam_meta_data[u'metadata'][u'dec_pointing']

    # change the time to half way though if after 2019
    if obsid > 1230339618:
        obsid = obsid + (beam_meta_data[u'stoptime'] - beam_meta_data[u'starttime']) / 2 
    obstime = Time(float(obsid),format='gps')
    
    # Convert meta Alt Az to ra and dec
    beam_atlaz = SkyCoord(az, alt, unit=(u.deg, u.deg), frame='altaz',
                          obstime=obstime, location=earth_location)
    tile_radec = beam_atlaz.transform_to(ICRS())
    conv_ra = tile_radec.ra.deg
    conv_dec = tile_radec.dec.deg

    # Convert meta ra and dec to Alt Az
    sky_posn = SkyCoord(ra, dec, unit=(u.deg,u.deg))
    altaz = sky_posn.transform_to(AltAz(obstime=obstime, location=earth_location))
    Alt = altaz.alt.deg
    Az  = altaz.az.deg 
    
    print(f'ObsID: {obsid}')
    print(f'meta ra:  {ra:6.2f}    dec: {dec:6.2f}')
    print(f'conv ra:  {conv_ra:6.2f}    dec: {conv_dec:6.2f}')
    print(f'meta alt: {alt:6.2f}    az:  {az:6.2f}')
    print(f'conv alt: {Alt:6.2f}    az:  {Az:6.2f}')
    print('')
