import numpy as np
def get_map_ticks(lat_ll, lon_ll, lat_tr, lon_tr):
    delta_lat = lat_tr - lat_ll
    delta_lon = lon_tr - lon_ll
    pot_lon = np.floor(np.log10(delta_lon))
    lon_tick = np.floor(delta_lon / 10 ** pot_lon) * 10 ** pot_lon / 4

    pot_lat = np.floor(np.log10(delta_lat))
    lat_tick = np.floor(delta_lat / 10 ** pot_lat) * 10 ** pot_lat / 3

    lon_tick_array = np.arange(
        lon_tick*int((lon_ll - delta_lon * 0.3)/lon_tick),
        lon_tick*int((lon_tr+delta_lon * 0.3)/lon_tick),
        lon_tick)

    lat_tick_array = np.arange(
        lat_tick*int((lat_ll - delta_lat * 0.3)/lat_tick),
        lat_tick*int((lat_tr+delta_lat * 0.3)/lat_tick),
        lat_tick)
    return (lat_tick_array, lon_tick_array)
