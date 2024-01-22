def _check_requirements():
    BASEMAP_AVAILABLE = True
    CV2_AVAILABLE = True
    NETCDF_AVAILABLE = True
    try:
        from mpl_toolkits.basemap import Basemap
    except ModuleNotFoundError:
        BASEMAP_AVAILABLE = False

    try:
        from cv2 import pyrUp
    except ModuleNotFoundError: # pragma: no cover
        CV2_AVAILABLE = False
    try:
        from netCDF4 import Dataset
    except ModuleNotFoundError: # pragma: no cover
        NETCDF_AVAILABLE = False

    return (BASEMAP_AVAILABLE,
            CV2_AVAILABLE,
            NETCDF_AVAILABLE)


def _init_local_topodir():
    try:
        import os
        home = os.path.expanduser('~')
        LOCAL_TOPO_DIR = os.path.join(home, '.geonum')
        if not os.path.exists(LOCAL_TOPO_DIR): # pragma: no cover
            os.mkdir(LOCAL_TOPO_DIR)
        TOPO_INFO_FILE = os.path.join(LOCAL_TOPO_DIR,  "LOCAL_TOPO_PATHS")
        if not os.path.exists(TOPO_INFO_FILE): # pragma: no cover
            with open(TOPO_INFO_FILE,'w') as f:
                f.write(f'{LOCAL_TOPO_DIR}\n')
    except Exception: # pragma: no cover
        print('Failed to create local topo directory for geonum')
        LOCAL_TOPO_DIR, TOPO_INFO_FILE = None, None
    return (LOCAL_TOPO_DIR, TOPO_INFO_FILE)


def _init_dir_and_version():
    import os
    path = os.path.abspath(os.path.dirname(__file__))
    try:
        import importlib.metadata
        version = importlib.metadata.version("geonum")
    except ModuleNotFoundError:
        from pkg_resources import get_distribution
        version = get_distribution('geonum').version
    return (path, version)
