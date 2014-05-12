import csv
import tables
import os


class HisparcClusters(tables.IsDescription):
    station_id = tables.UInt32Col(pos=1)
    cluster_id = tables.StringCol(40, pos=2)
    password = tables.StringCol(20, pos=3)
    description = tables.StringCol(20, pos=4)


class HisparcEvent(tables.IsDescription):
    # DISCUSS: use of signed (dflt -1) vs unsigned (labview code)
    event_id = tables.UInt32Col(pos=0)
    timestamp = tables.Time32Col(pos=2)
    nanoseconds = tables.UInt32Col(pos=3)
    ext_timestamp = tables.UInt64Col(pos=4)
    data_reduction = tables.BoolCol(pos=5)
    trigger_pattern = tables.UInt32Col(pos=6)
    baseline = tables.Int16Col(shape=4, dflt=-1, pos=7)
    std_dev = tables.Int16Col(shape=4, dflt=-1, pos=8)
    n_peaks = tables.Int16Col(shape=4, dflt=-1, pos=9)
    pulseheights = tables.Int16Col(shape=4, dflt=-1, pos=10)
    integrals = tables.Int32Col(shape=4, dflt=-1, pos=11)
    traces = tables.Int32Col(shape=4, dflt=-1, pos=12)
    event_rate = tables.Float32Col(pos=13)


class HisparcError(tables.IsDescription):
    event_id = tables.UInt32Col(pos=0)
    timestamp = tables.Time32Col(pos=2)
    messages = tables.Int32Col(pos=3)


class HisparcComparatorData(tables.IsDescription):
    event_id = tables.UInt32Col(pos=0)
    timestamp = tables.Time32Col(pos=2)
    nanoseconds = tables.UInt32Col(pos=3)
    ext_timestamp = tables.UInt64Col(pos=4)
    device = tables.UInt8Col(pos=5)
    comparator = tables.UInt8Col(pos=6)
    count = tables.UInt16Col(pos=7)


class HisparcConfiguration(tables.IsDescription):
    event_id = tables.UInt32Col()
    timestamp = tables.Time32Col()
    gps_latitude = tables.Float64Col()
    gps_longitude = tables.Float64Col()
    gps_altitude = tables.Float64Col()
    mas_version = tables.Int32Col(dflt=-1)
    slv_version = tables.Int32Col(dflt=-1)
    trig_low_signals = tables.UInt32Col()
    trig_high_signals = tables.UInt32Col()
    trig_external = tables.UInt32Col()
    trig_and_or = tables.BoolCol()
    precoinctime = tables.Float64Col()
    coinctime = tables.Float64Col()
    postcoinctime = tables.Float64Col()
    detnum = tables.UInt16Col()
    password = tables.Int32Col(dflt=-1)
    spare_bytes = tables.UInt8Col()
    use_filter = tables.BoolCol()
    use_filter_threshold = tables.BoolCol()
    reduce_data = tables.BoolCol()
    buffer = tables.Int32Col(dflt=-1)
    startmode = tables.BoolCol()
    delay_screen = tables.Float64Col()
    delay_check = tables.Float64Col()
    delay_error = tables.Float64Col()
    mas_ch1_thres_low = tables.Float64Col()
    mas_ch1_thres_high = tables.Float64Col()
    mas_ch2_thres_low = tables.Float64Col()
    mas_ch2_thres_high = tables.Float64Col()
    mas_ch1_inttime = tables.Float64Col()
    mas_ch2_inttime = tables.Float64Col()
    mas_ch1_voltage = tables.Float64Col()
    mas_ch2_voltage = tables.Float64Col()
    mas_ch1_current = tables.Float64Col()
    mas_ch2_current = tables.Float64Col()
    mas_comp_thres_low = tables.Float64Col()
    mas_comp_thres_high = tables.Float64Col()
    mas_max_voltage = tables.Float64Col()
    mas_reset = tables.BoolCol()
    mas_ch1_gain_pos = tables.UInt8Col()
    mas_ch1_gain_neg = tables.UInt8Col()
    mas_ch2_gain_pos = tables.UInt8Col()
    mas_ch2_gain_neg = tables.UInt8Col()
    mas_ch1_offset_pos = tables.UInt8Col()
    mas_ch1_offset_neg = tables.UInt8Col()
    mas_ch2_offset_pos = tables.UInt8Col()
    mas_ch2_offset_neg = tables.UInt8Col()
    mas_common_offset = tables.UInt8Col()
    mas_internal_voltage = tables.UInt8Col()
    mas_ch1_adc_gain = tables.Float64Col()
    mas_ch1_adc_offset = tables.Float64Col()
    mas_ch2_adc_gain = tables.Float64Col()
    mas_ch2_adc_offset = tables.Float64Col()
    mas_ch1_comp_gain = tables.Float64Col()
    mas_ch1_comp_offset = tables.Float64Col()
    mas_ch2_comp_gain = tables.Float64Col()
    mas_ch2_comp_offset = tables.Float64Col()
    slv_ch1_thres_low = tables.Float64Col()
    slv_ch1_thres_high = tables.Float64Col()
    slv_ch2_thres_low = tables.Float64Col()
    slv_ch2_thres_high = tables.Float64Col()
    slv_ch1_inttime = tables.Float64Col()
    slv_ch2_inttime = tables.Float64Col()
    slv_ch1_voltage = tables.Float64Col()
    slv_ch2_voltage = tables.Float64Col()
    slv_ch1_current = tables.Float64Col()
    slv_ch2_current = tables.Float64Col()
    slv_comp_thres_low = tables.Float64Col()
    slv_comp_thres_high = tables.Float64Col()
    slv_max_voltage = tables.Float64Col()
    slv_reset = tables.BoolCol()
    slv_ch1_gain_pos = tables.UInt8Col()
    slv_ch1_gain_neg = tables.UInt8Col()
    slv_ch2_gain_pos = tables.UInt8Col()
    slv_ch2_gain_neg = tables.UInt8Col()
    slv_ch1_offset_pos = tables.UInt8Col()
    slv_ch1_offset_neg = tables.UInt8Col()
    slv_ch2_offset_pos = tables.UInt8Col()
    slv_ch2_offset_neg = tables.UInt8Col()
    slv_common_offset = tables.UInt8Col()
    slv_internal_voltage = tables.UInt8Col()
    slv_ch1_adc_gain = tables.Float64Col()
    slv_ch1_adc_offset = tables.Float64Col()
    slv_ch2_adc_gain = tables.Float64Col()
    slv_ch2_adc_offset = tables.Float64Col()
    slv_ch1_comp_gain = tables.Float64Col()
    slv_ch1_comp_offset = tables.Float64Col()
    slv_ch2_comp_gain = tables.Float64Col()
    slv_ch2_comp_offset = tables.Float64Col()


class HisparcWeather(tables.IsDescription):
    event_id = tables.UInt32Col(pos=0)
    timestamp = tables.Time32Col(pos=1)
    temp_inside = tables.Float32Col(pos=2)
    temp_outside = tables.Float32Col(pos=3)
    humidity_inside = tables.Int16Col(pos=4)
    humidity_outside = tables.Int16Col(pos=5)
    barometer = tables.Float32Col(pos=6)
    wind_dir = tables.Int16Col(pos=7)
    wind_speed = tables.Int16Col(pos=8)
    solar_rad = tables.Int16Col(pos=9)
    uv = tables.Int16Col(pos=10)
    evapotranspiration = tables.Float32Col(pos=11)
    rain_rate = tables.Float32Col(pos=12)
    heat_index = tables.Int16Col(pos=13)
    dew_point = tables.Float32Col(pos=14)
    wind_chill = tables.Float32Col(pos=15)


def open_or_create_file(data_dir, date):
    """Open an existing file or create a new one

    This function opens an existing PyTables file according to the event
    date.  If the file does not yet exist, a new one is created.

    :param data_dir: the directory containing all data files
    :param date: the event date

    """
    dir = os.path.join(data_dir, '%d/%d' % (date.year, date.month))
    file = os.path.join(dir, '%d_%d_%d.h5' % (date.year, date.month,
                                              date.day))

    if not os.path.exists(dir):
        # create dir and parent dirs with mode rwxr-xr-x
        os.makedirs(dir, 0755)

    return tables.open_file(file, 'a')


def get_or_create_station_group(file, cluster, station_id):
    """Get an existing station group or create a new one

    :param file: the PyTables data file
    :param cluster: the name of the cluster
    :param station_id: the station number

    """
    cluster = get_or_create_cluster_group(file, cluster)
    node_name = 'station_%d' % station_id
    try:
        station = file.get_node(cluster, node_name)
    except tables.NoSuchNodeError:
        station = file.create_group(cluster, node_name,
                                   'HiSPARC station %d data' % station_id)
        file.flush()

    return station


def get_or_create_cluster_group(file, cluster):
    """Get an existing cluster group or create a new one

    :param file: the PyTables data file
    :param cluster: the name of the cluster

    """
    try:
        hisparc = file.get_node('/', 'hisparc')
    except tables.NoSuchNodeError:
        hisparc = file.create_group('/', 'hisparc', 'HiSPARC data')
        file.flush()

    node_name = 'cluster_' + cluster.lower()
    try:
        cluster = file.get_node(hisparc, node_name)
    except tables.NoSuchNodeError:
        cluster = file.create_group(hisparc, node_name,
                                   'HiSPARC cluster %s data' % cluster)
        file.flush()

    return cluster


def get_or_create_node(file, cluster, node):
    """Get an existing node or create a new one

    :param file: the PyTables data file
    :param cluster: the parent (cluster) node
    :param node: the node (e.g. events, blobs)

    """
    try:
        node = file.get_node(cluster, node)
    except tables.NoSuchNodeError:
        if node == 'events':
            node = file.create_table(cluster, 'events', HisparcEvent,
                                   'HiSPARC coincidences table')
        elif node == 'errors':
            node = file.create_table(cluster, 'errors', HisparcError,
                                    'HiSPARC error messages')
        elif node == 'comparator':
            node = file.create_table(cluster, 'comparator',
                                    HisparcComparatorData,
                                    'HiSPARC comparator messages')
        elif node == 'blobs':
            node = file.create_vlarray(cluster, 'blobs',
                                      tables.VLStringAtom(),
                                      'HiSPARC binary data')
        elif node == 'config':
            node = file.create_table(cluster, 'config',
                                    HisparcConfiguration,
                                    'HiSPARC configuration messages')
        elif node == 'weather':
            node = file.create_table(cluster, 'weather',
                                    HisparcWeather,
                                    'HiSPARC weather data')
        file.flush()

    return node
