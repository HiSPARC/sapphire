import tables


DATAFILE = 'data.h5'


class CoreReconstruction(object):
    def __init__(self, data, simulation_path):
        self.data = data
        self.simulation = data.getNode(simulation_path)
        self.cluster = self.simulation._v_attrs.cluster

    def draw_cluster(self):
        for station in self.cluster.stations:
            for detector in station.detectors:
                x, y = detector.get_xy_coordinates()
                scatter(x, y, c='r', s=5, edgecolor='none')
            x, y, alpha = station.get_xyalpha_coordinates()
            scatter(x, y, c='orange', s=10, edgecolor='none')

    def plot_coincidence(self, index=0, multiplicity=3):
        clf()
        coincidence = self.get_coincidence_with_multiplicity(index,
                                                             multiplicity)
        scatter(coincidence['x'], coincidence['y'], c='b', s=10)

        for event in self.get_events_from_coincidence(coincidence):
            self.plot_event_on_map(event)

    def get_coincidence_with_multiplicity(self, index, multiplicity):
        coincidences = self.simulation.coincidences.read()
        sel = coincidences.compress(coincidences[:]['N'] >= multiplicity)
        return sel[index]

    def get_events_from_coincidence(self, coincidence):
        events = []
        id = coincidence['id']

        for index in self.simulation.c_index[id]:
            events.append(self.simulation.observables[index])

        return events

    def plot_event_on_map(self, event):
        station_id = event['station_id']
        station = self.cluster.stations[station_id - 1]

        detectors = station.detectors
        num_particles = [event[u] for u in ['n1', 'n2', 'n3', 'n4']]

        for detector, num in zip(detectors, num_particles):
            self.plot_detector_on_map(detector, num)

    def plot_detector_on_map(self, detector, num_particles):
        x, y = detector.get_xy_coordinates()
        size = num_particles * 10
        scatter(x, y, c='r', s=size)


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'r')

    c = CoreReconstruction(data, '/ldfsim')
