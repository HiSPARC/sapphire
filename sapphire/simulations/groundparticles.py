from base import BaseSimulation


class GroundParticlesSimulation(BaseSimulation):

    def __init__(self, corsikafile, *args, **kwargs):
        super(GroundParticlesSimulation, self).__init__(*args, **kwargs)
        
        self.corsikafile = corsikafile

    def generate_shower_parameters(self):
        """Generate shower parameters like core position, energy, etc."""

        pass

    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower."""

        # implement this!
        observables = None

        return observables

    def simulate_trigger(self, station_observables):
        """Simulate a trigger response."""

        return True
