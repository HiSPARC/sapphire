from lazy import lazy
from numpy import degrees, radians, log10
import tables

from .particles import name, particle_id


class CorsikaQuery(object):

    def __init__(self, data, simulations_group='/simulations'):
        """Setup variables to point to the tables

        :param data: either a PyTables file or path to a HDF5 file
        :param simulations_group: path to the simulations group

        """
        if not isinstance(data, tables.File):
            self.data = tables.open_file(data, 'r')
        else:
            self.data = data
        self.sims = self.data.get_node(simulations_group)

    def finish(self):
        """Clean-up after using

        Do not use if you opened the file yourself and still intent to
        use it.

        """
        self.data.close()

    def all_simulations(self, iterator=False):
        """Get all simulations

        :return: all simulations.

        """
        if iterator:
            simulations = self.sims.iterrows()
        else:
            simulations = self.sims.read()
        return simulations

    def seeds(self, simulations, iterator=False):
        """Get combined seeds for set of simulations

        :return: combined seed1 and seed2.

        """
        seeds = ('%d_%d' % (sim['seed1'], sim['seed2']) for sim in simulations)
        if not iterator:
            seeds = list(seeds)
        return seeds

    def get_info(self, seeds):
        """Get info for a simulation

        :param seeds: combined string with seed1 and seed2.
        :return: row matching the seeds.

        """
        seed1, seed2 = seeds.split('_')
        queries = []
        queries.append(self.filter('seed1', seed1))
        queries.append(self.filter('seed2', seed2))
        query = ' & '.join(queries)
        simulation = self.perform_query(query, iterator=False)[0]

        return simulation

    @lazy
    def all_energies(self):
        """All available energies

        :return: set of available simulation energies in log10(eV).

        """
        return set(log10(self.sims.col('energy')))

    @lazy
    def all_particles(self):
        """All available particles

        :return: set of available simulation particles.

        """
        return {name(p_id) for p_id in set(self.sims.col('particle_id'))}

    @lazy
    def all_azimuths(self):
        """All available azimuths

        :return: set of available simulation azimuths.

        """
        return {degrees(azimuth) for azimuth in set(self.sims.col('azimuth'))}

    @lazy
    def all_zeniths(self):
        """All available zeniths

        :return: set of available simulation zeniths.

        """
        return {degrees(zenith) for zenith in set(self.sims.col('zenith'))}

    def simulations(self, particle='proton', energy=None, zenith=None,
                    azimuth=None, iterator=False):
        """Set of available energies given the requirements

        :param particle: primary particle must be this kind, name of particle.
                         Defaults to proton.
        :param energy: primary energy must be this value, in log10(eV).
        :param zenith: shower zenith must be this value, in degrees.
        :param azimuth: shower azimuth must be this value, in degrees.
        :return: simulations matching the query.

        """
        queries = []
        if particle is not None:
            if particle not in self.all_particles:
                raise RuntimeError('Particle not available')
            queries.append(self.filter('particle_id', particle_id(particle)))
        if energy is not None:
            if energy not in self.all_energies:
                raise RuntimeError('Energy not available')
            queries.append(self.float_filter('log10(energy)', energy))
        if zenith is not None:
            queries.append(self.float_filter('zenith', radians(zenith)))
        if azimuth is not None:
            queries.append(self.float_filter('azimuth', radians(azimuth)))
        query = ' & '.join(queries)

        filtered_simulations = self.perform_query(query, iterator)

        return filtered_simulations

    def available_parameters(self, parameter, *args, **kwargs):
        """Get set of available values of type parameter for a subset

        :param parameter: name of the parameter for which the available
                          values are returned.
        :param: use the arguments available to :meth:`simulations` to filter
                the simulations.
        :return: set of available values.

        """
        sims = self.simulations(*args, **kwargs)
        available = set(sims[parameter])
        if parameter == 'energy':
            return {log10(energy) for energy in available}
        elif parameter in ['zenith', 'azimuth']:
            return {degrees(angle) for angle in available}
        elif parameter == 'particle_id':
            return {name(particle) for particle in available}
        else:
            return available

    def filter(self, type, value):
        """Filter to be in a range

        :param type: variable to filter.
        :param value: value to match.
        :return: query.

        """
        query = '(%s == %s)' % (type, value)

        return query

    def float_filter(self, type, value):
        """Filter float values

        Take into account that the values are likely not a perfect match.

        :param type: variable to filter.
        :param value: value to match.
        :return: query.

        """
        query = '(abs(%s - %s) < 1e-4)' % (type, value)

        return query

    def range_filter(self, type, min=None, max=None):
        """Filter to be in a range

        :param type: variable to filter.
        :param min,max: limits on the value.
        :return: query.

        """
        queries = []
        if min is not None:
            queries.append('(%s >= %s)' % (type, min))
        if max is not None:
            queries.append('(%s <= %s)' % (type, max))
        query = ' & '.join(queries)

        return query

    def perform_query(self, query, iterator=False):
        """Perform a query on the simulations table

        :param query: a valid PyTables query string for the simulations table.
        :return: simulations matching the query.

        """
        if query:
            if iterator:
                filtered_simulations = self.sims.where(query)
            else:
                filtered_simulations = self.sims.read_where(query)
        else:
            filtered_simulations = self.all_simulations(iterator)

        return filtered_simulations
