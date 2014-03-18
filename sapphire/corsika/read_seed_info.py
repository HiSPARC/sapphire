import tables


DAT_URL = '/data/hisparc/corsika/seed_info.h5'


def get_seeds():
    """Get information from seed_info.h5 and print the seeds."""

    seed_info = tables.open_file(DAT_URL, 'r')
    seeds = seed_info.root.simulations
    results = seeds.read_where('energy > 1E13')
    for result in results:
        print '%d_%d/corsika.h5' % (result['seed1'], result['seed2'])
    seed_info.close()

    return results

if __name__ == '__main__':
    results = get_seeds()
