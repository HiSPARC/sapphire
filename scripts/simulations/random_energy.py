flux = lambda x: x ** -2.7


def random_energy(a, b, size=1):
    y0 = flux(a)
    y1 = flux(b)

    energies = []
    for i in range(size):
        while True:
            x = random.uniform(a, b)
            y = random.uniform(y0, y1)

            if y <= flux(x):
                break

        energies.append(x)

    if size == 1:
        return energies[0]
    else:
        return array(energies)
