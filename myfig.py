figure()
for shift in [-12, -13, -13.1802, -14]:
    c = (hisparc.analysis.kascade_coincidences
                .search_coincidences(h, k, shift, None))
    hist([abs(x[0] * 1e-9) for x in c], bins=linspace(0, 1, 100),
         histtype='step', label="Shift %s s" % (shift))
legend()
xlabel("Time difference (s)")
ylabel("Counts")
title("Waiting time between triggered and analyzed data")
