px, py = [], []
mu, me = [], []
for m in matches:
    i = m[1]
    ki = m[2]
    #if d['ph'][i] > 400 and (kyy[ki] > 50 or ky[ki] > .1):
    #if kyy[ki] > 50:
    #if ky[ki] > .2:
    if ky[ki] > .1 and kyy[ki] < 10:
        px.append(200*ky[ki]+kyy[ki])
        py.append(d['ph'][i])
        mu.append(ky[ki])
        me.append(kyy[ki])
