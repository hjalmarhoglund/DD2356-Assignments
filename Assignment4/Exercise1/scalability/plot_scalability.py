import matplotlib.pyplot as plt

t1 = None
p, t = [], []

with open('times.txt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        exe, ranks, secs = line.split()
        ranks, secs = int(ranks), float(secs)

        if exe == 'original.out' and ranks == 1:
            t1 = secs
        elif exe == 'halo.out':
            p.append(ranks)
            t.append(secs)

# speed-up & efficiency
S = [t1 / x for x in t]
E = [s / r for s, r in zip(S, p)]

# plot speed-up
plt.plot(p, S, 'o-')
plt.xscale('log', base=2)
plt.grid(True)
plt.xlabel('MPI processes $p$')
plt.ylabel('Speed-up $S(p)$')
plt.title('Scaling')
plt.tight_layout()
plt.savefig('speedup.png', dpi=150)

