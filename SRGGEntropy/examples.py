from functions import square_entropy_hard, square_entropy_rayleigh, square_entropy_triangular, torus_1d_entropy_rayleigh, p_bar_rayleigh, line_entropy_rayleigh, line_entropy_hard, entropy_sim
import numpy as np
import matplotlib.pyplot as plt

# plot the entropy curve

n = 3
r0s = np.linspace(0.38, 0.40, 20)
lim = 1_000_000

entropies = [entropy_sim(3, "line", "rayleigh", r0, 1, 50_000_000) for r0 in r0s]
print(entropies)
r0max = r0s[np.argmax(entropies)]
p_bar_max = p_bar_rayleigh("line", r0max, 1, lim)
print(f"r0max: {r0max}, p_bar_max: {p_bar_max}")

plt.plot(r0s, entropies)
plt.xlabel('r0')
plt.ylabel('Entropy')
plt.title('Entropy curves for n=3')
plt.show()