from functions import entropy_sim, p_bar_rayleigh
import numpy as np

n = 3
r0s = np.linspace(0.38, 0.40, 20)
entropy_lim = 100_000_000

entropies = [entropy_sim(3, "line", "rayleigh", r0, 1, entropy_lim) for r0 in r0s]
print(entropies)

r0max = r0s[np.argmax(entropies)]
p_bar_max = p_bar_rayleigh("line", r0max, 1, 1_000_000)
print(f"r0max: {r0max}, p_bar_max: {p_bar_max}")