import numpy as np
import matplotlib.pyplot as plt

labels = ["Glucosa", "Gluconato", "Glicerol", "Malato", "Malato+Glucosa",
          "Piruvato", "Succinato+L-Glutamato", "Fructosa"]
experimental = np.array([0.59, 0.42, 0.40, 0.57, 0.75, 0.17, 0.22, 0.53])
iyo844 = np.array([0.6242, 0.4091, 0.3083, 0.6116, 0.6242, 0.2741, 0.2926, 0.5067])#pan
# Ajuste lineal (solo para calcular R²)
coeffs = np.polyfit(experimental, iyo844, 1)
poly1d_fn = np.poly1d(coeffs)

y_pred = poly1d_fn(experimental)
ss_res = np.sum((iyo844 - y_pred) ** 2)
ss_tot = np.sum((iyo844 - np.mean(iyo844)) ** 2)
r2 = 1 - (ss_res / ss_tot)

fig, ax = plt.subplots(figsize=(8,8))
ax.set_facecolor('#e9f6fb')

# Línea de identidad
xline = np.linspace(0, 0.85, 100)
ax.plot(xline, xline, linestyle='--', color='gray', linewidth=1.5)

# Colores distintos para cada fuente
colors = plt.cm.tab10(np.linspace(0, 1, len(labels)))

sizes = np.full_like(experimental, 180.0)


# Offsets personalizados para evitar solapamientos
offsets = {
    "Glucosa": (0.01, -0.02),
    "Gluconato": (0.01, -0.03),
    "Glicerol": (0.01, 0.02),
    "Malato": (-0.08, 0.02),
    "Malato+Glucosa": (-0.06, 0.03),
    "Piruvato": (-0.06, -0.03),
    "Succinato+L-Glutamato": (0.01, 0.0),
    "Fructosa": (0.01, -0.02)
}

# Dibujar bolitas sólidas con etiquetas desplazadas
for i, lab in enumerate(labels):
    ax.scatter(experimental[i], iyo844[i], s=sizes[i], color=colors[i],
               edgecolors='k', linewidth=0.8, marker='o', zorder=3)
    dx, dy = offsets.get(lab, (0.01, -0.02))
    ax.text(experimental[i] + dx, iyo844[i] + dy, lab,
            fontsize=10, fontfamily='serif', zorder=5)

# Recuadro con la ecuación y R²
eq_text = f"$R^2$ = {r2:.4f}"
bbox_props = dict(boxstyle="round,pad=0.6", fc="white", ec="gray", alpha=0.9)
ax.text(0.05, 0.75, eq_text, transform=ax.transAxes, fontsize=12,
        va='top', ha='left', bbox=bbox_props)

# Estética
ax.set_xlim(-0.02, 0.82)
ax.set_ylim(-0.02, 0.82)
ax.set_xlabel("Velocidad de crecimiento experimental ($h^{-1}$) ", fontsize=12)
ax.set_ylabel("Velocidad de crecimiento simulada ($h^{-1}$)", fontsize=12)
ax.set_title("Pan genomic model", fontsize=16, pad=12)

ax.set_xticks(np.arange(0, 0.9, 0.2))
ax.set_yticks(np.arange(0, 0.9, 0.2))
ax.grid(which='major', linestyle='-', linewidth=1.0, color='white')
ax.grid(which='minor', linestyle=':', linewidth=0.5, color='gray', alpha=0.3)

fig.tight_layout()
fig.savefig('grafica_paridad.png', dpi=250)
plt.show()
