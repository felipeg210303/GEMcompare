import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from cobra.io import read_sbml_model
from pathlib import Path
from IPython.display import display  # Importar la función


# Configuración estética
plt.style.use('seaborn')
sns.set_palette("husl")
plt.rcParams['figure.figsize'] = (14, 8)

# 1. Cargar modelos (ejemplo con rutas - ajusta a tus archivos)
model_paths = {
    "Modelo 1": "B.-subtilis-FBA\ecBSU1\ecBSU1_bigg.xml",
    "Modelo 2": "B.-subtilis-FBA\iBsu1103\iBsu1103_bigg.xml",
    # ... completa con los 9 modelos
    "Modelo 9": "B.-subtilis-FBA\iBsu1147\iBsu1147_bigg.xml"
}

# 2. Función para extraer estadísticas
def get_model_stats(model):
    return {
        "Reacciones": len(model.reactions),
        "Metabolitos": len(model.metabolites),
        "Genes": len(model.genes)
    }

# 3. Procesar todos los modelos
stats_data = []
for name, path in model_paths.items():
    try:
        model = read_sbml_model(path)
        stats = get_model_stats(model)
        stats["Nombre"] = name
        stats_data.append(stats)
    except Exception as e:
        print(f"Error cargando {name}: {str(e)}")

# Convertir a DataFrame
df = pd.DataFrame(stats_data).set_index("Nombre")

# 4. Visualización interactiva
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Gráfico de reacciones
sns.barplot(data=df.reset_index(), x="Nombre", y="Reacciones", ax=axes[0])
axes[0].tick_params(axis='x', rotation=45)
axes[0].set_title("Comparación de Reacciones")

# Gráfico de metabolitos
sns.barplot(data=df.reset_index(), x="Nombre", y="Metabolitos", ax=axes[1])
axes[1].tick_params(axis='x', rotation=45)
axes[1].set_title("Comparación de Metabolitos")

# Gráfico de genes
sns.barplot(data=df.reset_index(), x="Nombre", y="Genes", ax=axes[2])
axes[2].tick_params(axis='x', rotation=45)
axes[2].set_title("Comparación de Genes")

plt.tight_layout()

# 5. Exportar resultados (opcional)
output_dir = Path("resultados_comparacion")
output_dir.mkdir(exist_ok=True)

# Guardar gráficos
plt.savefig(output_dir / "comparacion_modelos.png", dpi=300, bbox_inches='tight')

# Guardar datos en Excel
df.to_excel(output_dir / "estadisticas_modelos.xlsx")

# Mostrar tabla resumen
display(df.style.background_gradient(cmap='Blues'))

print("✔ Análisis completado. Resultados guardados en:", output_dir.resolve())