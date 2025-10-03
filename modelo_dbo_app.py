pip install matplotlib
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

st.set_page_config(page_title="Modelo de DBO en cuerpos conectados", layout="centered")

# ----------------------------------------
# TÍTULO Y DESCRIPCIÓN
# ----------------------------------------
st.title("Simulación de DBO en cuerpos de agua interconectados")
st.markdown("""
Este modelo simula la evolución de la **Demanda Bioquímica de Oxígeno (DBO)** en dos cuerpos de agua conectados, considerando:

- Transporte dispersivo
- Degradación biológica (1er orden)
- **Descarga puntual solo durante el primer día**
- Condición de frontera con la bahía

🧪 Modifica los parámetros y observa los resultados.
""")

# ----------------------------------------
# ENTRADAS INTERACTIVAS
# ----------------------------------------

st.sidebar.header("🔧 Parámetros modificables")

L_input = st.sidebar.slider("Carga de DBO durante el primer día [mg/día]", 0, 5_000_000, 1_000_000, step=100_000)

CB = st.sidebar.slider("Concentración de DBO en la bahía [mg/L]", 0.0, 10.0, 1.0, step=0.1)

t_descarga = 1.0  # [días] duración de la descarga

# ----------------------------------------
# PARÁMETROS FIJOS DEL MODELO
# ----------------------------------------

# Volúmenes [m³]
V1 = 1e5
V2 = 2e5

# Conexión bahía ↔ Cuerpo 1
A1B = 25      # m²
L1B = 80      # m
E1B = 0.6     # m²/s
D1B = E1B * A1B / L1B * 86400  # m³/día

# Conexión Cuerpo 1 ↔ Cuerpo 2
A12 = 18
L12 = 120
E12 = 0.4
D12 = E12 * A12 / L12 * 86400

# Remoción de DBO
k1 = 0.2  # 1/día
k2 = 0.1  # 1/día

# Tiempo de simulación
Tmax = 30  # días
dt = 0.1
N = int(Tmax / dt)

# Condiciones iniciales
C1 = 0.0
C2 = 0.0

# ----------------------------------------
# FUNCIONES DEL MODELO
# ----------------------------------------

def derivadas(C1, C2, t):
    # Descarga puntual solo durante el primer día
    L = L_input if t < t_descarga else 0.0

    dC1_dt = (1/V1)*(D1B*(CB - C1) + D12*(C2 - C1) + L) - k1*C1
    dC2_dt = (1/V2)*(D12*(C1 - C2)) - k2*C2
    return dC1_dt, dC2_dt

# ----------------------------------------
# RUNGE-KUTTA 4
# ----------------------------------------

t_list = [0]
C1_list = [C1]
C2_list = [C2]

t = 0

for i in range(N):
    k1_1, k2_1 = derivadas(C1, C2, t)
    k1_2, k2_2 = derivadas(C1 + dt*k1_1/2, C2 + dt*k2_1/2, t + dt/2)
    k1_3, k2_3 = derivadas(C1 + dt*k1_2/2, C2 + dt*k2_2/2, t + dt/2)
    k1_4, k2_4 = derivadas(C1 + dt*k1_3, C2 + dt*k2_3, t + dt)

    C1 += (dt/6)*(k1_1 + 2*k1_2 + 2*k1_3 + k1_4)
    C2 += (dt/6)*(k2_1 + 2*k2_2 + 2*k2_3 + k2_4)

    t += dt
    t_list.append(t)
    C1_list.append(C1)
    C2_list.append(C2)

# ----------------------------------------
# GRÁFICA DE RESULTADOS
# ----------------------------------------

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(t_list, C1_list, label="Cuerpo 1 (con descarga)", linewidth=2)
ax.plot(t_list, C2_list, label="Cuerpo 2 (ciénaga)", linewidth=2)
ax.axhline(CB, color='gray', linestyle='--', label='Bahía (frontera)', linewidth=1.5)
ax.set_xlabel("Tiempo [días]")
ax.set_ylabel("Concentración de DBO [mg/L]")
ax.set_title("Evolución de la DBO con descarga temporal")
ax.grid(True)
ax.legend()
st.pyplot(fig)
