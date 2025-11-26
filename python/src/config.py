"""
Configuración central del proyecto Ising híbrido.
Unidades naturales: k_B = 1, ℏ = 1
"""

from pathlib import Path
from enum import Enum

# ============================================================
# DIRECTORIOS
# ============================================================
from pathlib import Path

# Obtener la carpeta donde está este archivo (src/)
SRC_DIR = Path(__file__).parent.absolute()

# Subir dos niveles: src/ → python/ → PROYECTO_ISING_HIBRIDO/
PYTHON_DIR = SRC_DIR.parent.absolute()
PROJECT_ROOT = PYTHON_DIR.parent.absolute()

# Los datos están en PROYECTO_ISING_HIBRIDO/cpp/data/
DATA_DIR = PROJECT_ROOT / 'cpp' / 'data'

# Directorios de salida en PROYECTO_ISING_HIBRIDO/graficas/
OUTPUT_DIR = PROJECT_ROOT / 'graficas'

# Subdirectorios de entrada (en DATA_DIR)
PARA_DATA_DIR = DATA_DIR / 'paramagnetismo'
HYST_DATA_DIR = DATA_DIR / 'ferromagnetismo'
TRANS_DATA_DIR = DATA_DIR / 'transicion'
SNAP_DATA_DIR = DATA_DIR / 'snapshots'
RELAX_DATA_DIR = DATA_DIR / 'relax'
TC_DATA_DIR = DATA_DIR / 'temperatura_critica'

# Subdirectorios de salida
PARA_OUT_DIR = OUTPUT_DIR / 'paramagnetismo'
HYST_OUT_DIR = OUTPUT_DIR / 'ferromagnetismo'
TRANS_OUT_DIR = OUTPUT_DIR / 'transicion'
SNAP_OUT_DIR = OUTPUT_DIR / 'snapshots'
COMP_OUT_DIR = OUTPUT_DIR / 'comparacion_z'
RELAX_OUT_DIR = OUTPUT_DIR / 'relajacion'
ANALYSIS_OUT_DIR = OUTPUT_DIR / 'analisis'
TC_OUT_DIR = OUTPUT_DIR / 'temperatura_critica'
M_VS_T_OUT_DIR = OUTPUT_DIR / 'm_vs_T'

# ============================================================
# PARÁMETROS FÍSICOS Y DE SIMULACIÓN
# ============================================================

class Lattice(Enum):
    """Tipos de redes cristalinas"""
    CHAIN = ('chain', 2)
    HONEYCOMB = ('honeycomb', 3)
    SQUARE = ('square', 4)
    BCC = ('bcc', 8)

    @property
    def name_str(self) -> str:
        return self.value[0]

    @property
    def coordination_number(self) -> int:
        """Número de coordinación z"""
        return self.value[1]


LATTICES = [Lattice.CHAIN, Lattice.HONEYCOMB, Lattice.SQUARE, Lattice.BCC]
Q_VALUES = [0.5, 0.8, 1.0]
T_VALUES_SNAPSHOT = [1.0, 3.0, 5.0]

# Unidades naturales: k_B = 1
BOLTZMANN_CONSTANT = 1.0  # Sistema de unidades naturales

# ============================================================
# PARÁMETROS DE GRAFICADO Y ESTÉTICA
# ============================================================

FIGURE_DPI = 150
FIGURE_DPI_HIGH = 300
DEFAULT_FIGSIZE = (10, 6)
WIDE_FIGSIZE = (15, 4)
GRID_FIGSIZE = (12, 10)

# Colores y estilos
CMAP_SPIN = 'RdBu'  # Colormap para visualización de spins
MARKER_STYLE = 'o-'
MARKER_SIZE = 6
LINE_WIDTH = 1.5
LINE_WIDTH_THICK = 2.5

# ============================================================
# ETIQUETAS Y UNIDADES (Sistema natural: k_B = 1)
# ============================================================

LABELS = {
    'temperature': r'$T$ (en unidades de $J/k_B$)',
    'field': r'$H$ (en unidades de $J$)',
    'magnetization': r'Magnetización $m$ (adimensional)',
    'energy': r'Energía $E$ (en unidades de $J$)',
    'susceptibility': r'Susceptibilidad $\chi$ (adimensional)',
    'concentration': r'Concentración $q$',
    'coordination': r'Número de coordinación $z$',
    'mc_steps': 'Pasos de Monte Carlo',
    'reduced_field': r'$H/T$ (adimensional)',
}

TITLES = {
    'paramagnetism': 'Paramagnetismo',
    'reduced_magnetization': 'Ley de Estados Correspondientes',
    'susceptibility': 'Susceptibilidad Magnética',
    'hysteresis': 'Histéresis Magnética',
    'phase_transition': 'Transición de Fase',
    'energy_transition': 'Energía vs Temperatura',
    'relaxation': 'Relajación Dinámica',
    'tanh_fit': 'Ajuste a Ley de Tanh',
    'critical_temp_vs_q': r'Temperatura Crítica $T_c$ vs Concentración $q$',
    'critical_temp_vs_z': r'Temperatura Crítica $T_c$ vs Número de Coordinación $z$',
}

# ============================================================
# CONSTANTES DE ANÁLISIS
# ============================================================

TANH_FIT_MAXFEV = 5000
TANH_FIT_P0 = [1.0, 1.0]
SMOOTHING_WINDOW_FRAC = 1/50  # Fracción de datos para suavizado

# Tolerancia para encontrar Tc (50% de magnetización máxima)
TC_TOLERANCE_FRACTION = 0.5

# Temperaturas para histéresis (ferromagnetismo)
T_VALUES_HYSTERESIS = [2.0, 5.0]

# Temperaturas para relajación
T_VALUES_RELAXATION = [1.0, 3.0, 5.0]

# Valores de J para relajación
J_VALUES_RELAXATION = [0, 1]