"""
Utilidades y funciones auxiliares para el proyecto Ising.
"""

from pathlib import Path
from typing import Optional, List
import numpy as np
from scipy.optimize import curve_fit


# ============================================================
# UTILIDADES DE ARCHIVOS
# ============================================================

def read_exact(f, n: int) -> bytes:
    """Lee exactamente n bytes o lanza EOFError si no hay suficientes."""
    buf = f.read(n)
    if len(buf) != n:
        raise EOFError(f"Se esperaban {n} bytes pero se recibieron {len(buf)}")
    return buf


def try_multiple_filenames(base_candidates: List[Path]) -> Optional[Path]:
    """Devuelve la primera ruta existente de la lista o None."""
    for p in base_candidates:
        if p.exists():
            return p
    return None


def format_q_str(q: float) -> str:
    """Convierte q=0.8 a formato 'q80' consistente con C++"""
    return f"q{int(round(q*100)):02d}"


def ensure_dir(path: Path) -> Path:
    """Crea directorio si no existe y devuelve la ruta."""
    path.mkdir(parents=True, exist_ok=True)
    return path


# ============================================================
# FUNCIONES DE AJUSTE
# ============================================================

def tanh_model(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """
    Modelo de magnetización: m(H/T) = a * tanh(b * H/T)
    
    Representa el comportamiento paramagnético esperado.
    """
    return a * np.tanh(b * x)


def fit_tanh(x_data: np.ndarray, y_data: np.ndarray, 
             p0: tuple = (1.0, 1.0),
             maxfev: int = 5000) -> dict:
    """
    Ajusta datos a modelo tanh con análisis estadístico.
    
    Parameters
    ----------
    x_data : np.ndarray
        Variable independiente (H/T)
    y_data : np.ndarray
        Magnetización observada
    p0 : tuple
        Parámetros iniciales [a, b]
    maxfev : int
        Máximo número de evaluaciones
        
    Returns
    -------
    Dict con: popt, pcov, r_squared, y_pred, residuals, success
    """
    try:
        popt, pcov = curve_fit(tanh_model, x_data, y_data, 
                               p0=p0, maxfev=maxfev)
        a_opt, b_opt = popt
        
        # Predicciones
        y_pred = tanh_model(x_data, a_opt, b_opt)
        
        # R²
        ss_res = np.sum((y_data - y_pred) ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
        
        # Residuales
        residuals = y_data - y_pred
        
        return {
            'popt': popt,
            'pcov': pcov,
            'a': a_opt,
            'b': b_opt,
            'r_squared': r_squared,
            'y_pred': y_pred,
            'residuals': residuals,
            'success': True
        }
    except Exception as e:
        return {
            'success': False,
            'error': str(e)
        }


def estimate_critical_temperature(T: np.ndarray, m: np.ndarray,
                                   tolerance: float = 0.5) -> float:
    """
    Estima temperatura crítica como punto donde m = tolerance * m_max
    
    Parameters
    ----------
    T : np.ndarray
        Array de temperaturas
    m : np.ndarray
        Array de magnetizaciones
    tolerance : float
        Fracción de m_max para definir Tc (default: 0.5)
        
    Returns
    -------
    float
        Temperatura crítica estimada
    """
    if len(T) < 2 or len(m) < 2:
        return 0.0
    
    m_max = np.max(np.abs(m))
    target_m = m_max * tolerance
    idx = np.argmin(np.abs(np.abs(m) - target_m))
    
    return float(T[idx])


def calculate_critical_exponent(T: np.ndarray, m: np.ndarray, 
                                 Tc: float, 
                                 temp_range: float = 0.2) -> dict:
    """
    Estima exponente crítico β usando ley de potencias.
    
    m ~ |T - Tc|^β cerca de Tc
    
    Parameters
    ----------
    T : np.ndarray
        Temperaturas
    m : np.ndarray
        Magnetizaciones
    Tc : float
        Temperatura crítica
    temp_range : float
        Rango (fracción de Tc) alrededor de Tc para ajuste
        
    Returns
    -------
    Dict con parámetros del ajuste de ley de potencias
    """
    # Seleccionar región cerca de Tc
    mask = np.abs(T - Tc) < Tc * temp_range
    if np.sum(mask) < 3:
        return {'success': False, 'message': 'Insuficientes puntos cerca de Tc'}
    
    T_near = T[mask]
    m_near = m[mask]
    
    try:
        def power_law(t, coeff, beta):
            return coeff * np.abs(t - Tc) ** beta
        
        popt, _ = curve_fit(power_law, T_near, np.abs(m_near),
                           p0=[1.0, 0.5], maxfev=1000)
        
        return {
            'success': True,
            'coefficient': popt[0],
            'beta': popt[1],
            'Tc': Tc
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}


# ============================================================
# UTILIDADES DE VISUALIZACIÓN
# ============================================================

def smooth_array(arr: np.ndarray, window_frac: float = 1/50) -> np.ndarray:
    """Suaviza un array usando media móvil"""
    window = max(1, int(len(arr) * window_frac))
    if window < 2:
        return arr
    return np.convolve(arr, np.ones(window)/window, mode='valid')

# ============================================================
# ANÁLISIS ESTADÍSTICO
# ============================================================

def calculate_blocking_error(time_series: np.ndarray) -> float:
    """
    Calcula el error estándar usando el método de bloques (Blocking Method).
    Corrige la subestimación del error debida a la autocorrelación en Monte Carlo.
    """
    n = len(time_series)
    if n < 50:
        # Fallback para series muy cortas: error estándar simple
        return float(np.std(time_series, ddof=1) / np.sqrt(n))

    # Comenzamos con el error naive (bloque tamaño 1)
    max_error = np.std(time_series, ddof=1) / np.sqrt(n)

    # Probamos bloques de tamaño 2, 4, 8...
    block_size = 2
    while True:
        num_blocks = n // block_size
        if num_blocks < 10:  # Necesitamos suficientes bloques para tener estadística fiable
            break

        # Recortar los datos para que encajen perfectamente en bloques enteros
        trimmed_len = num_blocks * block_size
        blocks = time_series[:trimmed_len].reshape(num_blocks, block_size)

        # 1. Calcular el promedio de cada bloque
        block_means = np.mean(blocks, axis=1)

        # 2. Calcular el error estándar de esos promedios
        # Error = Desviación(medias) / Raíz(Num_bloques)
        current_error = np.std(block_means, ddof=1) / np.sqrt(num_blocks)

        # Si el error crece (se revela la correlación), actualizamos el máximo
        if current_error > max_error:
            max_error = current_error

        block_size *= 2

    return float(max_error)


def calculate_hysteresis_area(H_down: np.ndarray, m_down: np.ndarray, 
                               H_up: np.ndarray, m_up: np.ndarray) -> float:
    """
    Calcula el área encerrada en un ciclo de histéresis.
    """
    # Compatibilidad NumPy 1.x y 2.x
    trapz_func = getattr(np, 'trapezoid', np.trapz)
    
    # Ordenar ambas ramas por H (ascendente)
    idx_down = np.argsort(H_down)
    idx_up = np.argsort(H_up)
    
    H_down_sorted = H_down[idx_down]
    m_down_sorted = m_down[idx_down]
    H_up_sorted = H_up[idx_up]
    m_up_sorted = m_up[idx_up]
    
    area_down = float(trapz_func(m_down_sorted, H_down_sorted))
    area_up = float(trapz_func(m_up_sorted, H_up_sorted))
    
    return abs(area_down - area_up)

def build_filename_candidates(base_dir: Path, pattern: str, **kwargs) -> List[Path]:
    """
    Genera lista de candidatos de nombres de archivo basados en un patrón.
    
    Parameters
    ----------
    base_dir : Path
        Directorio base donde buscar
    pattern : str
        Patrón del nombre con placeholders {param}
        Ejemplo: 'hysteresis_{lattice}_q{q}_T{T}.bin'
    **kwargs : dict
        Valores para los placeholders
        
    Returns
    -------
    List[Path]
        Lista de rutas candidatas (más específica primero)
        
    Examples
    --------
    >>> candidates = build_filename_candidates(
    ...     HYST_DATA_DIR, 
    ...     'hysteresis_{lattice}_q{q}_T{T}.bin',
    ...     lattice='square', q=80, T=20
    ... )
    """
    candidates = []
    
    # Intentar con todos los parámetros
    try:
        full_name = pattern.format(**kwargs)
        candidates.append(base_dir / full_name)
    except KeyError:
        pass
    
    # Intentar sin temperatura (retrocompatibilidad)
    if 'T' in kwargs:
        kwargs_no_T = {k: v for k, v in kwargs.items() if k != 'T'}
        try:
            no_T_pattern = pattern.replace('_T{T}', '')
            name_no_T = no_T_pattern.format(**kwargs_no_T)
            candidates.append(base_dir / name_no_T)
        except KeyError:
            pass
    
    # Intentar sin q (retrocompatibilidad)
    if 'q' in kwargs:
        kwargs_no_q = {k: v for k, v in kwargs.items() if k != 'q'}
        try:
            no_q_pattern = pattern.replace('_q{q}', '')
            name_no_q = no_q_pattern.format(**kwargs_no_q)
            candidates.append(base_dir / name_no_q)
        except KeyError:
            pass
    
    return candidates