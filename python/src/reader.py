"""
Lectores robustos de archivos binarios generados por simulaciones C++.
"""

from pathlib import Path
from typing import Optional, Dict, Any, List
import struct
import numpy as np


# Funciones auxiliares locales
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


class IsingDataReader:
    """Lectores especializados para diferentes tipos de datos Ising"""

    @staticmethod
    def read_paramagnetism(filename: str) -> Optional[Dict[str, Any]]:
        """Lee archivo de paramagnetismo: L, q, J, T[], H[], m[][], chi[][]"""
        filepath = Path(filename)
        candidates = [filepath, filepath.with_suffix('.bin')]
        found = try_multiple_filenames(candidates)
        
        if not found:
            print(f"⚠ Archivo no encontrado: {filename}")
            return None

        try:
            with open(found, 'rb') as f:
                L = struct.unpack('<I', read_exact(f, 4))[0]
                q = struct.unpack('<f', read_exact(f, 4))[0]
                J = struct.unpack('<f', read_exact(f, 4))[0]
                n_T = struct.unpack('<I', read_exact(f, 4))[0]
                n_H = struct.unpack('<I', read_exact(f, 4))[0]

                lattice = read_exact(f, 32).decode('utf-8', errors='ignore').rstrip('\x00')
                temperatures = np.frombuffer(read_exact(f, n_T * 4), dtype=np.float32).copy()
                fields = np.frombuffer(read_exact(f, n_H * 4), dtype=np.float32).copy()

                m_data = np.zeros((n_T, n_H), dtype=np.float32)
                for i in range(n_T):
                    m_data[i, :] = np.frombuffer(read_exact(f, n_H * 4), dtype=np.float32).copy()

                chi_data = np.zeros((n_T, n_H), dtype=np.float32)
                for i in range(n_T):
                    chi_data[i, :] = np.frombuffer(read_exact(f, n_H * 4), dtype=np.float32).copy()

            return {
                'L': L, 'q': q, 'J': J, 'lattice': lattice,
                'T': temperatures, 'H': fields, 'm': m_data, 'chi': chi_data
            }
        except EOFError as e:
            print(f"❌ EOF leyendo {found}: {e}")
            return None
        except Exception as e:
            print(f"❌ Error leyendo {found}: {e}")
            return None

    @staticmethod
    def read_hysteresis(filename: str) -> Optional[Dict[str, Any]]:
        """Lee archivo de histéresis: n, H_down[], m_down[], H_up[], m_up[]"""
        filepath = Path(filename)
        candidates = [filepath, filepath.with_suffix('.bin')]
        found = try_multiple_filenames(candidates)
        
        if not found:
            print(f"⚠ Archivo no encontrado: {filename}")
            return None

        try:
            with open(found, 'rb') as f:
                n = struct.unpack('<I', read_exact(f, 4))[0]
                H_down = np.frombuffer(read_exact(f, n * 4), dtype=np.float32).copy()
                m_down = np.frombuffer(read_exact(f, n * 4), dtype=np.float32).copy()
                H_up = np.frombuffer(read_exact(f, n * 4), dtype=np.float32).copy()
                m_up = np.frombuffer(read_exact(f, n * 4), dtype=np.float32).copy()
            
            return {'H_down': H_down, 'm_down': m_down, 'H_up': H_up, 'm_up': m_up}
        except Exception as e:
            print(f"❌ Error leyendo {found}: {e}")
            return None

    @staticmethod
    def read_phase_transition(filename: str) -> Optional[Dict[str, Any]]:
        """Lee archivo de transición: n, T[], m[], E[]"""
        filepath = Path(filename)
        candidates = [filepath, filepath.with_suffix('.bin')]
        found = try_multiple_filenames(candidates)
        
        if not found:
            print(f"⚠ Archivo no encontrado: {filename}")
            return None

        try:
            with open(found, 'rb') as f:
                n = struct.unpack('<I', read_exact(f, 4))[0]
                T = np.frombuffer(read_exact(f, n * 4), dtype=np.float32).copy()
                m = np.frombuffer(read_exact(f, n * 4), dtype=np.float32).copy()
                E = np.frombuffer(read_exact(f, n * 4), dtype=np.float32).copy()
            
            return {'T': T, 'm': m, 'E': E}
        except Exception as e:
            print(f"❌ Error leyendo {found}: {e}")
            return None

    @staticmethod
    def read_snapshot(filename: str) -> Optional[Dict[str, Any]]:
        """Lee snapshot de configuración: dims[], T, spins[], coordinates[]"""
        filepath = Path(filename)
        candidates = [filepath, filepath.with_suffix('.bin')]
        found = try_multiple_filenames(candidates)
        
        if not found:
            print(f"⚠ Snapshot no encontrado: {filename}")
            return None

        try:
            with open(found, 'rb') as f:
                # Leer dimensiones y shape
                n_dims = struct.unpack('<I', read_exact(f, 4))[0]
                dims: List[int] = []
                for _ in range(n_dims):
                    dims.append(struct.unpack('<I', read_exact(f, 4))[0])
                
                shape = tuple(dims)
                
                # Leer temperatura y número de sitios ocupados
                temp = struct.unpack('<f', read_exact(f, 4))[0]
                n_occupied = struct.unpack('<I', read_exact(f, 4))[0]
                
                # Leer valores de spins
                spins = np.frombuffer(read_exact(f, n_occupied), dtype=np.int8).copy()
                
                # =====================================================
                # NUEVO: Leer coordenadas de sitios ocupados
                # =====================================================
                coordinates: List[tuple] = []
                
                if n_dims <= 2:
                    # Para 1D y 2D: leer pares (i, j)
                    for _ in range(n_occupied):
                        i = struct.unpack('<I', read_exact(f, 4))[0]
                        j = struct.unpack('<I', read_exact(f, 4))[0]
                        coordinates.append((i, j))
                else:
                    # Para 3D: leer triplas (i, j, k)
                    for _ in range(n_occupied):
                        i = struct.unpack('<I', read_exact(f, 4))[0]
                        j = struct.unpack('<I', read_exact(f, 4))[0]
                        k = struct.unpack('<I', read_exact(f, 4))[0]
                        coordinates.append((i, j, k))

            return {
                'shape': shape,
                'T': temp,
                'n_occupied': n_occupied,
                'spins': spins,
                'coordinates': coordinates
            }
        except Exception as e:
            print(f"❌ Error leyendo snapshot {found}: {e}")
            return None
        

    @staticmethod
    def read_relaxation(filename: str) -> Optional[Dict[str, Any]]:
        """Lee archivo de relajación: n, E[]"""
        filepath = Path(filename)
        candidates = [filepath, filepath.with_suffix('.bin'), filepath.parent / (filepath.stem + '.bin')]
        found = try_multiple_filenames(candidates)
        
        if not found:
            print(f"⚠ Archivo de relajación no encontrado: {filename}")
            return None

        try:
            with open(found, 'rb') as f:
                n = struct.unpack('<I', read_exact(f, 4))[0]
                energies = np.frombuffer(read_exact(f, n * 4), dtype=np.float32).copy()
            
            return {'E': energies}
        except Exception as e:
            print(f"❌ Error leyendo relajación {found}: {e}")
            return None

    @staticmethod
    def read_tc_analysis(filename: str) -> Optional[Dict[str, Any]]:
        """Lee archivo de análisis Tc: lattices[], q_values[], tc_matrix[][]"""
        filepath = Path(filename)
        candidates = [filepath, filepath.with_suffix('.bin')]
        found = try_multiple_filenames(candidates)
        
        if not found:
            print(f"⚠ Archivo de análisis Tc no encontrado: {filename}")
            return None

        try:
            with open(found, 'rb') as f:
                n_lattices = struct.unpack('<I', read_exact(f, 4))[0]
                n_q = struct.unpack('<I', read_exact(f, 4))[0]

                lattices = []
                for _ in range(n_lattices):
                    lattice = read_exact(f, 32).decode('utf-8', errors='ignore').rstrip('\x00')
                    lattices.append(lattice)

                q_values = np.frombuffer(read_exact(f, n_q * 4), dtype=np.float32).copy()

                tc_matrix = []
                for _ in range(n_lattices):
                    row = np.frombuffer(read_exact(f, n_q * 4), dtype=np.float32).copy()
                    tc_matrix.append(row)

            return {'lattices': lattices, 'q_values': q_values, 'tc_matrix': tc_matrix}
        except Exception as e:
            print(f"❌ Error leyendo análisis Tc {found}: {e}")
            return None
        

    @staticmethod
    def read_phase_transition_raw(filename: str) -> Optional[Dict[str, Any]]:
        """
        Lee archivo de transición con DATOS CRUDOS (series de tiempo) para calcular barras de error.
        
        Formato esperado:
        [n_T (uint32)] [n_samples (uint32)]
        [Vector T (float32, size=n_T)]
        [Matriz m (float32, size=n_T*n_samples)]
        [Matriz E (float32, size=n_T*n_samples)]
        """
        filepath = Path(filename)
        candidates = [filepath, filepath.with_suffix('.bin')]
        found = try_multiple_filenames(candidates)
        
        if not found:
            return None

        try:
            with open(found, 'rb') as f:
                # 1. Leer dimensiones
                n_T = struct.unpack('<I', read_exact(f, 4))[0]
                n_samples = struct.unpack('<I', read_exact(f, 4))[0]
                
                # 2. Leer vector de Temperaturas
                T = np.frombuffer(read_exact(f, n_T * 4), dtype=np.float32).copy()
                
                # Calcular tamaño total de las matrices planas
                total_elements = n_T * n_samples
                
                # 3. Leer Matriz de Magnetización (aplanada) y redimensionar
                m_flat = np.frombuffer(read_exact(f, total_elements * 4), dtype=np.float32)
                # Se asume que C++ guardó fila por fila: m_raw[i] es la historia de T[i]
                m_raw = m_flat.reshape((n_T, n_samples)).copy()
                
                # 4. Leer Matriz de Energía (aplanada) y redimensionar
                e_flat = np.frombuffer(read_exact(f, total_elements * 4), dtype=np.float32)
                e_raw = e_flat.reshape((n_T, n_samples)).copy()
            
            return {
                'T': T, 
                'm_raw': m_raw, 
                'E_raw': e_raw,
                'n_samples': n_samples
            }
            
        except Exception as e:
            # Si falla (por ejemplo, si intentamos leer un archivo antiguo con este lector),
            # simplemente devolvemos None y dejamos que el código use el lector antiguo.
            return None
        
    @staticmethod
    def read_paramagnetism_raw(filename: str) -> Optional[Dict[str, Any]]:
        """
        Lee archivo de paramagnetismo con DATOS CRUDOS (3D).
        Estructura: [Temp][Campo][Muestras]
        
        Header:
            n_T (uint32), n_H (uint32), n_samples (uint32)
        Data:
            T[] (n_T floats)
            H[] (n_H floats)
            m_raw[][][] (n_T * n_H * n_samples floats, aplanado)
        """
        filepath = Path(filename)
        candidates = [filepath, filepath.with_suffix('.bin')]
        found = try_multiple_filenames(candidates)
        
        if not found:
            return None

        try:
            with open(found, 'rb') as f:
                # 1. Leer dimensiones
                n_T = struct.unpack('<I', read_exact(f, 4))[0]
                n_H = struct.unpack('<I', read_exact(f, 4))[0]
                n_samples = struct.unpack('<I', read_exact(f, 4))[0]
                
                # 2. Leer ejes
                temperatures = np.frombuffer(read_exact(f, n_T * 4), dtype=np.float32).copy()
                fields = np.frombuffer(read_exact(f, n_H * 4), dtype=np.float32).copy()
                
                # 3. Leer la gran matriz de datos
                total_elements = n_T * n_H * n_samples
                raw_buffer = read_exact(f, total_elements * 4)
                m_flat = np.frombuffer(raw_buffer, dtype=np.float32)
                
                # Redimensionar a (n_T, n_H, n_samples)
                # C++ guardó: loop T -> loop H -> vector samples
                m_raw = m_flat.reshape((n_T, n_H, n_samples)).copy()
            
            return {
                'T': temperatures,
                'H': fields,
                'm_raw': m_raw,  # Matriz 3D
                'n_samples': n_samples
            }
            
        except Exception as e:
            # Si falla, devolvemos None para que main.py intente el método antiguo
            return None