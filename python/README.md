# Simulación Monte Carlo del Modelo de Ising en Diferentes Tipos de Redes

Implementación híbrida C++/Python para simulaciones Monte Carlo del modelo de Ising con dinámica de Metropolis, aplicado a múltiples topologías de red cristalina.

## 👤 Autor

**Jorge Steven Sarmiento Arboleda**  
Universidad de Antioquia, Medellín, Colombia  
Curso: Física Estadística (0302470)

## 🔬 Descripción

Este proyecto implementa simulaciones del modelo de Ising para estudiar:

- **Paramagnetismo** (J=0): Curvas m vs H, ley de estados correspondientes
- **Ferromagnetismo** (J=1): Histéresis magnética, transiciones de fase
- **Temperatura crítica**: Estimación de Tc para diferentes topologías y concentraciones

### Redes Soportadas

| Red | Dimensión | Coordinación (z) |
|-----|-----------|------------------|
| Cadena | 1D | 2 |
| Honeycomb (Panal) | 2D | 3 |
| Cuadrada | 2D | 4 |
| BCC | 3D | 8 |

## 📁 Estructura del Proyecto

```text
proyecto_ising_hibrido/
├── cpp/
│   ├── ising_model.hpp        # Definición de la clase IsingModel
│   ├── ising_model.cpp        # Implementación del modelo
│   ├── main.cpp               # Funciones de simulación y main
│   ├── ising_sim              # Ejecutable compilado
│   └── data/                  # Datos binarios generados
│       ├── paramagnetismo/
│       ├── ferromagnetismo/
│       ├── transicion/
│       ├── temperatura_critica/
│       ├── snapshots/
│       └── relax/
├── python/
│   ├── src/
│   │   ├── config.py          # Configuración y constantes
│   │   ├── reader.py          # Lectores de datos binarios
│   │   ├── plotter.py         # Sistema de graficado
│   │   ├── utils.py           # Funciones auxiliares
│   │   └── __init__.py
│   ├── main.py                # Punto de entrada para gráficas
│   └── graficas/              # Gráficas generadas
│       ├── paramagnetismo/
│       ├── ferromagnetismo/
│       ├── transicion/
│       ├── snapshots/
│       └── relajacion/
├── requirements.txt
├── README.md
└── venv/
```

## 🚀 Instalación y Uso

### 1. Compilar C++

```bash
cd cpp
g++ -o ising_sim main.cpp ising_model.cpp -std=c++17 -O2
```

### 2. Ejecutar Simulaciones

```bash
./ising_sim
```

Esto genera los datos binarios en `cpp/data/`.

### 3. Configurar Python

```bash
cd ../python
python3 -m venv venv
source venv/bin/activate  # Linux/Mac
pip install -r requirements.txt
```

### 4. Generar Gráficas

```bash
python main.py
```

Las gráficas se guardan en `python/graficas/`.

## 📦 Dependencias

**C++:**

- Compilador con soporte C++17 (g++ 7+ o clang++ 5+)
- Librería estándar (`<filesystem>`, `<random>`, `<vector>`)

**Python (`requirements.txt`):**

```text
numpy
matplotlib
scipy
```

## 📐 Sistema de Unidades

El código usa unidades naturales:

- $k_B = 1$
- $J = 1$ (constante de acoplamiento)
- Temperatura en unidades de $J/k_B$
- Campo magnético en unidades de $J$

## 🧮 Métodos Implementados

### Simulaciones (C++)

| Simulación | Descripción | Termalización |
|------------|-------------|---------------|
| Relajación | Evolución temporal de E desde estado inicial | No |
| Snapshots | Configuraciones de espines en equilibrio | Sí (escala con z) |
| Paramagnetismo | Curvas m vs H para J=0 | Mínima |
| Histéresis | Ciclos m vs H para J=1 | Solo al inicio |
| Transición de fase | Curvas m vs T, estimación de Tc | Sí (escala con z) |

### Análisis (Python)

- **Ley de Estados Correspondientes**: Colapso de datos m vs H/T
- **Ajuste Tanh**: $m = a \cdot \tanh(b \cdot H/T)$ con bondad de ajuste $\chi^2$
- **Susceptibilidad**: $\chi = \text{Var}(m)/T$ con errores por blocking
- **Temperatura Crítica**: Estimada por máxima pendiente de m(T)

## 📊 Ejemplos de Resultados

### Paramagnetismo

- Curvas m vs H para diferentes T
- Verificación de $m = \tanh(H/T)$

### Ferromagnetismo

- Ciclos de histéresis con área no nula
- Transición de fase con identificación de Tc

### Efecto de la Topología

- Mayor z → Mayor Tc
- Mayor z → Histéresis más pronunciada

## 📝 Notas Técnicas

- Datos almacenados en formato binario little-endian (float32)
- Condiciones de frontera periódicas en todas las redes
- Algoritmo de Metropolis con selección aleatoria de sitios
- Errores estadísticos calculados por método de blocking

## 📄 Licencia

MIT License - Ver archivo LICENSE para más detalles.

---

**Versión**: 1.0  
**Fecha**: 2025
