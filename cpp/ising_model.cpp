#include "ising_model.hpp"
#include <algorithm>
#include <numeric>
#include <iostream>
#include <unordered_map>

// ====================================================================
// CONSTRUCTOR
// ====================================================================

IsingModel::IsingModel(int L_, float q_, const std::string& lattice_str, 
                       float J_, uint32_t seed)
    : L(L_), q(q_), J(J_), T(1.0f), H(0.0f), beta(1.0f),
      uniform_dist(0.0f, 1.0f), spin_dist(0, 1) {
    
    validate_parameters();
    
    if (lattice_str == "chain") {
        lattice_type = LatticeType::CHAIN;
        dimension = 1;
        z = 2;
    } else if (lattice_str == "square") {
        lattice_type = LatticeType::SQUARE;
        dimension = 2;
        z = 4;
    } else if (lattice_str == "honeycomb") {
        lattice_type = LatticeType::HONEYCOMB;
        dimension = 2;
        z = 3;
    } else if (lattice_str == "bcc") {
        lattice_type = LatticeType::BCC;
        dimension = 3;
        z = 8;
    } else {
        throw std::invalid_argument("Unknown lattice_type: " + lattice_str);
    }
    
    setup_rng(seed);
    apply_dilution();
    initialize_random();
}

void IsingModel::validate_parameters() {
    if (L <= 0) throw std::invalid_argument("L must be positive");
    if (q <= 0.0f || q > 1.0f) throw std::invalid_argument("q must be in (0, 1]");
    if (J < 0.0f) throw std::invalid_argument("J must be non-negative");
}

void IsingModel::setup_rng(uint32_t seed) {
    if (seed == 0) {
        std::random_device rd;
        rng.seed(rd());
    } else {
        rng.seed(seed);
    }
}

// ====================================================================
// DILUCIÓN: SELECTIONAR SITIOS OCUPADOS ALEATORIAMENTE
// ====================================================================

void IsingModel::apply_dilution() {
    if (dimension == 1) {
        apply_dilution_chain();
    } else if (dimension == 2) {
        apply_dilution_square();
    } else {
        apply_dilution_bcc();
    }
    
    // Inicializar array de spins con tamaño = sitios ocupados
    spin.assign(occupied_sites_2d.size() + occupied_sites_3d.size(), 0);
}

void IsingModel::apply_dilution_chain() {
    int n_total = L;
    int n_occupied = static_cast<int>(q * n_total);
    
    std::vector<int> all_sites(n_total);
    std::iota(all_sites.begin(), all_sites.end(), 0);
    std::shuffle(all_sites.begin(), all_sites.end(), rng);
    
    occupied_sites_2d.clear();
    for (int i = 0; i < n_occupied; ++i) {
        occupied_sites_2d.push_back({all_sites[i], 0});
    }
}

void IsingModel::apply_dilution_square() {
    int n_total = L * L;
    int n_occupied = static_cast<int>(q * n_total);
    
    std::vector<Site> all_sites;
    all_sites.reserve(n_total);
    
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            all_sites.push_back({i, j});
        }
    }
    
    std::shuffle(all_sites.begin(), all_sites.end(), rng);
    occupied_sites_2d.assign(all_sites.begin(), all_sites.begin() + n_occupied);
}

void IsingModel::apply_dilution_bcc() {
    int n_total = L * L * L;
    int n_occupied = static_cast<int>(q * n_total);
    
    std::vector<Site3D> all_sites;
    all_sites.reserve(n_total);
    
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            for (int k = 0; k < L; ++k) {
                all_sites.push_back({i, j, k});
            }
        }
    }
    
    std::shuffle(all_sites.begin(), all_sites.end(), rng);
    occupied_sites_3d.assign(all_sites.begin(), all_sites.begin() + n_occupied);
}

// ====================================================================
// PARÁMETROS
// ====================================================================

void IsingModel::set_T(float temperature) {
    if (temperature <= 0.0f) throw std::invalid_argument("Temperature must be positive");
    if (std::isinf(temperature) || std::isnan(temperature))
        throw std::invalid_argument("Invalid temperature value");
    T = temperature;
    beta = 1.0f / T;
}

void IsingModel::set_H(float magnetic_field) {
    if (std::isinf(magnetic_field) || std::isnan(magnetic_field))
        throw std::invalid_argument("Invalid magnetic field value");
    H = magnetic_field;
}

// ====================================================================
// INICIALIZACIÓN
// ====================================================================

void IsingModel::initialize_random() {
    for (size_t i = 0; i < spin.size(); ++i) {
        spin_at(i) = (spin_dist(rng) == 0) ? -1 : 1;
    }
}

void IsingModel::initialize_ordered(int spin_value) {
    if (spin_value != -1 && spin_value != 1)
        throw std::invalid_argument("spin_value must be ±1");
    
    for (size_t i = 0; i < spin.size(); ++i) {
        spin_at(i) = spin_value;
    }
}

// ====================================================================
// VECINOS
// ====================================================================

void IsingModel::get_neighbors_chain(int i, std::vector<int>& neighbors) const {
    neighbors.clear();
    neighbors.push_back((i - 1 + L) % L);
    neighbors.push_back((i + 1) % L);
}

void IsingModel::get_neighbors_square(int i, int j, std::vector<Site>& neighbors) const {
    neighbors.clear();
    neighbors.push_back({(i - 1 + L) % L, j});
    neighbors.push_back({(i + 1) % L, j});
    neighbors.push_back({i, (j - 1 + L) % L});
    neighbors.push_back({i, (j + 1) % L});
}

void IsingModel::get_neighbors_honeycomb(int i, int j, std::vector<Site>& neighbors) const {
    neighbors.clear();
    
    // Red hexagonal bipartita con z=3 vecinos por sitio
    // Determinado por paridad (i+j) mod 2
    if ((i + j) % 2 == 0) {  // Sublattice A
        neighbors.push_back({i, (j + 1) % L});
        neighbors.push_back({(i - 1 + L) % L, j});
        neighbors.push_back({(i - 1 + L) % L, (j + 1) % L});
    } else {  // Sublattice B
        neighbors.push_back({i, (j - 1 + L) % L});
        neighbors.push_back({(i + 1) % L, j});
        neighbors.push_back({(i + 1) % L, (j - 1 + L) % L});
    }
}

void IsingModel::get_neighbors_bcc(int i, int j, int k, std::vector<Site3D>& neighbors) const {
    neighbors.clear();
    
    // BCC: 8 vecinos diagonales
    int di[] = {-1, -1, -1, -1, 1, 1, 1, 1};
    int dj[] = {-1, -1, 1, 1, -1, -1, 1, 1};
    int dk[] = {-1, 1, -1, 1, -1, 1, -1, 1};
    
    for (int n = 0; n < 8; ++n) {
        neighbors.push_back({
            (i + di[n] + L) % L,
            (j + dj[n] + L) % L,
            (k + dk[n] + L) % L
        });
    }
}

// ====================================================================
// ENERGÍA Y MAGNETIZACIÓN
// ====================================================================

float IsingModel::energy() const {
    float E_interaction = 0.0f;
    float E_field = 0.0f;
    
    std::vector<Site> neighbors_2d;
    std::vector<Site3D> neighbors_3d;
    
    // Mapa: coordenada global -> índice local
    std::unordered_map<int64_t, int> coord_to_local;
    
    if (dimension <= 2) {
        for (size_t local_idx = 0; local_idx < occupied_sites_2d.size(); ++local_idx) {
            const auto& site = occupied_sites_2d[local_idx];
            int64_t key = ((int64_t)site.first << 32) | (site.second & 0xFFFFFFFF);
            coord_to_local[key] = local_idx;
        }
    } else {
        for (size_t local_idx = 0; local_idx < occupied_sites_3d.size(); ++local_idx) {
            const auto& site = occupied_sites_3d[local_idx];
            int64_t key = ((int64_t)site[0] << 40) | ((int64_t)site[1] << 20) | site[2];
            coord_to_local[key] = local_idx;
        }
    }
    
    if (dimension == 1) {
        for (size_t local_idx = 0; local_idx < occupied_sites_2d.size(); ++local_idx) {
            int i = occupied_sites_2d[local_idx].first;
            int8_t spin_i = spin_at(local_idx);
            
            std::vector<int> neighbors_1d;
            get_neighbors_chain(i, neighbors_1d);
            
            for (int n_i : neighbors_1d) {
                int64_t key = ((int64_t)n_i << 32);
                auto it = coord_to_local.find(key);
                if (it != coord_to_local.end()) {
                    int8_t spin_j = spin_at(it->second);
                    E_interaction -= J * spin_i * spin_j;
                }
            }
            
            E_field -= H * spin_i;
        }
    } else if (dimension == 2) {
        for (size_t local_idx = 0; local_idx < occupied_sites_2d.size(); ++local_idx) {
            int i = occupied_sites_2d[local_idx].first;
            int j = occupied_sites_2d[local_idx].second;
            int8_t spin_i = spin_at(local_idx);
            
            neighbors_2d.clear();
            if (lattice_type == LatticeType::SQUARE) {
                get_neighbors_square(i, j, neighbors_2d);
            } else {
                get_neighbors_honeycomb(i, j, neighbors_2d);
            }
            
            for (const auto& n : neighbors_2d) {
                int64_t key = ((int64_t)n.first << 32) | (n.second & 0xFFFFFFFF);
                auto it = coord_to_local.find(key);
                if (it != coord_to_local.end()) {
                    int8_t spin_j = spin_at(it->second);
                    E_interaction -= J * spin_i * spin_j;
                }
            }
            
            E_field -= H * spin_i;
        }
    } else {
        for (size_t local_idx = 0; local_idx < occupied_sites_3d.size(); ++local_idx) {
            int i = occupied_sites_3d[local_idx][0];
            int j = occupied_sites_3d[local_idx][1];
            int k = occupied_sites_3d[local_idx][2];
            int8_t spin_i = spin_at(local_idx);
            
            get_neighbors_bcc(i, j, k, neighbors_3d);
            
            for (const auto& n : neighbors_3d) {
                int64_t key = ((int64_t)n[0] << 40) | ((int64_t)n[1] << 20) | n[2];
                auto it = coord_to_local.find(key);
                if (it != coord_to_local.end()) {
                    int8_t spin_j = spin_at(it->second);
                    E_interaction -= J * spin_i * spin_j;
                }
            }
            
            E_field -= H * spin_i;
        }
    }
    
    E_interaction /= 2.0f;  // Cada par contado dos veces
    return E_interaction + E_field;
}

// ====================================================================
// METROPOLIS
// ====================================================================

inline void IsingModel::metropolis_step() {
    std::vector<size_t> indices(spin.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);
    
    std::vector<Site> neighbors_2d;
    std::vector<Site3D> neighbors_3d;
    
    // Mapa para búsqueda rápida
    std::unordered_map<int64_t, int> coord_to_local;
    
    if (dimension <= 2) {
        for (size_t local_idx = 0; local_idx < occupied_sites_2d.size(); ++local_idx) {
            const auto& site = occupied_sites_2d[local_idx];
            int64_t key = ((int64_t)site.first << 32) | (site.second & 0xFFFFFFFF);
            coord_to_local[key] = local_idx;
        }
    } else {
        for (size_t local_idx = 0; local_idx < occupied_sites_3d.size(); ++local_idx) {
            const auto& site = occupied_sites_3d[local_idx];
            int64_t key = ((int64_t)site[0] << 40) | ((int64_t)site[1] << 20) | site[2];
            coord_to_local[key] = local_idx;
        }
    }
    
    for (size_t local_idx : indices) {
        int8_t spin_i = spin_at(local_idx);
        float sum_neighbors = 0.0f;
        
        if (dimension == 1) {
            int i = occupied_sites_2d[local_idx].first;
            std::vector<int> neighbors_1d;
            get_neighbors_chain(i, neighbors_1d);
            
            for (int n_i : neighbors_1d) {
                int64_t key = ((int64_t)n_i << 32);
                auto it = coord_to_local.find(key);
                if (it != coord_to_local.end()) {
                    sum_neighbors += spin_at(it->second);
                }
            }
        } else if (dimension == 2) {
            int i = occupied_sites_2d[local_idx].first;
            int j = occupied_sites_2d[local_idx].second;
            
            neighbors_2d.clear();
            if (lattice_type == LatticeType::SQUARE) {
                get_neighbors_square(i, j, neighbors_2d);
            } else {
                get_neighbors_honeycomb(i, j, neighbors_2d);
            }
            
            for (const auto& n : neighbors_2d) {
                int64_t key = ((int64_t)n.first << 32) | (n.second & 0xFFFFFFFF);
                auto it = coord_to_local.find(key);
                if (it != coord_to_local.end()) {
                    sum_neighbors += spin_at(it->second);
                }
            }
        } else {
            int i = occupied_sites_3d[local_idx][0];
            int j = occupied_sites_3d[local_idx][1];
            int k = occupied_sites_3d[local_idx][2];
            
            get_neighbors_bcc(i, j, k, neighbors_3d);
            
            for (const auto& n : neighbors_3d) {
                int64_t key = ((int64_t)n[0] << 40) | ((int64_t)n[1] << 20) | n[2];
                auto it = coord_to_local.find(key);
                if (it != coord_to_local.end()) {
                    sum_neighbors += spin_at(it->second);
                }
            }
        }
        
        float dE = 2.0f * spin_i * (J * sum_neighbors + H);
        
        if (dE <= 0.0f || uniform_dist(rng) < std::exp(-beta * dE)) {
            spin_at(local_idx) *= -1;
        }
    }
}

void IsingModel::thermalize(int n_steps) {
    for (int step = 0; step < n_steps; ++step) {
        metropolis_step();
    }
}

// ====================================================================
// GETTERS
// ====================================================================

std::string IsingModel::get_lattice_type_str() const {
    switch (lattice_type) {
        case LatticeType::CHAIN: return "chain";
        case LatticeType::SQUARE: return "square";
        case LatticeType::HONEYCOMB: return "honeycomb";
        case LatticeType::BCC: return "bcc";
        default: return "unknown";
    }
}

size_t IsingModel::get_n_occupied() const {
    return (dimension == 3) ? occupied_sites_3d.size() : occupied_sites_2d.size();
}

size_t IsingModel::get_n_total() const {
    if (dimension == 1) return L;
    if (dimension == 2) return L * L;
    return L * L * L;
}

std::vector<int> IsingModel::get_shape() const {
    if (dimension == 1) return {L};
    if (dimension == 2) return {L, L};
    return {L, L, L};
}