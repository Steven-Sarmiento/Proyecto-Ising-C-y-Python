#ifndef ISING_MODEL_HPP
#define ISING_MODEL_HPP

#include <vector>
#include <array>
#include <string>
#include <cstdint>
#include <random>
#include <cmath>
#include <stdexcept>
#include <algorithm>

/**
 * @class IsingModel
 * @brief Simulación optimizada del modelo de Ising con dilución magnética
 * 
 * Topologías soportadas:
 * - chain (1D): z=2
 * - square (2D): z=4
 * - honeycomb (2D): z=3
 * - bcc (3D): z=8
 * 
 * DISEÑO: Solo se almacenan spins OCUPADOS en array contiguo.
 * occupied_sites mapea índice local -> coordenada global.
 */
class IsingModel {
public:
    using Site = std::pair<int, int>;
    using Site3D = std::array<int, 3>;
    
    enum class LatticeType { CHAIN, SQUARE, HONEYCOMB, BCC };
    
    IsingModel(int L, float q, const std::string& lattice_str, float J, uint32_t seed = 0);
    
    IsingModel(const IsingModel&) = delete;
    IsingModel& operator=(const IsingModel&) = delete;
    ~IsingModel() = default;
    
private:
    // Parámetros del modelo
    int L;
    float q;
    LatticeType lattice_type;
    float J, T, H, beta;
    int dimension, z;
    
    // Spins: solo sitios OCUPADOS (array contiguo eficiente)
    std::vector<int8_t> spin;
    
    // Mapeo: índice local -> coordenada global
    // Para 1D: {i}
    // Para 2D: {i, j}
    // Para 3D: {i, j, k}
    std::vector<Site> occupied_sites_2d;
    std::vector<Site3D> occupied_sites_3d;
    
    // RNG
    std::mt19937 rng;
    std::uniform_real_distribution<float> uniform_dist;
    std::uniform_int_distribution<int> spin_dist;
    
    // Métodos privados
    void validate_parameters();
    void setup_rng(uint32_t seed);
    void apply_dilution();
    void apply_dilution_chain();
    void apply_dilution_square();
    void apply_dilution_bcc();
    
    void get_neighbors_chain(int i, std::vector<int>& neighbors) const;
    void get_neighbors_square(int i, int j, std::vector<Site>& neighbors) const;
    void get_neighbors_honeycomb(int i, int j, std::vector<Site>& neighbors) const;
    void get_neighbors_bcc(int i, int j, int k, std::vector<Site3D>& neighbors) const;
    
public:
    void set_T(float temperature);
    void set_H(float magnetic_field);
    void initialize_random();
    void initialize_ordered(int spin_value = 1);
    
    void metropolis_step();
    void thermalize(int n_steps);
    
    inline int8_t& spin_at(size_t local_idx) { return spin[local_idx]; }
    inline const int8_t& spin_at(size_t local_idx) const { return spin[local_idx]; }
    
    float magnetization() const {
        if (spin.empty()) return 0.0f;
        int sum_spins = 0;
        for (int8_t s : spin) sum_spins += s;
        return static_cast<float>(sum_spins) / spin.size();
    }
    
    float total_magnetization() const {
        int sum_spins = 0;
        for (int8_t s : spin) sum_spins += s;
        return static_cast<float>(sum_spins);
    }
    
    float energy() const;
    
    int get_L() const { return L; }
    float get_q() const { return q; }
    float get_J() const { return J; }
    float get_T() const { return T; }
    float get_H() const { return H; }
    int get_z() const { return z; }
    int get_dimension() const { return dimension; }
    std::string get_lattice_type_str() const;
    
    size_t get_n_occupied() const;
    size_t get_n_total() const;
    
    std::vector<int8_t> get_configuration() const { return spin; }
    std::vector<int> get_shape() const;
    
    // Para debugging
    const std::vector<Site>& get_occupied_2d() const { return occupied_sites_2d; }
    const std::vector<Site3D>& get_occupied_3d() const { return occupied_sites_3d; }
};

#endif // ISING_MODEL_HPP