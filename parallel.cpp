#include "/opt/homebrew/Cellar/mpich/4.1/include/mpi.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>

constexpr int ROOT = 0;
constexpr double x_0 = -1.0f;
constexpr double y_0 = -1.0f;
constexpr double z_0 = -1.0f;
constexpr double Dx = 2.0f;
constexpr double Dy = 2.0f;
constexpr double Dz = 2.0f;
constexpr int N_x = 4;
constexpr int N_y = 4;
constexpr int N_z = 4;
constexpr double a = 1e5;
constexpr double epsilon = 1e-8;
constexpr double h_x = Dx / (N_x - 1);
constexpr double h_y = Dy / (N_y - 1);
constexpr double h_z = Dz / (N_z - 1);
constexpr double alpha = 1.0f / ((2.0f / h_x * h_x) + (2.0f / h_y * h_y) + (2.0f / h_z * h_z) + a);


using namespace std;

class Grid {
private:
    vector<double> phi_;
    vector<double> phi_plus_one_;
    vector<double> &current_phi_ = phi_;
    int volume_;
    int x_size_;
    int y_size_;
    int z_size_;
    int layer_number_;
    int layers_count_;

    static double phi(vector<double> coords) {
        return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2];
    }

    static double rho(vector<double> coords) {
        return 6 - a * phi(coords);
    }

    int get_row_index(int x, int y, int z) const {
        return x * this->y_size_ * this->z_size_ + z * this->y_size_ + y;
    }


public:
    Grid(int x_size, int y_size, int z_size, int layer_number, int layers_count) {
        this->x_size_ = x_size;
        this->y_size_ = y_size;
        this->z_size_ = z_size;
        this->volume_ = this->x_size_ * this->y_size_ * this->z_size_;
        this->phi_.resize(this->volume_);
        this->phi_plus_one_.resize(this->volume_);
        this->layer_number_ = layer_number;
        this->layers_count_ = layers_count;
        std::fill(this->phi_.begin(), this->phi_.end(), 0);
        std::fill(this->phi_plus_one_.begin(), this->phi_plus_one_.end(), 0);
    }

    void set_borders_phi() {
        if (this->layer_number_ == ROOT) {
            // заполняем верхнюю и нижнюю строку нижнего слоя
            for (int i = 0; i < this->x_size_; ++i) {
                this->phi_[this->get_row_index(i, 0, 0)] = phi(this->get_coords_from_matrix_indices(i, 0, 0));
                this->phi_[this->get_row_index(i, this->y_size_ - 1, 0)] = phi(
                        this->get_coords_from_matrix_indices(i, this->y_size_ - 1, 0));
            }
            // заполняем левую и правую строку нижнего слоя
            for (int j = 0; j < this->y_size_; ++j) {
                this->phi_[this->get_row_index(0, j, 0)] = phi(this->get_coords_from_matrix_indices(0, j, 0));
                this->phi_[this->get_row_index(this->x_size_ - 1, j, 0)] = phi(
                        this->get_coords_from_matrix_indices(this->x_size_ - 1, j, 0));
            }
        }
        if (this->layer_number_ == layers_count_ - 1) {
            // заполняем верхнюю и нижнюю строку верхнего слоя
            for (int i = 0; i < this->x_size_; ++i) {
                this->phi_[this->get_row_index(i, 0, this->z_size_ - 1)] = phi(
                        this->get_coords_from_matrix_indices(i, 0, this->z_size_ - 1));
                this->phi_[this->get_row_index(i, this->y_size_ - 1, this->z_size_ - 1)] = phi(
                        this->get_coords_from_matrix_indices(i, this->y_size_ - 1, this->z_size_ - 1));
            }
            // заполняем левую и правую строку верхнего слоя
            for (int j = 0; j < this->y_size_; ++j) {
                this->phi_[this->get_row_index(0, j, this->z_size_ - 1)] = phi(
                        this->get_coords_from_matrix_indices(0, j, this->z_size_ - 1));
                this->phi_[this->get_row_index(this->x_size_ - 1, j, this->z_size_ - 1)] = phi(
                        this->get_coords_from_matrix_indices(this->x_size_ - 1, j, this->z_size_ - 1));
            }
        }
        // заполняем углы
        for (int k = 0; k < this->z_size_; ++k) {
            this->phi_[this->get_row_index(0, 0, k)] = phi(this->get_coords_from_matrix_indices(0, 0, k));
            this->phi_[this->get_row_index(0, this->y_size_ - 1, k)] = phi(
                    this->get_coords_from_matrix_indices(0, this->y_size_ - 1, k));
            this->phi_[this->get_row_index(this->x_size_ - 1, 0, k)] = phi(
                    this->get_coords_from_matrix_indices(this->x_size_ - 1, 0, k));
            this->phi_[this->get_row_index(this->x_size_ - 1, this->y_size_ - 1, k)] = phi(
                    this->get_coords_from_matrix_indices(this->x_size_ - 1, this->y_size_ - 1, k));
        }

    }

    // this->layer_number_ - сколько слоев нужно пропустить
    // z_size_ - размер слоя
    // h_z - размер отступа в слое
    vector<double> get_coords_from_matrix_indices(int i, int j, int k) const {
        return vector<double>{x_0 + i * h_x, y_0 + j * h_y, z_0 + k * h_z + this->layer_number_ * z_size_ * h_z};
    }

    friend std::ostream &operator<<(std::ostream &os, const Grid grid) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(3);
        for (int k = 0; k < grid.z_size_; k++) {
            os << "Layer #" << k + 1 << endl;
            for (int j = 0; j < grid.y_size_; j++) {
                for (int i = 0; i < grid.x_size_; i++) {
                    ss << '(' << grid.get_coords_from_matrix_indices(i, j, k)[0] << ","
                       << grid.get_coords_from_matrix_indices(i, j, k)[1]
                       << ", " << grid.get_coords_from_matrix_indices(i, j, k)[2] << ") "
                       << grid.phi_[grid.get_row_index(i, j, k)];
                    os << ss.str() << ' ';
                    ss.str("");
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }
        return os;
    }

};


int main(int argc, char **argv) {
    int world_rank, world_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    assert(N_z % world_size == 0);
    Grid grid(N_x, N_y, N_z / world_size, world_rank, world_size);
//    cout << "rank: " << world_rank << endl;
    grid.set_borders_phi();
    if (world_rank != ROOT)
        cout << grid;
    MPI_Finalize();
    return 0;
}