#include <mpi.h>
#include <utility>
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
constexpr int upper_border_send_tag = 1;
constexpr int lower_border_send_tag = 2;


using namespace std;

class Grid {
public:
    vector<double> phi_;
    vector<double> phi_plus_one_;
    int volume_;
    int x_size_;
    int y_size_;
    int z_size_;
    int layer_number_;
    int layers_count_;
    double delta_ = 0.0f;
    double global_delta_ = std::numeric_limits<double>::infinity();;

    static double phi(vector<double> coords) {
        return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2];
    }

    static double rho(vector<double> coords) {
        return 6 - a * phi(std::move(coords));
    }


//public:
    vector<double> upper_neighbor_;
    vector<double> lower_neighbor_;

    int get_row_index_from_matrix_indices(int x, int y, int z) const {
        return x * this->y_size_ * this->z_size_ + z * this->y_size_ + y;
    }

    Grid(int x_size, int y_size, int z_size, int layer_number, int layers_count) {

        this->x_size_ = x_size;
        this->y_size_ = y_size;
        this->z_size_ = z_size;
        this->volume_ = this->x_size_ * this->y_size_ * this->z_size_;
        this->phi_.resize(this->volume_);
        this->phi_plus_one_.resize(this->volume_);
        this->layer_number_ = layer_number;
        this->layers_count_ = layers_count;
        this->lower_neighbor_.resize(this->x_size_ * this->y_size_);
        this->upper_neighbor_.resize(this->x_size_ * this->y_size_);
        std::fill(this->phi_.begin(), this->phi_.end(), 0);
        std::fill(this->phi_plus_one_.begin(), this->phi_plus_one_.end(), 0);
        std::fill(this->lower_neighbor_.begin(), this->lower_neighbor_.end(), 0);
        std::fill(this->upper_neighbor_.begin(), this->upper_neighbor_.end(), 0);
    }

    void set_borders_phi() {
        for (int k = 0; k < this->z_size_; ++k) {
            // заполняем верхнюю и нижнюю строку слоя
            for (int i = 0; i < this->x_size_; ++i) {
                this->phi_[this->get_row_index_from_matrix_indices(i, 0, k)] = phi(
                        this->get_coords_from_matrix_indices(i, 0, k));
                this->phi_[this->get_row_index_from_matrix_indices(i, this->y_size_ - 1, k)] = phi(
                        this->get_coords_from_matrix_indices(i, this->y_size_ - 1, k));
            }
            // заполняем левую и правую строку слоя
            for (int j = 0; j < this->y_size_; ++j) {
                this->phi_[this->get_row_index_from_matrix_indices(0, j, k)] = phi(
                        this->get_coords_from_matrix_indices(0, j, k));
                this->phi_[this->get_row_index_from_matrix_indices(this->x_size_ - 1, j, k)] = phi(
                        this->get_coords_from_matrix_indices(this->x_size_ - 1, j, k));
            }
        }
    }

    // this->layer_number_ - сколько слоев нужно пропустить
    // z_size_ - размер слоя
    // h_z - размер отступа в слое
    vector<double> get_coords_from_matrix_indices(int i, int j, int k) const {
        return vector<double>{x_0 + i * h_x, y_0 + j * h_y, z_0 + k * h_z + this->layer_number_ * z_size_ * h_z};
    }

    // получаем верхнюю границу
    vector<double> get_upper_border() const {
        vector<double> border(this->x_size_ * this->y_size_);
        for (int j = 0; j < this->y_size_; ++j) {
            for (int i = 0; i < this->x_size_; ++i) {
                border[i * this->y_size_ + j] = this->phi_[this->get_row_index_from_matrix_indices(i, j,
                                                                                                   this->z_size_ - 1)];
            }
        }
        return border;
    }

    // получаем нижнюю границу
    vector<double> get_lower_border() const {
        vector<double> border(this->x_size_ * this->y_size_);
        for (int i = 0; i < this->x_size_; ++i) {
            for (int j = 0; j < this->y_size_; ++j) {
                border[i * this->y_size_ + j] = this->phi_[this->get_row_index_from_matrix_indices(i, j, 0)];
            }
        }
        return border;
    }

    void iteration() {
        this->delta_ = 0.0f;
        MPI_Request upper_border_send_request, lower_border_send_request;
        if (layers_count_ != 1) {
            if (this->layer_number_ == ROOT) {
                // нижний слой должен отправить только верхнюю границу
                MPI_Isend(this->get_upper_border().data(), this->x_size_ * this->y_size_, MPI_DOUBLE, ROOT + 1,
                          upper_border_send_tag, MPI_COMM_WORLD, &upper_border_send_request);
            } else if (this->layer_number_ == layers_count_ - 1) {
                // верхний слой должен отправить только нижнюю границу
                MPI_Isend(this->get_lower_border().data(), this->x_size_ * this->y_size_, MPI_DOUBLE,
                          this->layer_number_ - 1, 1, MPI_COMM_WORLD, &lower_border_send_request);
            } else if (layers_count_ != 2) {
                // все остальные отправляют и верхнюю и нижнюю границу
                MPI_Isend(this->get_upper_border().data(), this->x_size_ * this->y_size_, MPI_DOUBLE,
                          this->layer_number_ + 1,
                          1, MPI_COMM_WORLD, &upper_border_send_request);
                MPI_Isend(this->get_lower_border().data(), this->x_size_ * this->y_size_, MPI_DOUBLE,
                          this->layer_number_ - 1, 1, MPI_COMM_WORLD, &lower_border_send_request);
            }

            if (this->layer_number_ == ROOT) {
                // нижний слой должен получить только верхнего соседа
                MPI_Irecv(this->upper_neighbor_.data(), this->x_size_ * this->y_size_, MPI_DOUBLE, ROOT + 1,
                          1, MPI_COMM_WORLD, &lower_border_send_request);
            } else if (this->layer_number_ == layers_count_ - 1) {
                // верхний слой должен получить только нижнего соседа
                MPI_Irecv(this->lower_neighbor_.data(), this->x_size_ * this->y_size_, MPI_DOUBLE,
                          this->layer_number_ - 1,
                          1, MPI_COMM_WORLD, &upper_border_send_request);
            } else if (layers_count_ != 2) {
                // все остальные получают обоих соседей
                MPI_Irecv(this->upper_neighbor_.data(), this->x_size_ * this->y_size_, MPI_DOUBLE,
                          this->layer_number_ + 1,
                          1, MPI_COMM_WORLD, &lower_border_send_request);
                MPI_Irecv(this->lower_neighbor_.data(), this->x_size_ * this->y_size_, MPI_DOUBLE,
                          this->layer_number_ - 1,
                          1, MPI_COMM_WORLD, &upper_border_send_request);
            }
        }
        // так как вычисление в центре не зависит от соседей, то итерируемся по центру
        this->center_iteration();
        // а вот вычисление на границе зависит от соседей, поэтому нужно дождаться получения этих соседей
        if (layers_count_ != 1) {
            if (layer_number_ == ROOT) {
                MPI_Wait(&lower_border_send_request, MPI_STATUS_IGNORE);
            } else if (this->layer_number_ == layers_count_ - 1) {
                MPI_Wait(&upper_border_send_request, MPI_STATUS_IGNORE);
            } else if (layers_count_ != 2){
                MPI_Wait(&lower_border_send_request, MPI_STATUS_IGNORE);
                MPI_Wait(&upper_border_send_request, MPI_STATUS_IGNORE);
            }
        }
        // дождались, итерируемся по границе
        this->borders_iteration();
        // считаем ∆
        this->update_delta();
        MPI_Allreduce(&this->delta_, &global_delta_, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Bcast(&global_delta_, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
        // новая фи становится старой
        swap(this->phi_, this->phi_plus_one_);
    }

    void update_delta() {
        for (int k = 0; k < this->z_size_; ++k) {
            for (int j = 0; j < this->y_size_; ++j) {
                for (int i = 0; i < this->x_size_; ++i) {
                    double current_delta = std::abs(
                            this->phi_plus_one_[this->get_row_index_from_matrix_indices(i, j, k)] -
                            this->phi_[this->get_row_index_from_matrix_indices(i, j, k)]);
                    this->delta_ = std::max(this->delta_, current_delta);
                }
            }
        }
    }

    double calculate_phi_plus_one(int i, int j, int k) {
        double tmp = (this->get_current_phi(i + 1, j, k) + this->get_current_phi(i - 1, j, k)) / (h_x * h_x) +
                     (this->get_current_phi(i, j + 1, k) + this->get_current_phi(i, j - 1, k)) / (h_y * h_y) +
                     (this->get_current_phi(i, j, k + 1) + this->get_current_phi(i, j, k - 1)) / (h_z * h_z) -
                     rho(this->get_coords_from_matrix_indices(i, j, k));
        return alpha * tmp;
    }

    double get_current_phi(int i, int j, int k) {
        // вылезли влево или вправо за слой
        if (i < 0 || j < 0 || i >= this->x_size_ || j >= this->y_size_ && k >= 0 && k < this->z_size_) {
            return phi(this->get_coords_from_matrix_indices(i, j, k));
        }
        // вылезли вниз из нижнего слоя
        if (k < 0 && this->layer_number_ == ROOT) {
            return phi(this->get_coords_from_matrix_indices(i, j, k));
        }
        // вылезли вверх из верхнего слоя
        if (k >= this->z_size_ && this->layer_number_ == layers_count_ - 1) {
            return phi(this->get_coords_from_matrix_indices(i, j, k));
        }
        // вылезли вниз из какого-то внутреннего слоя
        if (k < 0) {
            return this->lower_neighbor_[i * this->y_size_ + j];
        }
        // вылезли вверх из какого-то внутреннего слоя
        if (k >= this->z_size_) {
            return this->upper_neighbor_[i * this->y_size_ + j];
        }
        return this->phi_[this->get_row_index_from_matrix_indices(i, j, k)];
    }

    void center_iteration() {
        // идем по каждому слою отступая от краёв на один элемент и вычисляем фи
        for (int k = 0; k < this->z_size_; ++k) {
            for (int j = 1; j < this->y_size_ - 1; ++j) {
                for (int i = 1; i < this->x_size_ - 1; ++i) {
                    this->phi_plus_one_[this->get_row_index_from_matrix_indices(i, j,
                                                                                k)] = this->calculate_phi_plus_one(i, j,
                                                                                                                   k);
                }
            }
        }
    }

    void borders_iteration() {
        for (int k = 0; k < this->z_size_; ++k) {
            // заполняем верхнюю и нижнюю строку слоя
            for (int i = 0; i < this->x_size_; ++i) {
                this->phi_plus_one_[this->get_row_index_from_matrix_indices(i, 0, k)] =
                        this->calculate_phi_plus_one(i, 0, k);
                this->phi_plus_one_[this->get_row_index_from_matrix_indices(i, this->y_size_ - 1, k)] =
                        this->calculate_phi_plus_one(i, this->y_size_ - 1, k);
            }
            // заполняем левую и правую строку слоя
            for (int j = 0; j < this->y_size_; ++j) {
                this->phi_plus_one_[this->get_row_index_from_matrix_indices(0, j, k)] =
                        this->calculate_phi_plus_one(0, j, k);
                this->phi_plus_one_[this->get_row_index_from_matrix_indices(this->x_size_ - 1, j, k)] =
                        this->calculate_phi_plus_one(this->x_size_ - 1, j, k);
            }
        }
    }

    void check_result() {
        double this_delta = 0.0f;
        double current_delta;
        for (int k = 0; k < this->z_size_; ++k) {
            for (int j = 0; j < this->y_size_; ++j) {
                for (int i = 0; i < this->x_size_; ++i) {
                    current_delta = std::abs(this->phi_[this->get_row_index_from_matrix_indices(i, j, k)] -
                                             phi(this->get_coords_from_matrix_indices(i, j, k)));
                    this_delta = std::max(this_delta, current_delta);
                }
            }
        }
        MPI_Allreduce(&this_delta, &global_delta_, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (this->layer_number_ == ROOT) {
            cout << global_delta_ << endl;
        }
    }


    friend std::ostream &operator<<(std::ostream &os, const Grid &grid) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(5);
        for (int k = 0; k < grid.z_size_; k++) {
            os << "Layer #" << k + 1 << endl;
            for (int j = 0; j < grid.y_size_; j++) {
                for (int i = 0; i < grid.x_size_; i++) {
                    ss << '(' << grid.get_coords_from_matrix_indices(i, j, k)[0] << ","
                       << grid.get_coords_from_matrix_indices(i, j, k)[1]
                       << ", " << grid.get_coords_from_matrix_indices(i, j, k)[2] << ") "
                       << grid.phi_[grid.get_row_index_from_matrix_indices(i, j, k)];
                    os << ss.str() << ' ';
                    ss.str("");
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }
        return os;
    }

    double get_delta() {
        return global_delta_;
    }

};


int main(int argc, char **argv) {
    int world_rank, world_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    assert(N_z % world_size == 0);
    Grid grid(N_x, N_y, N_z / world_size, world_rank, world_size);
    grid.set_borders_phi();
    for (int i = 0; i < 5; ++i) {
        grid.iteration();
        if (world_rank == 0) {
            cout << grid.phi_.data() << ' ' << grid.phi_plus_one_.data() << endl;
        }
    }
    if (world_rank == 0) {
        cout << grid;
    }
    grid.check_result();
    MPI_Finalize();
    return 0;
}