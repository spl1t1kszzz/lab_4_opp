#include <mpi.h>
#include <utility>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>

constexpr int ROOT = 0;
constexpr double x_0 = -1.0f;
constexpr double y_0 = -1.0f;
constexpr double z_0 = -1.0f;
constexpr double Dx = 2.0f;
constexpr double Dy = 2.0f;
constexpr double Dz = 2.0f;
constexpr int N_x = 36;
constexpr int N_y = 36;
constexpr int N_z = 36;
constexpr double a = 100000.0;
constexpr double epsilon = 1e-8;
constexpr double h_x = Dx / (N_x - 1);
constexpr double h_y = Dy / (N_y - 1);
constexpr double h_z = Dz / (N_z - 1);
constexpr double alpha = 1.0f / ((2.0f / (h_x * h_x)) + (2.0f / (h_y * h_y)) + (2.0f / (h_z * h_z)) + a);


using namespace std;

class Grid {
private:
    vector<double> phi_;
    vector<double> phi_plus_one_;
    int volume_;
    int x_size_;
    int y_size_;
    int z_size_;
    int layer_number_;
    int layers_count_;
    double global_delta_ = std::numeric_limits<double>::infinity();
    vector<double> upper_neighbor_;
    vector<double> lower_neighbor_;

    static double phi(vector<double> coords) {
        return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2];
    }

    double rho(int i, int j, int k) const {
        return 6 - a * phi(this->get_coords_from_matrix_indices(i, j, k));
    }

public:

    int get_row_index_from_matrix_indices(int x, int y, int z) const {
        return x * this->y_size_ * this->z_size_ + y * this->z_size_ + z;
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
            // заполняем заднюю и переднюю строку слоя
            for (int i = 0; i < this->x_size_; ++i) {
                this->phi_[this->get_row_index_from_matrix_indices(i, 0, k)] = phi(
                        this->get_coords_from_matrix_indices(i, 0, k));
                this->phi_[this->get_row_index_from_matrix_indices(i, this->y_size_ - 1, k)] = phi(
                        this->get_coords_from_matrix_indices(i, this->y_size_ - 1, k));
                this->phi_plus_one_[this->get_row_index_from_matrix_indices(i, 0,
                                                                            k)] = this->phi_[this->get_row_index_from_matrix_indices(
                        i, 0, k)];
                this->phi_plus_one_[this->get_row_index_from_matrix_indices(i, this->y_size_ - 1,
                                                                            k)] = this->phi_[this->get_row_index_from_matrix_indices(
                        i, this->y_size_ - 1, k)];

            }
            // заполняем левую и правую строку слоя
            for (int j = 0; j < this->y_size_; ++j) {
                this->phi_[this->get_row_index_from_matrix_indices(0, j, k)] = phi(
                        this->get_coords_from_matrix_indices(0, j, k));
                this->phi_plus_one_[this->get_row_index_from_matrix_indices(0, j,
                                                                            k)] = this->phi_[this->get_row_index_from_matrix_indices(
                        0, j, k)];
                this->phi_[this->get_row_index_from_matrix_indices(this->x_size_ - 1, j, k)] = phi(
                        this->get_coords_from_matrix_indices(this->x_size_ - 1, j, k));
                this->phi_plus_one_[this->get_row_index_from_matrix_indices(this->x_size_ - 1, j,
                                                                            k)] = this->phi_[this->get_row_index_from_matrix_indices(
                        this->x_size_ - 1, j, k)];

            }
        }
        for (int i = 0; i < this->x_size_; ++i) {
            for (int j = 0; j < this->y_size_; ++j) {
                if (layer_number_ == ROOT) {
                    this->phi_[this->get_row_index_from_matrix_indices(i, j, 0)] = phi(
                            this->get_coords_from_matrix_indices(i, j, 0));
                    this->phi_plus_one_[this->get_row_index_from_matrix_indices(i, j, 0)] = phi(
                            this->get_coords_from_matrix_indices(i, j, 0));
                }
                if (layer_number_ == layers_count_ - 1) {
                    this->phi_[this->get_row_index_from_matrix_indices(i, j, this->z_size_ - 1)] = phi(
                            this->get_coords_from_matrix_indices(i, j, this->z_size_ - 1));
                    this->phi_plus_one_[this->get_row_index_from_matrix_indices(i, j, this->z_size_ - 1)] = phi(
                            this->get_coords_from_matrix_indices(i, j, this->z_size_ - 1));
                }
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
        double global_delta = std::numeric_limits<double>::infinity();
        while (global_delta >= epsilon) {
            MPI_Request upper_border_send_request, lower_border_send_request, upper_border_recv_request, lower_border_recv_request;
            if (layers_count_ != 1) {
                if (this->layer_number_ == ROOT) {
                    // нижний слой должен отправить только верхнюю границу
                    MPI_Isend(this->get_upper_border().data(), this->x_size_ * this->y_size_, MPI_DOUBLE, ROOT + 1,
                              1, MPI_COMM_WORLD, &upper_border_send_request);
                } else if (this->layer_number_ == (layers_count_ - 1)) {
                    // верхний слой должен отправить только нижнюю границу
                    MPI_Isend(this->get_lower_border().data(), this->x_size_ * this->y_size_, MPI_DOUBLE,
                              this->layer_number_ - 1, 1, MPI_COMM_WORLD, &lower_border_send_request);
                } else {
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
                              1, MPI_COMM_WORLD, &upper_border_recv_request);
                } else if (this->layer_number_ == (layers_count_ - 1)) {
                    // верхний слой должен получить только нижнего соседа
                    MPI_Irecv(this->lower_neighbor_.data(), this->x_size_ * this->y_size_, MPI_DOUBLE,
                              this->layer_number_ - 1,
                              1, MPI_COMM_WORLD, &lower_border_recv_request);
                } else {
                    // все остальные получают обоих соседей
                    MPI_Irecv(this->upper_neighbor_.data(), this->x_size_ * this->y_size_, MPI_DOUBLE,
                              this->layer_number_ + 1,
                              1, MPI_COMM_WORLD, &upper_border_recv_request);
                    MPI_Irecv(this->lower_neighbor_.data(), this->x_size_ * this->y_size_, MPI_DOUBLE,
                              this->layer_number_ - 1,
                              1, MPI_COMM_WORLD, &lower_border_recv_request);
                }
            }
            // так как вычисление в центре не зависит от соседей, то итерируемся по центру
            this->center_iteration();
            // а вот вычисление на границе зависит от соседей, поэтому нужно дождаться получения этих соседей
            if (layers_count_ > 1) {
                if (layer_number_ == ROOT) {
                    MPI_Wait(&upper_border_recv_request, MPI_STATUS_IGNORE);
                } else if (layer_number_ == (layers_count_ - 1)) {
                    MPI_Wait(&lower_border_recv_request, MPI_STATUS_IGNORE);
                } else {
                    MPI_Wait(&lower_border_recv_request, MPI_STATUS_IGNORE);
                    MPI_Wait(&upper_border_recv_request, MPI_STATUS_IGNORE);
                }
            }
            // дождались, итерируемся по границе
            this->borders_iteration();
            // считаем ∆
            double local_delta = calculate_local_delta();

            // собираем максимальное значение delta в главном процессе
            global_delta = calculate_global_delta(local_delta);

            cout << layer_number_ << ": " << global_delta << endl;

            // новая фи становится старой
            std::swap(this->phi_, this->phi_plus_one_);
        }
    }

    double calculate_local_delta() {
        double local_delta = 0.0;
        for (int i = 0; i < volume_; ++i) {
            double current_delta = std::fabs(phi_plus_one_[i] - phi_[i]);
            local_delta = std::max(local_delta, current_delta);
        }
        return local_delta;
    }

    static double calculate_global_delta(double local_delta) {
        double global_delta;
        MPI_Allreduce(&local_delta, &global_delta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return global_delta;
    }

    double calculate_phi_plus_one(int i, int j, int k) {
        double ro = rho(i, j, k);
        double right = this->get_current_phi(i + 1, j, k);
        double left = this->get_current_phi(i - 1, j, k);
        double upY = this->get_current_phi(i, j + 1, k);
        double downY = this->get_current_phi(i, j - 1, k);
        double upZ = this->get_current_phi(i, j, k + 1);
        double downZ = this->get_current_phi(i, j, k - 1);
        double tmp = (right + left) / (h_x * h_x) +
                     (upY + downY) / (h_y * h_y) +
                     (upZ + downZ) / (h_z * h_z) -
                     ro;
        return alpha * tmp;
    }

    double get_current_phi(int i, int j, int k) {
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
        // идем по каждому слою отступая от краёв на два элемента и вычисляем фи
        for (int k = 1; k < this->z_size_ - 1; ++k) {
            for (int j = 1; j < this->y_size_ - 1; ++j) {
                for (int i = 1; i < this->x_size_ - 1; ++i) {
                    double var = this->calculate_phi_plus_one(i, j, k);
                    this->phi_plus_one_[this->get_row_index_from_matrix_indices(i, j, k)] = var;
                }
            }
        }

    }

    void borders_iteration() {
        for (int j = 1; j < this->y_size_ - 1; ++j) {
            for (int i = 1; i < this->x_size_ - 1; ++i) {
                if (layer_number_ != ROOT)
                    this->phi_plus_one_[this->get_row_index_from_matrix_indices(i, j, 0)] = calculate_phi_plus_one(i, j,
                                                                                                                   0);
                if (layer_number_ != layers_count_ - 1)
                    this->phi_plus_one_[this->get_row_index_from_matrix_indices(i, j, this->z_size_ -
                                                                                      1)] = calculate_phi_plus_one(i, j,
                                                                                                                   this->z_size_ -
                                                                                                                   1);
            }
        }
    }

    double check_result() {
        double this_delta = 0.0f;
        double current_delta;
        for (int k = 0; k < this->z_size_; ++k) {
            for (int j = 0; j < this->y_size_; ++j) {
                for (int i = 0; i < this->x_size_; ++i) {
                    current_delta = std::fabs(this->phi_[this->get_row_index_from_matrix_indices(i, j, k)] -
                                              phi(this->get_coords_from_matrix_indices(i, j, k)));
                    this_delta = std::max(this_delta, current_delta);
                }
            }
        }
        MPI_Allreduce(&this_delta, &global_delta_, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return global_delta_;
    }


    friend std::ostream &operator<<(std::ostream &os, const Grid &grid) {
        std::stringstream ss;
        for (int k = 0; k < grid.z_size_; k++) {
            os << "Layer #" << k + 1 << endl;
            for (int j = 0; j < grid.y_size_; j++) {
                for (int i = 0; i < grid.x_size_; i++) {
                    ss << std::fixed << std::setprecision(6);
                    ss << '(' << grid.get_coords_from_matrix_indices(i, j, k)[0] << ","
                       << grid.get_coords_from_matrix_indices(i, j, k)[1]
                       << ", " << grid.get_coords_from_matrix_indices(i, j, k)[2] << ") " <<
                       std::fixed << std::setprecision(6) <<
                       grid.phi_[grid.get_row_index_from_matrix_indices(i, j, k)];
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
    grid.set_borders_phi();
    double start = MPI_Wtime();
    grid.iteration();
    double end = MPI_Wtime();
    if (world_rank == 0) {
        cout << end - start << endl;
    }
    assert(grid.check_result() < 1e-4);
    MPI_Finalize();
    return 0;
}
