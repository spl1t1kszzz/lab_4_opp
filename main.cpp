#include <iostream>
#include <vector>
#include <sstream>
#include <limits>
#include <iomanip>

constexpr double Dx = 2.0f;
constexpr double Dy = 2.0f;
constexpr double Dz = 2.0f;
constexpr double N_x = 6;
constexpr double N_y = 6;
constexpr double N_z = 6;
constexpr double a = 1e5;
constexpr double epsilon = 1e-8;
constexpr double h_x = Dx / (N_x - 1);
constexpr double h_y = Dy / (N_y - 1);
constexpr double h_z = Dz / (N_z - 1);
constexpr double alpha = 1.0f / ((2.0f / h_x * h_x) + (2.0f / h_y * h_y) + (2.0f / h_z * h_z) + a);


class Grid {
private:
    std::vector<double> phi_;
    std::vector<double> phi_plus_one_;
    std::vector<double> &current_phi_ = phi_;
    int volume_;
    int x_size_;
    int y_size_;
    int z_size_;
    double delta_ = 0.0f;
public:
    Grid(int x_size, int y_size, int z_size) {
        this->x_size_ = x_size;
        this->y_size_ = y_size;
        this->z_size_ = z_size;
        this->volume_ = this->x_size_ * this->y_size_ * this->z_size_;
        std::fill(this->phi_.begin(), this->phi_.end(), 0.0f);
        this->phi_.resize(this->volume_);
        this->phi_plus_one_.resize(this->volume_);
    }


    static double phi(double x, double y, double z) {
        return (x * x + y * y + z * z);
    }

    double rho(double x, double y, double z) {
        return 6 - a * phi(x, y, z);
    }


    int get_index_from_coords(int x, int y, int z) const {
        return x * this->y_size_ * this->z_size_ + z * this->y_size_ + y;
    }

    void set_borders_phi() {
        for (int i = 0; i < this->z_size_; i++) {
            for (int j = 0; j < this->y_size_; j++) {
                for (int k = 0; k < this->x_size_; k++) {
                    if (i == 0 || i == this->z_size_ - 1 || j == 0 || j == this->y_size_ - 1 || k == 0 ||
                        k == this->x_size_ - 1) {
                        if ((i == 0 || i == this->z_size_ - 1) && (j == 0 || j == this->y_size_ - 1) ||
                            (i == 0 || i == this->z_size_ - 1) && (k == 0 || k == this->x_size_ - 1) ||
                            (j == 0 || j == this->z_size_ - 1) && (k == 0 || k == this->x_size_ - 1)) {
                            this->current_phi_[this->get_index_from_coords(i, j, k)] = this->phi(i, j, k);
                        }
                    }
                }
            }
        }
    }

    double calculate_phi_plus_one(int x, int y, int z) {
        double tmp = (this->get_current_phi(x + 1, y, z) + this->get_current_phi(x - 1, y, z)) / (h_x * h_x) +
                     (this->get_current_phi(x, y + 1, z) + this->get_current_phi(x, y - 1, z)) / (h_y * h_y) +
                     (this->get_current_phi(x, y, z + 1) + this->get_current_phi(x, y, z - 1)) / (h_z * h_z) -
                     this->rho(x, y, z);
        return alpha * tmp;
    }

    void iteration() {
        this->delta_ = 0.0f;
        for (int i = 0; i < this->z_size_; ++i) {
            for (int j = 0; j < this->y_size_; ++j) {
                for (int k = 0; k < this->x_size_; ++k) {
                    this->phi_plus_one_[this->get_index_from_coords(k, j, i)] = this->calculate_phi_plus_one(k, j, i);
                    double current_delta = this->phi_plus_one_[this->get_index_from_coords(k, j, i)] -
                                           this->phi_[this->get_index_from_coords(k, j, i)];
                    this->delta_ = std::max(this->delta_, current_delta);
                }
            }
        }
        std::swap(this->phi_, this->phi_plus_one_);
    }

    double get_current_phi(int i, int j, int k) {
        if (i < 0 || j < 0 || k < 0 || i >= this->x_size_ || j >= this->y_size_ || k >= this->z_size_) {
            return this->phi(i, j, k);
        }
        return this->phi_[this->get_index_from_coords(i, j, k)];
    }

    friend std::ostream &operator<<(std::ostream &os, Grid obj) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(5);
        for (int i = 0; i < obj.z_size_; i++) {
            for (int j = 0; j < obj.y_size_; j++) {
                for (int k = 0; k < obj.x_size_; k++) {
                    ss << '(' << i << j << k << ") " << obj.phi_[obj.get_index_from_coords(i, j, k)];
                    os << ss.str() << ' ';
                    ss.str("");
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }
        return os;
    }

    void run() {
        double delta = std::numeric_limits<double>::infinity();
        while (delta >= epsilon) {
            this->iteration();
            delta = this->delta_;
            std::cout << delta << std::endl;
        }
    }
};


int main() {
    Grid grid(Dx / h_x, Dy / h_y, Dz / h_z);
    grid.set_borders_phi();
    grid.run();
    std::cout << grid;
    return 0;
}

