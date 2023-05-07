#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip>


constexpr double Dx = 2.0f;
constexpr double Dy = 2.0f;
constexpr double Dz = 2.0f;
constexpr double N_x = 4;
constexpr double N_y = 4;
constexpr double N_z = 4;
constexpr double a = 1e5;
constexpr double epsilon = 1e-8;
constexpr double h_x = Dx / (N_x - 1);
constexpr double h_y = Dx / (N_y - 1);
constexpr double h_z = Dx / (N_z - 1);
constexpr double alpha = 1.0f / ((2.0f / h_x * h_x) + (2.0f / h_y * h_y) + (2.0f / h_z * h_z) + a);


class Grid {
private:
    std::vector<double> phi_;
    std::vector<double> phi_plus_one_;
    int volume_;
    int x_size_;
    int y_size_;
    int z_size_;
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

    double phi(double x, double y, double z) {
        return (x * x + y * y + z * z);
    }

    double rho(double x, double y, double z) {
        return 6 - a * phi(x, y, z);
    }


    int get_index_from_coords(int x, int y, int z) {
        return x * this->y_size_ * this->z_size_ + z * this->y_size_ + y;
    }

    void set_borders_phi() {
        for (int i = 0; i < this->z_size_; i++) {
            for (int j = 0; j < this->y_size_; j++) {
                for (int k = 0; k < this->x_size_; k++) {
                    if (i == 0 || i == this->z_size_-1 || j == 0 || j == this->y_size_-1 || k == 0 || k ==this->x_size_-1) {
                        if ((i == 0 || i == this->z_size_-1) && (j == 0 || j == this->y_size_-1) ||
                            (i == 0 || i == this->z_size_-1) && (k == 0 || k == this->x_size_-1) ||
                            (j == 0 || j == this->z_size_-1) && (k == 0 || k == this->x_size_-1)) {
                           this->phi_[this->get_index_from_coords(i,j,k)] = this->phi(i,j,k);
                        }
                    }
                }
            }
        }

    }

    friend std::ostream &operator<<(std::ostream &os, Grid obj) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(12);
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
};



int main() {
    Grid grid(N_x, N_y, N_z);
    grid.set_borders_phi();
    std::cout << grid;
    return 0;
}

