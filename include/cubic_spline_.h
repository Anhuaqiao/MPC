#include <vector>
#include <Eigen/Dense>
#include <algorithm> // for std::find

typedef std::vector<double> Vec_d;
typedef std::vector<float> Vec_f;

template <typename T>
T vec_diff(T x){
    T h;
    for(size_t i=0;i<x.size()-1;i++){
        h.push_back(x[i+1]-x[i]);
    }
    return h;
}

class cubic_spline_
{
private:
    Vec_d x;
    Vec_d y;
    Vec_d h;
    Vec_d a,b,c,d;
    Eigen::MatrixXd H;
    Eigen::MatrixXd Identity;
    Eigen::MatrixXd Zeros;
    Eigen::MatrixXd J;
    size_t nx,nf;

public:
    cubic_spline_(){};
    cubic_spline_(Vec_d x_, Vec_d y_);
    ~cubic_spline_();
};

cubic_spline_::cubic_spline_(Vec_d x_, Vec_d y_):x(x_), y(y_), nx(x_.size()),h(vec_diff(x_)){
    nf = nx-1;
    double * h_pointer = h.data();
    Eigen::Map<Eigen::VectorXd> h_vec(h_pointer, 4);
    H = h_vec.asDiagonal();
    Zeros.resize(nf,nf);
    Zeros.setZero();
    Identity.resize(nf,nf);
    Identity.setIdentity();
    J.resize(nf,nf);
    J.setIdentity();
    J.diagonal(1) = -Eigen::VectorXd::Ones(nf - 1);

    Eigen::MatrixXd A1(Identity.rows(), Identity.cols()+3*Zeros.cols());
    Eigen::MatrixXd H2 = H*H;
    Eigen::MatrixXd H3 = H2*H;
    Eigen::MatrixXd E = Identity.block(0,0,nf-1,nf);
    std::cout<< E << std::endl;

};

cubic_spline_::~cubic_spline_()
{

}
