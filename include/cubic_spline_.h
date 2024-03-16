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
    Eigen::VectorXd hvec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(h.data(), h.size());
    H = hvec.asDiagonal();
    Zeros.resize(nf,nf);
    Zeros.setZero();
    Identity.resize(nf,nf);
    Identity.setIdentity();
    J.resize(nf,nf);
    J.setIdentity();
    J.diagonal(1) = -Eigen::VectorXd::Ones(nf - 1);

    Eigen::MatrixXd H2 = H*H;
    Eigen::MatrixXd H3 = H2*H;
    Eigen::MatrixXd E = Identity.block(0,0,nf-1,nf);
    Eigen::MatrixXd Zero_2row;
    Zero_2row.resize(2,Identity.cols()+3*Zeros.cols());
    Zero_2row.setZero();
    Zero_2row(0,3*nf) = 1.0;
    Zero_2row(1,Identity.cols()+3*Zeros.cols()-1) = 1.0;

    Eigen::MatrixXd A1(Identity.rows(), Identity.cols()+3*Zeros.cols());
    A1<<Identity,Zeros,Zeros,Zeros;
    Eigen::VectorXd b1 = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(y.data(), y.size()-1);
    Eigen::MatrixXd A2(Identity.rows(), Identity.cols()+3*Zeros.cols());
    A2<<Identity,H,H2,H3;
    Eigen::VectorXd b2 = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(y.data()+1, y.size()-1);
    Eigen::MatrixXd A3(Identity.rows(), Identity.cols()+3*Zeros.cols());
    A3<<Zeros,J,2*H,3*H2;
    A3 = E*A3;
    Eigen::VectorXd b3;
    b3.resize(nf-1);
    b3.setZero();
    Eigen::MatrixXd A4(Identity.rows(), Identity.cols()+3*Zeros.cols());
    A4<<Zeros,Zeros,J,3*H;
    A4 = E*A4;
    Eigen::VectorXd b4;
    b4.resize(nf-1+2);
    b4.setZero();


    Eigen::MatrixXd A(A1.rows()+A2.rows()+A3.rows()+A4.rows()+Zero_2row.rows(), A1.cols());
    A<<A1,
        A2,
        A3,
        A4,
        Zero_2row;
    Eigen::VectorXd B(b1.size() + b2.size() + b3.size() + b4.size());
    B<<b1,
        b2,
        b3,
        b4;
    Eigen::VectorXd c_eigen = A.colPivHouseholderQr().solve(B);
    double * c_pointer = c_eigen.data();
    //Eigen::Map<Eigen::VectorXf>(c, c_eigen.rows(), 1) = c_eigen;
    c.assign(c_pointer, c_pointer+c_eigen.rows());

    for(auto const tmp:c){
      std::cout<< tmp <<std::endl;
    }

    std::cout << A << std::endl;
    std::cout<< B<< std::endl;
};

cubic_spline_::~cubic_spline_()
{

}
