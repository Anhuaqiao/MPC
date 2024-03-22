#include <vector>
#include <Eigen/Dense>
#include <algorithm> // for std::find

typedef std::vector<double> Vec_d;
typedef std::vector<float> Vec_f;
typedef std::array<double, 2> Poi_d;

template <typename T>
T vec_diff(T x){
    T h;
    for(size_t i=0;i<x.size()-1;i++){
        h.push_back(x[i+1]-x[i]);
    }
    return h;
}

template <typename T>
T cumsum(T x){
    T h{0};
    for(size_t i=0;i<x.size();i++){
        h.push_back(h[i]+x[i]);
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
    double calc(double t){
        if(t>x.back() || t<x.front()){
            std::cout << "received value is out of th predefined range" << std::endl;
        }
        int ind = bisect(t, 0, nf);
        double dx = t - x[ind];
        return a[ind] + dx*b[ind] + dx*dx*c[ind] + dx*dx*dx*d[ind];
    }
    double calc_d(double t){
        if(t>x.back() || t<x.front()){
            std::cout << "received value is out of th predefined range" << std::endl;
        }
        int ind = bisect(t, 0, nf);
        double dx = t - x[ind];
        return b[ind] + 2*dx*c[ind] + 3*dx*dx*d[ind];
    }

    double calc_dd(double t){
        if(t>x.back() || t<x.front()){
            std::cout << "received value is out of th predefined range" << std::endl;
        }
        int ind = bisect(t, 0, nf);
        double dx = t - x[ind];
        return 2*c[ind] + 6*dx*d[ind];
    }
private:

    int bisect(double t, int start, int end){
        int mid = (end+start)/2;
        if(t == x[mid]||(end-start)<=1){
            return mid;
        }else if(t>x[mid]){
            return bisect(t, mid, end);
        }else if(t<x[mid]){
            return bisect(t, start, mid);
        }
    }
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
    // Zero_2row(0,3*nf) = 1.0;
    // Zero_2row(0,3*nf+1) = -1.0;
    // Zero_2row(1,4*nf-1) = 1.0;
    // Zero_2row(1,4*nf-2) = -1.0; //not a knot spline

    Zero_2row(0,2*nf) = 1;
    Zero_2row(0,3*nf-1) = 1;

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
    Eigen::VectorXd abcd = A.colPivHouseholderQr().solve(B);
    double * abcd_pointer = abcd.data();
    //Eigen::Map<Eigen::VectorXf>(c, c_eigen.rows(), 1) = c_eigen;
    a.assign(abcd_pointer, abcd_pointer+nf);
    b.assign(abcd_pointer+nf, abcd_pointer+2*nf);
    c.assign(abcd_pointer+2*nf, abcd_pointer+3*nf);
    d.assign(abcd_pointer+3*nf, abcd_pointer+4*nf);
    std::cout<<b.size()<<std::endl;
    for(auto const tmp:c){
      std::cout<< tmp <<std::endl;
    }

    std::cout << A << std::endl;
    std::cout<< B<< std::endl;
};

cubic_spline_::~cubic_spline_()
{

}

class cubic_spline_2D{
public:
    cubic_spline_ sx;
    cubic_spline_ sy;
    Vec_d s;

    cubic_spline_2D(Vec_d x, Vec_d y ){
        s = calc_s(x, y);
        sx = cubic_spline_(s,x);
        sy = cubic_spline_(s,y);
    }
    Poi_d calc_position(double s_t){
        double x = sx.calc(s_t);
        double y = sy.calc(s_t);

        return {{x,y}};
    }
    double calc_yaw(double s_t){
        double dx = sx.calc_d(s_t);
        double dy = sy.calc_d(s_t);
        return atan2(dy,dx);
    }
    double calc_curvature(float s_t){
    double dx = sx.calc_d(s_t);
    double ddx = sx.calc_dd(s_t);
    double dy = sy.calc_d(s_t);
    double ddy = sy.calc_dd(s_t);
    return (ddy * dx - ddx * dy)/(dx * dx + dy * dy);
  };

private:
    Vec_d calc_s(Vec_d x, Vec_d y){
        Vec_d ds;
        Vec_d out_s;
        Vec_d dx = vec_diff(x);
        Vec_d dy = vec_diff(y);

        for(size_t i = 0; i<dx.size();i++){
            ds.push_back(std::sqrt(dx[i]*dx[i] + dy[i]*dy[i]));
        }
        out_s = cumsum(ds);
        return out_s;
    };

};