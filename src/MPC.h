
#ifndef MPC_H
#define MPC_H

#include <iostream>
#include "cubic_spline_.h"
#include "matplotlibcpp.h"
#include <vector>
#include <Eigen/Dense>
#include <algorithm> // for std::find
#include<cppad/cppad.hpp>
#include<cppad/ipopt/solve.hpp>

#define MAX_STEER 45.0/180*M_PI
#define MAX_DSTEER  30.0/180*M_PI
#define MAX_SPEED   55.0/3.6
#define MIN_SPEED  -20.0/3.6
#define MAX_ACCEL 1.0

#define WB 2.5
#define LENGTH  4.5
#define WIDTH 2.0
#define BACKTOWHEEL 1.0
#define WHEEL_LEN 0.3
#define WHEEL_WIDTH 0.2
#define TREAD 0.7

#define DT 0.2
#define PI 3.1415926

#define N_IND_SEARCH 20

namespace plt = matplotlibcpp;

namespace mpctxh{
struct State{
  double x;
  double y;
  double yaw;
  double v;
  State(double x_, double y_, double yaw_, double v_){
    x = x_;
    y = y_;
    yaw = yaw_;
    v = v_;
  };
};
};

using namespace mpctxh;

template <typename T>
T PI2PI(T angle){
    if(angle>PI) angle -= 2*PI;
    if(angle<-PI) angle += 2*PI;
    return angle;
};

class FG_EVAL{
public:
    Eigen::MatrixXd ref;
    int nT;
    int nx;
    int nTu;
    int xstart;
    int ystart;
    int yawstart;
    int vstart;
    int delta_start;
    int a_start;
    double Q,R;
    FG_EVAL(Eigen::MatrixXd ref_, int nT_, int nx_, int nTu_,double Q_, double R_, State x0_):nT(nT_),nx(nx_),nTu(nTu_),Q(Q_),R(R_){
    this->xstart = 0;
    this->ystart = this->xstart+nT+1;
    this->yawstart = this->ystart+nT+1;
    this->vstart = this->yawstart+nT+1;
    this->delta_start = this->vstart+nT+1;
    this->a_start = this->delta_start+nTu;
    ref = ref_;
};
    typedef CPPAD_TESTVECTOR(CppAD::AD<double>) ADvector;
    void operator()(ADvector& fg, ADvector& vars){
        fg[0]=0;
        for(int i = 0;i<nTu;i++){
            fg[0] += Q*CppAD::pow(vars[delta_start+i],2);
            fg[0] += Q*CppAD::pow(vars[a_start],2);
        }
        //for(int i = 0;i<nTu-1;i++){
        //    fg[0] += Q*CppAD::pow(vars[delta_start+i+1]-vars[delta_start+i],2);
        //    fg[0] += Q*CppAD::pow(vars[a_start+i+1]-vars[a_start+i],2);
        //}
        fg[1+xstart] = vars[xstart];
        fg[1+ystart] = vars[ystart];
        fg[1+yawstart] = vars[yawstart];
        fg[1+vstart] = vars[vstart];

        for(int i=0;i<nT;i++){
            CppAD::AD<double> x1 = vars[xstart+i+1];
            CppAD::AD<double> y1 = vars[ystart+i+1];
            CppAD::AD<double> yaw1 = vars[yawstart+i+1];
            CppAD::AD<double> v1 = vars[vstart+i+1];

            CppAD::AD<double> x0 = vars[xstart+i];
            CppAD::AD<double> y0 = vars[ystart+i];
            CppAD::AD<double> yaw0 = vars[yawstart+i];
            CppAD::AD<double> v0 = vars[vstart+i];
            CppAD::AD<double> delta;
            CppAD::AD<double> a;
            if(i>=nTu){
                delta = vars[delta_start+nTu-1];
                a = vars[a_start+nTu-1];
            }else{
                delta = vars[delta_start+i];
                a = vars[a_start+i];
            }


            fg[2+xstart+i] = x1 - x0 + v0*cos(yaw0)*DT;
            fg[2+ystart+i] = y1 - y0 + v0*sin(yaw0)*DT;
            fg[2+yawstart+i] = yaw1 - yaw0 + v0/WB*tan(delta)*DT;
            fg[2+vstart+i] = v1 - v0+a*DT;

            fg[0] += CppAD::pow(ref(0,i)-x1, 2);
            fg[0] += CppAD::pow(ref(1,i)-y1, 2);
            fg[0] += CppAD::pow(ref(2,i)-yaw1, 2);
            fg[0] += CppAD::pow(ref(3,i)-v1, 2);
        }

        int tmp = 0;
        for(const auto a:fg){
            std::cout<< "fg " << tmp << ": " << a << std::endl;
            tmp++;
        }
    };
};


class MPC{
public:
    int nT;
    int nx;
    int nTu;
    int xstart;
    int ystart;
    int yawstart;
    int vstart;
    int delta_start;
    int a_start;
    double Q,R;
    int target_int;
    MPC(void);
    MPC(int nT_, int nx_, int nu_,double Q_, double R_);
    void calc_ref(Vec_d cx, Vec_d cy, Vec_d cyaw, Vec_d ck, Vec_d speed_profile, Eigen::MatrixXd& ref);
    Vec_d solve(State x0, double delta, double a, Eigen::MatrixXd RefPath);
    void simulation(Vec_d cx, Vec_d cy, Vec_d cyaw, Vec_d ck, Vec_d speed_profile, Poi_d goal);
    void update(State& state, float a, float delta);
    int calc_target_index(State state, Vec_d cx, Vec_d cy, int pind);
    void calc_ref(State state, Vec_d cx, Vec_d cy, Vec_d cyaw, Vec_d ck, Vec_d speed_profile, int target_int_, Eigen::MatrixXd& ref);
};

MPC::MPC(void){};
MPC::MPC(int nT_, int nx_, int nTu_,double Q_, double R_):nT(nT_),nx(nx_),nTu(nTu_),Q(Q_),R(R_){
    this->xstart = 0;
    this->ystart = this->xstart+nT+1;
    this->yawstart = this->ystart+nT+1;
    this->vstart = this->yawstart+nT+1;
    this->delta_start = this->vstart+nT+1;
    this->a_start = this->delta_start+nTu;

};
Vec_d MPC::solve(State x0, double delta, double a, Eigen::MatrixXd RefPath){
    typedef CPPAD_TESTVECTOR(double) Dvector;
    int nvar = 4*(nT+1)+2*nTu;
    int n_constraints = 4*(nT+1);
    Dvector var(nvar);
    for(int i=0;i<nvar;i++){
        var[i]=0;
    }

    //initial state
    var[xstart] = x0.x;
    var[ystart] = x0.y;
    var[yawstart] = x0.yaw;
    var[vstart] = x0.v;
    var[delta_start] = delta;
    var[a_start] = a;

    Dvector varl(nvar);
    Dvector varu(nvar);
    Dvector gl(n_constraints);
    Dvector gu(n_constraints);
    
    // variable constraints
    for(int i=xstart;i<vstart;i++){
        varl[i] = -10000000;
        varu[i] = 10000000;
    }
    for(int i=vstart;i<delta_start;i++){
        varl[i] = -10000000;
        varu[i] = 10000000;
    }
    for(int i=delta_start;i<a_start;i++){
        varl[i] = -MAX_STEER;
        varu[i] = MAX_STEER;
    }
    for(int i=a_start;i<nvar;i++){
        varl[i] = -MAX_ACCEL;
        varu[i] = MAX_ACCEL;
    }

    for(int i=0;i<n_constraints;i++){
        gl[i] = 0;
        gu[i] = 0;
    }
    gl[xstart] = x0.x;
    gl[ystart] = x0.y;
    gl[yawstart] = x0.yaw;
    gl[vstart] = x0.v;

    gu[xstart] = x0.x;
    gu[ystart] = x0.y;
    gu[yawstart] = x0.yaw;
    gu[vstart] = x0.v;
    int tmp = 0;
        for(const auto a:var){
        std::cout<< "var " << tmp << ": " << a << std::endl;
        tmp++;
    }
    tmp=0;
    for(const auto a:varl){
        std::cout<< "varl " << tmp << ": " << a << std::endl;
        tmp++;
    }
    tmp=0;
    for(const auto a:varu){
        std::cout<< "varu"  << tmp << ": "  << a << std::endl;
        tmp++;
    }
    tmp=0;
    for(const auto a:gl){
        std::cout<< "gl"  << tmp << ": "  << a << std::endl;
        tmp++;
    }
    tmp=0;
    for(const auto a:gu){
        std::cout<< "gu"  << tmp << ": " << a << std::endl;
        tmp++;
    }


    FG_EVAL fg_eval(RefPath, nT, nx, nTu, Q, R, x0);

    // options
    std::string options;
    options += "Integer print_level  0\n";
    // options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";
    options += "Integer max_iter      50\n";

    CppAD::ipopt::solve_result<Dvector> solution;

    CppAD::ipopt::solve<Dvector,FG_EVAL>(
        options, var, varl, varu, gl, gu, fg_eval, solution);
    
    bool ok = true;
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    Vec_d result;
    for(auto i=0;i<nvar;i++){
        result.push_back((double)solution.x[i]);
    }
    return result;
};

int MPC::calc_target_index(State state, Vec_d cx, Vec_d cy, int pind){
    double mind = std::numeric_limits<double>::max();
    double ind = 0;
    for(int i=pind;i<pind+N_IND_SEARCH;i++){
        double idx = cx[pind] - state.x;
        double idy = cy[pind] - state.y;
        double d_e = idx*idx + idy*idy;

        if(d_e<mind){
            mind = d_e;
            ind = i;
        }
    }
    return ind;

};

void MPC::calc_ref(State state, Vec_d cx, Vec_d cy, Vec_d cyaw, Vec_d ck, Vec_d speed_profile, int target_int_, Eigen::MatrixXd& ref){
    ref.resize(4,nT);
    ref.setZero();

    int ind = calc_target_index(state, cx, cy, target_int_);
    target_int = ind;
    int ncourse = cx.size();

    float travel = 0.0;

    for(int i=0;i<nT;i++){
        travel += std::abs(state.v)*DT;
        int dind = (int)std::round(travel/1);

        if((ind+dind)<ncourse){
            ref(0,i) = cx[ind+dind];
            ref(1,i) = cy[ind+dind];
            ref(2,i) = cyaw[ind+dind];
            ref(3,i) = speed_profile[ind+dind];
        }else{
            ref(0,i) = cx[ncourse-1];
            ref(1,i) = cy[ncourse-1];
            ref(2,i) = cyaw[ncourse-1];
            ref(3,i) = speed_profile[ncourse-1];
        }
    }
};

void MPC::update(State& state, float a, float delta){
  if (delta >= MAX_STEER) delta = MAX_STEER;
  if (delta <= - MAX_STEER) delta = - MAX_STEER;

  state.x = state.x + state.v * std::cos(state.yaw) * DT;
  state.y = state.y + state.v * std::sin(state.yaw) * DT;
  state.yaw = state.yaw + state.v / WB * CppAD::tan(delta) * DT;
  state.v = state.v + a * DT;

  if (state.v > MAX_SPEED) state.v = MAX_SPEED;
  if (state.v < MIN_SPEED) state.v = MIN_SPEED;

};

void MPC::simulation(Vec_d cx, Vec_d cy, Vec_d cyaw, Vec_d ck, Vec_d speed_profile, Poi_d goal){
    State state(cx[0],cy[0],cyaw[0],speed_profile[0]);
    double goal_dist = 0.5;
    int iter_count=0;
    target_int=0;
    Vec_d x_h;
    Vec_d y_h;

    Eigen::MatrixXd ref;
    double delta = 0;
    double a = 0;
    while (5000>=iter_count)
    {
        calc_ref(state, cx, cy, cyaw, ck, speed_profile, target_int, ref);
        std::cout << "ref matrix: " << ref<< std::endl;
        Vec_d output = solve(state, delta, a, ref);
        a = output[a_start];
        delta = output[delta_start];
        update(state, output[a_start], output[delta_start]);
        float steer = output[delta_start];
        float dx = state.x - goal[0];
        float dy = state.y - goal[1];
        if (std::sqrt(dx*dx + dy*dy) <= goal_dist) {
        std::cout<<("Goal")<<std::endl;
        break;
        }
        x_h.push_back(state.x);
        y_h.push_back(state.y);
        iter_count++;
    }
  // Set the size of output image to 1200x780 pixels
    plt::figure_size(1200, 780);
    // Plot line from given x and y data. Color is selected automatically.
    plt::plot(x_h, y_h);
    // Add graph title
    plt::title("Route Path");
    // Save the image (file format is determined by the extension)
    plt::save("./vehicle-path.png");



};
#endif //MPC_H