#include "MPC.h"

int main(){
  // Vec_d wx({0.0, 60.0, 125.0,  50.0,   75.0,  35.0,  -10.0});
  // Vec_f wy({0.0,  4.0,  -4.0,  4.0,  -4.0,   4.0,  0.0});
  // Vec_d wy({0.0,  0.0,  50.0,  65.0,   30.0,  50.0,  -20.0});
  Vec_d wx({0.0, -10.0, -20.0, -40.0, -50.0, -60.0, -70.0});
  Vec_d wy({0.0, -1.0, 1.0, 0.0, -1.0, 1.0, 0.0});
  cubic_spline_2D csp_obj(wx, wy);
  Vec_d r_x;
  Vec_d r_y;
  Vec_d ryaw;
  Vec_d rcurvature;
  Vec_d rs;
  for(double i=0; i<csp_obj.s.back(); i+=1.0){
    std::array<double, 2> point_ = csp_obj.calc_position(i);
    r_x.push_back(point_[0]);
    r_y.push_back(point_[1]);
    ryaw.push_back(csp_obj.calc_yaw(i));
    rcurvature.push_back(csp_obj.calc_curvature(i));
    rs.push_back(i);
  }
    int nT = 6;
    int nTu = 6;
    int nx = 4;
    double Q = 0.1;
    double R = 0.01;
    double target_speed = 10.0 / 3.6;
    Vec_d speed_profile(ryaw.size(), target_speed);
    MPC mpc = MPC(nT,nx,nTu,Q,R);

    mpc.simulation(r_x, r_y, ryaw, rcurvature, speed_profile, {{wx.back(), wy.back()}});
        // save figure
    const char* filename = "./mpc_demo.png";
    std::cout << "Saving result to " << filename << std::endl;
    plt::save(filename);
    plt::show();
    return 0;
} 