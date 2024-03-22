#include "MPC.h"

Vec_d floatVecToDoubleVec(const Vec_f& floatVec) {
    Vec_d doubleVec(floatVec.size());
    std::transform(floatVec.begin(), floatVec.end(), doubleVec.begin(), [](float f) { return static_cast<double>(f); });
    return doubleVec;
}

int main(){
    Vec_f wx({0.0, 60.0, 125.0,  50.0,   75.0,  35.0,  -10.0});
  // Vec_f wy({0.0,  4.0,  -4.0,  4.0,  -4.0,   4.0,  0.0});
  Vec_f wy({0.0,  0.0,  50.0,  65.0,   30.0,  50.0,  -20.0});

  cpprobotics::Spline2D csp_obj(wx, wy);
  Vec_f r_x;
  Vec_f r_y;
  Vec_f ryaw;
  Vec_f rcurvature;
  Vec_f rs;
  for(float i=0; i<csp_obj.s.back(); i+=1.0){
    std::array<float, 2> point_ = csp_obj.calc_postion(i);
    r_x.push_back(point_[0]);
    r_y.push_back(point_[1]);
    ryaw.push_back(csp_obj.calc_yaw(i));
    rcurvature.push_back(csp_obj.calc_curvature(i));
    rs.push_back(i);
  }
  // Set the size of output image to 1200x780 pixels
    plt::figure_size(1200, 780);
    // Plot line from given x and y data. Color is selected automatically.
    plt::plot(r_x, r_y);
    // Add graph title
    plt::title("Route Path");
    // Save the image (file format is determined by the extension)
    plt::save("./Route_path.png");

    int nT = 6;
    int nTu = 6;
    int nx = 4;
    double Q = 0.1;
    double R = 0.01;
    double target_speed = 10.0 / 3.6;
    Vec_d speed_profile(ryaw.size(), target_speed);
    MPC mpc = MPC(nT,nx,nTu,Q,R);


    Vec_d r_x1 = floatVecToDoubleVec(r_x);
    Vec_d r_y1 = floatVecToDoubleVec(r_y);
    Vec_d ryaw1 = floatVecToDoubleVec(ryaw);    
    Vec_d rcurvature1 = floatVecToDoubleVec(rcurvature);

    mpc.simulation(r_x1, r_y1, ryaw1, rcurvature1, speed_profile, {{wx.back(), wy.back()}});
} 