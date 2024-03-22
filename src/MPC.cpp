#include "MPC.h"

int main(){
    Vec_d wx({0.0, 60.0, 125.0,  50.0,   75.0,  35.0,  -10.0});
    Vec_d wy({0.0,  0.0,  50.0,  65.0,   30.0,  50.0,  -20.0});

    // Vec_d wx({35.0, 10.0, 0.0, 0.0,0.0, 30.0, 6.0, 20.0, 35.0});
    // Vec_d wy({0.0, 0.0, 20.0, 35.0, 20.0, 20.0, 30.0, 5.0, 0.0});

    cubic_spline_2D c(wx, wy);
    Vec_d r_x;
    Vec_d r_y;
    Vec_d ryaw;
    Vec_d rcurvature;
    for(double i = 0;i<c.s.back();i+=1.0){
        std::array<double,2>point = c.calc_position(i);
        r_x.push_back(point[0]);
        r_y.push_back(point[1]);
        ryaw.push_back(c.calc_yaw(i));
        rcurvature.push_back(c.calc_curvature(i));
        std::cout<<point[0]<<" "<<point[1]<<std::endl;
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
    double Q = 0.01;
    double R = 1;
    double target_speed = 10.0 / 3.6;
    Vec_d speed_profile(ryaw.size(), target_speed);
    MPC mpc = MPC(nT,nx,nTu,Q,R);
    mpc.simulation(r_x, r_y, ryaw, rcurvature, speed_profile, {{wx.back(), wy.back()}});
} 