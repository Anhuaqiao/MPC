#include "MPC.h"

int main(){
    Vec_d wx({0.0, 60.0, 125.0,  50.0,   75.0,  35.0,  -10.0});
    Vec_d wy({0.0,  0.0,  50.0,  65.0,   30.0,  50.0,  -20.0});

    cubic_spline_ c(wx, wy);
}