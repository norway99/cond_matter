#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <chrono>
#include <random>

#define PI 3.1415926535879
#define MAX_STEPS 1.e6
#define ETA 0.89 // viscosity for water at 398 K
#define DIFF 0.01
#define PRINTERVAL 0.05

int main(int argc, char* argv[]){

  if (argc != 5) {
    std::cerr << "Usage: ./a.out particle_mass particle_radius alpha beta" << std::endl;
    std::cerr << "For starters try: ./a.out 0.00001 0.000001 0.0003 0.0005" << std::endl;
    return 1;
  }

  double mass = atof(argv[1]);
  double tau = mass/(6*PI*ETA*atof(argv[2]));
  double a = atof(argv[3])/mass;
  double b = atof(argv[4])/mass;
		    
  FILE *plotmsd, *plotrms;

  char command[15];
  sprintf(command, "mkdir ./traj");
  std::cout << "Making directory ./traj" << std::endl; 
  system(command);
  
  plotrms = fopen("./rms.txt","w");
  fprintf(plotrms, "step = \trms = \n");
  plotmsd = fopen("./msd.txt", "w");
  fprintf(plotmsd, "step = \tmsd = \n");
  
  double h[12]; // array of step sizes
  double coeff[4] = {1., 2.5, 5.0, 7.5};
  double len = sizeof(h)/sizeof(*h);
  int k = 0;
  
  for (int i = 3; i > 0; i--) 
    for (int j = 0; j < 4; j++){
      h[k] = coeff[j]*pow(10, -i); // populate h with log-scaled step sizes
      k++;
    }
  
  for (int i = 0; i < len; i++) {
    
    double x = 0., xdot = 0., t = 0., tprint = 0., sqx = 0., rms = 0.;
    double amp = sqrt(2*DIFF/h[i])/tau; 
    int numsteps = 0;

    char trj[30];
    sprintf(trj, "./traj/traj_%d.txt", i); 
    FILE *traj = fopen(trj, "w");
    fprintf(traj, "t = \tx = \txdot = \n");

    std::cout << "Calculating trajectory for step size " << h[i] << " & writing to file " << trj << std::endl; 
    
    while ((numsteps++) < MAX_STEPS) {
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::default_random_engine generator(seed);
      std::normal_distribution<double> dist(0., 1.);
      double noise = amp*dist(generator); // generate gaussian noise with rms amplitude given by einstein-stokes
      x += h[i]*xdot;                    
      sqx += x*x;
      xdot += h[i]*(noise - tau*xdot + 2*a*x - 4*b*pow(x, 3)); // langevin eqn. w/ U(x) = -alpha*x^2 + beta*x^4 
      rms += xdot*xdot;                                        // with coefficients scaled to a, b for particle mass
      tprint += h[i];
      if (tprint >= PRINTERVAL){ // only dump positions and velocities every t = 0.5 
	fprintf(traj, "%f   \t%f \t%f \n", t, x, xdot); 
	tprint = 0;
      }
      t += h[i];
    }

    sqx /= numsteps; // mean square position of particle over the course of the trajectory
    rms = sqrt(rms)/numsteps; // rms velocity of particle 
    std::cout << "Calculating mean square displacement and rms velocity" << std::endl; 
    fprintf(plotmsd, "%f \t%f\n", h[i], sqx);
    fprintf(plotrms, "%f \t%f\n", h[i], rms);
  }

  return 0;

}

    
	
 
    
  
  

  

