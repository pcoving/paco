#include "paco.hpp"

template <int NQUAD>
class MyPaco : public Paco<12> {
public:

  MyPaco() {
    cout << "MyPaco()" << endl;

    assert(NQUAD == 12);
  }
  
  void initialHook() {
    
        
    FOR_IP {
      mp[ip] = 0.0;
      FOR_I2 gp[ip][i] = 0.0;    
    }

    FOR_IT {
      triVec[it].updateVeX(xp);
      triVec[it].calcArea();
      
      FOR_IQ {
	double xquad[2] = {0.0, 0.0};
	FOR_I2 FOR_IVE xquad[i] += triVec[it].ve_x[ive][i]*quad.ve_wgt[iq][ive];
	
	int icell;
	int jcell;
	
	if (xquad[0] < 0.0 && xquad[1] >= 0.0) {
	  icell = 1;
	  jcell = 1;
	}
	else if (xquad[0] >= 0.0 && xquad[1] >= 0.0) {
	  icell = 2;
	  jcell = 1;
	}
	else if (xquad[0] < 0.0 && xquad[1] < 0.0) {
	  icell = 1;
	  jcell = 2;
	  }
	else {
	  assert(xquad[0] >= 0.0 && xquad[1] < 0.0);
	  icell = 2;
	  jcell = 2;
	}
	
	double this_rho, this_p;
	double this_u[2];
	calcEulerVortex(this_rho, this_u, this_p, xquad);
			
	const double lambda[2] = {triVec[it].lambda[iq][0],triVec[it].lambda[iq][1]};
	const double weight = quad.p_wgt[iq]*triVec[it].area/triVec[it].denom[iq];
	
	//#pragma omp parallel for 
	for (int ivar = -1; ivar <= 1; ++ivar) {
	  for (int jvar = -1; jvar <= 1; ++jvar) {
	    const int this_icell = icell + ivar;
	    const int this_jcell = jcell + jvar;
	    for (typename list<Point>::iterator il = cells[this_icell + ncells*this_jcell].begin(); il != cells[this_icell + ncells*this_jcell].end(); ++il) {
	      const double val = evalExp(xquad, il->x, lambda);
	      mp[il->ip] += this_rho*weight*val;
	      FOR_I2 up[il->ip][i] += this_u[i]*weight*val;
	    }
	  }
	}
	
      }
    }
    
  }
  
  void temporalHook() {
    
    /*
    for (int i = 0; i < ncells; ++i) {
      for (int j = 0; j < ncells; ++j) {
	for (typename list<Point>::iterator il = cells[i + ncells*j].begin(); il != cells[i + ncells*j].end(); ++il) {
	  il->tmp = 0.0;
        }
      }
    }
    
    
    FOR_IT {
      triVec[it].updateVeX(xp);
      triVec[it].calcArea();
      
      FOR_IQ {
	double xquad[2] = {0.0, 0.0};
	FOR_I2 FOR_IVE xquad[i] += triVec[it].ve_x[ive][i]*quad.ve_wgt[iq][ive];
	
	int icell, jcell;
	
	if (xquad[0] < 0.0 && xquad[1] >= 0.0) {
	  icell = 1;
	  jcell = 1;
	}
	else if (xquad[0] >= 0.0 && xquad[1] >= 0.0) {
	  icell = 2;
	  jcell = 1;
	}
	else if (xquad[0] < 0.0 && xquad[1] < 0.0) {
	  icell = 1;
	  jcell = 2;
	  }
	else {
	  assert(xquad[0] >= 0.0 && xquad[1] < 0.0);
	  icell = 2;
	  jcell = 2;
	}
	
	double this_rho, this_p;
	double this_u[2];
	calcEulerVortex(this_rho, this_u, this_p, xquad);
			
	const double lambda[2] = {triVec[it].lambda[iq][0],triVec[it].lambda[iq][1]};
	const double weight = quad.p_wgt[iq]*triVec[it].area/triVec[it].denom[iq];
	
        #pragma omp parallel for collapse(2) 
	for (int ivar = -1; ivar <= 1; ++ivar) {
	  for (int jvar = -1; jvar <= 1; ++jvar) {
	    const int this_icell = icell + ivar;
	    const int this_jcell = jcell + jvar;
	    for (typename list<Point>::iterator il = cells[this_icell + ncells*this_jcell].begin(); il != cells[this_icell + ncells*this_jcell].end(); ++il) {
	      const double val = evalExp(xquad, il->x, lambda);
	      il->tmp += this_rho*weight*val;
	    }
	  }
	}
	
      }
    }

   

    double * mp_exact = new double[np];
    FOR_IP mp_exact[ip] = 0.0;
    
    for (int i = 0; i < ncells; ++i) {
      for (int j = 0; j < ncells; ++j) {
        for (typename list<Point>::iterator il = cells[i + ncells*j].begin(); il != cells[i + ncells*j].end(); ++il) {
          mp_exact[il->ip] += il->tmp;
        }
      }
    }
    

    double vol_sum = 0.0;
    double l2_error = 0.0;
    FOR_IP {
      vol_sum += volp[ip];
      const double delta = (mp_exact[ip] - mp[ip])/volp[ip];
      // rho_error[ip] = delta;
      l2_error += delta*delta*volp[ip];
    }

    delete[] mp_exact;
    */

    #pragma omp parallel for 
    FOR_IP {
      int icell, jcell;
      
      if (xp[ip][0] < 0.0 && xp[ip][1] >= 0.0) {
	icell = 1;
	jcell = 1;
      }
      else if (xp[ip][0] >= 0.0 && xp[ip][1] >= 0.0) {
	icell = 2;
	jcell = 1;
      }
      else if (xp[ip][0] < 0.0 && xp[ip][1] < 0.0) {
	icell = 1;
	jcell = 2;
      }
      else {
	assert(xp[ip][0] >= 0.0 && xp[ip][1] < 0.0);
	icell = 2;
	jcell = 2;
      }

      double this_rho, this_p;
      double this_u[2];
      calcEulerVortex(this_rho, this_u, this_p, xp[ip]);
      rho_error[ip] = this_rho;

      double coeff = 1.0/denomp[ip];

      for (int ivar = -1; ivar <= 1; ++ivar) {
	for (int jvar = -1; jvar <= 1; ++jvar) {
	  const int this_icell = icell + ivar;
	  const int this_jcell = jcell + jvar;
	  for (typename list<Point>::iterator il = cells[this_icell + ncells*this_jcell].begin(); il != cells[this_icell + ncells*this_jcell].end(); ++il) {
	    rho_error[ip] -= coeff*rhop[il->ip]*evalExp(xp[ip], il->x, lambdap[ip]);
	  }
	}
      }
      
      rho_error[ip] *= volp[ip];
    }

    
    double vol_sum = 0.0;
    double l2_error = 0.0;
    FOR_IP {
      vol_sum += volp[ip];
      l2_error += rho_error[ip]*rho_error[ip]*volp[ip];
    }

    cout << " > time, L2: " << time << " " << sqrt(l2_error/vol_sum) << endl;

    
  }

  void finalHook() {
    
  }
  
  void updateU(double (*_u)[2]) {
    #pragma omp parallel for 
    FOR_IP {
      double rho, p;
      double u[2];
      calcEulerVortex(rho, u, p, xp[ip]);
      
      FOR_I2 _u[ip][i] = u[i];
    }
  }

  inline void calcEulerVortex(double &rho, double u[2], double &p, const double x[2]) {
    
    const double u_inf = 0.5;
    
    // x,y position of the vortex
    const double x0 = 0.0;
    const double y0 = 0.0;
  
    const double Ma_inf = u_inf/sqrt(gam);
    const double rc = 1.0;
  
    // circulation parameter...
    const double e_twopi = 0.08; // normal
    //const double e_twopi = 0.16;
    
    // setup...
    const double coeff = 0.5 * e_twopi*e_twopi * (gam-1.0) * Ma_inf*Ma_inf;
    
    double dx = x[0] - x0;
    double dy = x[1] - y0;
    
    const double f0 = 1.0 - (( dx*dx ) + ( dy*dy ))/( rc*rc );
    rho = pow( 1.0 - coeff * exp( f0 ) , 1.0/(gam-1.0) );
      
    u[0] = u_inf*(-e_twopi * ( dy )/rc * exp( f0 / 2.0 ));
    u[1] = u_inf*( e_twopi * ( dx )/rc * exp( f0 / 2.0 ));
    
    p = pow( 1.0 - coeff * exp( f0 ) , gam/(gam-1.0) );
        
  }

};

int main(int argc,char * argv[]) {

  MyPaco<12> * solver = new MyPaco<12>;
  solver->init();
  solver->run();

  return(0);
}
