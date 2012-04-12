#include <iostream>
#include <vector>
#include <list>
#include <math.h>
#include <assert.h>
#include <stdio.h>

using namespace std;

#include "quad.hpp"

#define XMIN -5.0
#define XMAX  5.0
#define YMIN -5.0
#define YMAX  5.0

#define XM_BIT 1
#define XP_BIT 2
#define YM_BIT 4
#define YP_BIT 8

#define FOR_IP for (int ip = 0; ip < np; ++ip)
#define FOR_IT for (int it = 0; it < triVec.size(); ++it)

#define FOR_I2 for (int i = 0; i < 2; ++i)
#define FOR_I3 for (int i = 0; i < 3; ++i)
#define FOR_J2 for (int j = 0; j < 2; ++j)
#define FOR_J3 for (int j = 0; j < 3; ++j)
#define FOR_IVE for (int ive = 0; ive < 3; ++ive) 
#define FOR_IQ for (int iq = 0; iq < NQUAD; ++iq) 

template <int NQUAD >
class Paco { 
public:  
  class Tri {
  public: 
    int ve_ip[3];
    int ve_bit[3];
    int ve_itnbr[3];
    int ve_ivenbr[3];
    double ve_x[3][2];
    double lambda[NQUAD][2];
    double denom[NQUAD];
    double g[NQUAD][2];
    double H[NQUAD][2][2];

    double area;

    Tri() {
      FOR_IQ FOR_I2 this->lambda[iq][i] = 0.0;
    }

    
    void calcArea() {
      area = 0.5*((ve_x[1][0]-ve_x[0][0])*(ve_x[2][1]-ve_x[0][1])-(ve_x[1][1]-ve_x[0][1])*(ve_x[2][0]-ve_x[0][0]));
      assert(area > 0.0);
    }
    
    void updateVeX(double (*xp)[2]) {
      FOR_IVE {
	if (ve_bit[ive] & XM_BIT) ve_x[ive][0] = xp[ve_ip[ive]][0] - (XMAX - XMIN);
	else if (ve_bit[ive] & XP_BIT) ve_x[ive][0] = xp[ve_ip[ive]][0] + (XMAX - XMIN);
	else ve_x[ive][0] = xp[ve_ip[ive]][0];
	
	if (ve_bit[ive] & YM_BIT) ve_x[ive][1] = xp[ve_ip[ive]][1] - (YMAX - YMIN);
	else if (ve_bit[ive] & YP_BIT) ve_x[ive][1] = xp[ve_ip[ive]][1] + (YMAX - YMIN);
	else ve_x[ive][1] = xp[ve_ip[ive]][1];
      }
    }
    
    void toggleXMBit(const int ive) {
      if (ve_bit[ive] & XP_BIT) 
	ve_bit[ive] &= ~XP_BIT;
      else {
	assert( ~(ve_bit[ive] & XM_BIT) );
	ve_bit[ive] |= XM_BIT;
      }
    }
    void toggleXPBit(const int ive) {
      if (ve_bit[ive] & XM_BIT) 
	ve_bit[ive] &= ~XM_BIT;
      else {
	assert( ~(ve_bit[ive] & XP_BIT) );
	ve_bit[ive] |= XP_BIT;
      }
    }
    void toggleYMBit(const int ive) {
      if (ve_bit[ive] & YP_BIT) 
	ve_bit[ive] &= ~YP_BIT;
      else {
	assert( ~(ve_bit[ive] & YM_BIT) );
	ve_bit[ive] |= YM_BIT;
      }
    }
    void toggleYPBit(const int ive) {
      if (ve_bit[ive] & YM_BIT) 
	ve_bit[ive] &= ~YM_BIT;
      else {
	assert( ~(ve_bit[ive] & YP_BIT) );
	ve_bit[ive] |= YP_BIT;
      }
    }

    void resetBits() {
                  
      // if all bits are set, then triangle has completely migrated across periodic boundary...
      if ((ve_bit[0] & XM_BIT) && (ve_bit[1] & XM_BIT) && (ve_bit[2] & XM_BIT)) 
	FOR_IVE ve_bit[ive] &= ~XM_BIT;
      else if ((ve_bit[0] & XP_BIT) && (ve_bit[1] & XP_BIT) && (ve_bit[2] & XP_BIT)) 
	FOR_IVE ve_bit[ive] &= ~XP_BIT;
      if ((ve_bit[0] & YM_BIT) && (ve_bit[1] & YM_BIT) && (ve_bit[2] & YM_BIT)) 
	FOR_IVE ve_bit[ive] &= ~YM_BIT;
      else if ((ve_bit[0] & YP_BIT) && (ve_bit[1] & YP_BIT) && (ve_bit[2] & YP_BIT)) 
	FOR_IVE ve_bit[ive] &= ~YP_BIT;
    }
    
  };


  Quad<NQUAD> quad;
  
  vector<Tri> triVec;

  struct Point {
    double x[2];
    int ip;
    double tmp;
  };

  vector<list <Point> > cells;  
  int ncells;

  int np;
  double (*xp)[2];
  double *mp;
  double *rhop;
  double (*gp)[2];
  double (*up)[2];
  double *volp;
  double (*lambdap)[2];
  double *denomp;
  int *ip_flag;
  
  double *rho_error;

  double beta;
  double alpha;
  double eta;

  double gam;

  int step;
  int nsteps;
  double time;
  double dt;
  int check_interval;
  int write_interval;
    
  Paco() {
    cout << "Paco()" << endl;
    
    mp        = NULL;
    rhop      = NULL;
    gp        = NULL;
    up        = NULL;
    volp      = NULL;
    ip_flag   = NULL;
  }
  
  virtual ~Paco() {
    
    if (mp != NULL)      delete[] mp;
    if (rhop != NULL)    delete[] rhop;
    if (gp != NULL)      delete[] gp;
    if (up != NULL)      delete[] up;
    if (volp != NULL)    delete[] volp;
    if (ip_flag != NULL) delete[] ip_flag;
    
  }

  virtual void initialHook() {
    cout << "initialHook()" << endl;
    FOR_IP {
      mp[ip] = 0.0;
      FOR_I2 gp[ip][i] = 0.0;    
    }
  }

  virtual void temporalHook() {
    cout << "temporalHook()" << endl;
  }

  virtual void finalHook() {
    cout << "finalHook()" << endl;
  }
  
  void init() {
    cout << "init()" << endl;

    nsteps = 10000;
    check_interval = 10;
    write_interval = 100;

    beta = 5;
    alpha = 0.25;
    eta = 0.5;
            
    gam = 1.4;

    ncells = 4;
    list<Point> cell_list;
    for (int i = 0; i < ncells; ++i) 
      for (int j = 0; j < ncells; ++j) 
	cells.push_back(cell_list);
    
    buildCartMesh(128,128);   
    assert(mp == NULL);        mp = new double[np];
    assert(rhop == NULL);      rhop = new double[np];
    assert(gp == NULL);        gp = new double[np][2]; 
    assert(up == NULL);        up = new double[np][2];
    assert(volp == NULL);      volp = new double[np];
    assert(ip_flag == NULL);   ip_flag = new int[np];
    assert(lambdap == NULL);   lambdap = new double[np][2];
    assert(denomp == NULL);    denomp = new double[np];
    
    FOR_IP FOR_I2 lambdap[ip][i] = 0.0;

    rho_error = new double[np];

    updateCells();
    repairTris();    
    updateLambda();
    updateVol();
    
    //char filename[32];
    //sprintf(filename,"lp.%08d.dat",step);
    //writeLp(filename);

    initialHook();
    updateRho();
            
    step = 0;
    time = 0.0;
    
  }

  void updateAfterAdvect() {
    updateCells();
    repairTris();
    updateLambda();
    updateVol();
    updateRho();
  }
  
  
  void updateCells() {
    
    for (int i = 0; i < ncells; ++i) 
      for (int j = 0; j < ncells; ++j) 
	cells[i + ncells*j].clear();
    
    FOR_IP {
      if (xp[ip][0] < 0.0 && xp[ip][1] >= 0.0) {
	Point this_point;
	this_point.ip = ip;
	FOR_I2 this_point.x[i] = xp[ip][i];
	cells[1 + ncells*1].push_back(this_point);
	FOR_I2 this_point.x[i] = xp[ip][i];
	this_point.x[0] += (XMAX - XMIN);
	cells[3 + ncells*1].push_back(this_point);
	FOR_I2 this_point.x[i] = xp[ip][i];
	this_point.x[1] -= (YMAX - YMIN);
	cells[1 + ncells*3].push_back(this_point);
	FOR_I2 this_point.x[i] = xp[ip][i];
	this_point.x[0] += (XMAX - XMIN);
	this_point.x[1] -= (YMAX - YMIN);
	cells[3 + ncells*3].push_back(this_point);
      }
      else if (xp[ip][0] >= 0.0 && xp[ip][1] >= 0.0) {
	Point this_point;
	this_point.ip = ip;
	FOR_I2 this_point.x[i] = xp[ip][i];
	cells[2 + ncells*1].push_back(this_point);
	FOR_I2 this_point.x[i] = xp[ip][i];
	this_point.x[0] -= (XMAX - XMIN);
	cells[0 + ncells*1].push_back(this_point);
	FOR_I2 this_point.x[i] = xp[ip][i];
	this_point.x[1] -= (YMAX - YMIN);
	cells[2 + ncells*3].push_back(this_point);
	FOR_I2 this_point.x[i] = xp[ip][i];
	this_point.x[0] -= (XMAX - XMIN);
	this_point.x[1] -= (YMAX - YMIN);
	cells[0 + ncells*3].push_back(this_point);
      }
      else if (xp[ip][0] < 0.0 && xp[ip][1] < 0.0) {
	Point this_point;
	this_point.ip = ip;
	FOR_I2 this_point.x[i] = xp[ip][i];
	cells[1 + ncells*2].push_back(this_point);
	FOR_I2 this_point.x[i] = xp[ip][i];
	this_point.x[0] += (XMAX - XMIN);
	cells[3 + ncells*2].push_back(this_point);
	FOR_I2 this_point.x[i] = xp[ip][i];
	this_point.x[1] += (YMAX - YMIN);
	cells[1 + ncells*0].push_back(this_point);
	FOR_I2 this_point.x[i] = xp[ip][i];
	this_point.x[0] += (XMAX - XMIN);
	this_point.x[1] += (YMAX - YMIN);
	cells[3 + ncells*0].push_back(this_point);
      }
      else {
	assert(xp[ip][0] >= 0.0 && xp[ip][1] < 0.0);
	Point this_point;
	this_point.ip = ip;
	FOR_I2 this_point.x[i] = xp[ip][i];
	cells[2 + ncells*2].push_back(this_point);
	FOR_I2 this_point.x[i] = xp[ip][i];
	this_point.x[0] -= (XMAX - XMIN);
	cells[0 + ncells*2].push_back(this_point);
	FOR_I2 this_point.x[i] = xp[ip][i];
	this_point.x[1] += (YMAX - YMIN);
	cells[2 + ncells*0].push_back(this_point);
	FOR_I2 this_point.x[i] = xp[ip][i];
	this_point.x[0] -= (XMAX - XMIN);
	this_point.x[1] += (YMAX - YMIN);
	cells[0 + ncells*0].push_back(this_point);
      }
    }
    
    /*
    char filename[32];
    sprintf(filename,"lp.%08d.dat",step);
    writeLp(filename);
    */
    
  }
  
  
  void writeLp(char * filename) {
    
    cout << "writing " << filename << endl;
    
    FILE * fp = fopen(filename,"w");
      
    int count = 0;
    FOR_IT {
      FOR_IQ {
	++count;
      }
    }
  
    fprintf(fp,"TITLE = \"%s\"\n",filename);
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"LAMX\"\n");
    fprintf(fp,"\"LAMY\"\n");
    
    fprintf(fp,"ZONE T=\"%s, np=%d\" I=%d, J=1, K=1, ZONETYPE=Ordered, DATAPACKING=POINT\n",filename,count,count);
    FOR_IT {
      triVec[it].updateVeX(xp);
      FOR_IQ {
	double xquad[2] = {0.0, 0.0};
	FOR_I2 FOR_IVE xquad[i] += triVec[it].ve_x[ive][i]*quad.ve_wgt[iq][ive];
	fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n", xquad[0], xquad[1], triVec[it].lambda[iq][0], triVec[it].lambda[iq][1]);
      }
    }
    
    fclose(fp);
    
  }


  void writeLp2(char * filename) {
    
    cout << "writing " << filename << endl;
    
    FILE * fp = fopen(filename,"w");
    
    int it = 100;
    int iq = 1;
    
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
    
    cout << icell << " " << jcell << endl;

    int count = 0;
    for (int ivar = -1; ivar <= 1; ++ivar) {
      for (int jvar = -1; jvar <= 1; ++jvar) {
	const int this_icell = icell + ivar;
	const int this_jcell = jcell + jvar;
	for (typename list<Point>::iterator il = cells[this_icell + ncells*this_jcell].begin(); il != cells[this_icell + ncells*this_jcell].end(); ++il) {
	  ++count;
	}
      }
    }
      
    fprintf(fp,"TITLE = \"%s\"\n",filename);
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"LAMX\"\n");
    fprintf(fp,"\"LAMY\"\n");
    
    fprintf(fp,"ZONE T=\"%s, np=%d\" I=%d, J=1, K=1, ZONETYPE=Ordered, DATAPACKING=POINT\n",filename,count+1,count+1);
        
    fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n", xquad[0], xquad[1], 1.0, 1.0);

    for (int ivar = -1; ivar <= 1; ++ivar) {
      for (int jvar = -1; jvar <= 1; ++jvar) {
	const int this_icell = icell + ivar;
	const int this_jcell = jcell + jvar;
	for (typename list<Point>::iterator il = cells[this_icell + ncells*this_jcell].begin(); il != cells[this_icell + ncells*this_jcell].end(); ++il) {
	  double tmp[2] = {0.0, 0.0};
	  const double val = evalExp(xquad, il->x, tmp); //triVec[it].lambda[iq]);
	  fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n", il->x[0], il->x[1], val, val);
	}
      }
    }
    
    fclose(fp);
    
  }


  void buildCartMesh(const int nx, const int ny) {
    
    const double dx = (XMAX - XMIN)/(double) nx;
    const double dy = (YMAX - YMIN)/(double) ny;
    
    // way to set to const?
    np = nx*ny;

    // allocate and populate coordinates... 
    xp = new double[np][2];
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
	xp[i*ny + j][0] = XMIN + i*dx + 0.5*dx - 0.5*(j%2)*dx;
	xp[i*ny + j][1] = YMIN + j*dy + 0.5*dy;
      }
    }
    
    // build tri's manually...
    // first interior...
    for (int i = 0; i < (nx-1); ++i) {
      for (int j = 0; j < (ny-1); ++j) {
	// "upward" pointing
	Tri tri1;
	// right handed tet..
	tri1.ve_ip[0] = i*ny + j;
	tri1.ve_ip[1] = (i+1)*ny + j;
	tri1.ve_ip[2] = (i+1)*ny + (j+1);
	
	// no periodicity...
	tri1.ve_bit[0] = 0;
	tri1.ve_bit[1] = 0;
	tri1.ve_bit[2] = 0;
	
	triVec.push_back(tri1);
	
	// "downward" pointing
	Tri tri2;
	// right handed tet..
	tri2.ve_ip[0] = i*ny + j;
	tri2.ve_ip[1] = (i+1)*ny + (j+1);
	tri2.ve_ip[2] = i*ny + (j+1);
	
	// no periodicity...
	tri2.ve_bit[0] = 0;
	tri2.ve_bit[1] = 0;
	tri2.ve_bit[2] = 0;
	
	triVec.push_back(tri2);     
	
      }
    }
    
    // now build periodic boundaries...
    // x boundaries..
    for (int j = 0; j < (ny-1); ++j) {
      // "upward" pointing
      Tri tri1;
      // right handed tet..
      tri1.ve_ip[0] = j;
      tri1.ve_ip[1] = j+1;
      tri1.ve_ip[2] = (nx-1)*ny + j;
      
      tri1.ve_bit[0] = 0;
      tri1.ve_bit[1] = 0;
      tri1.ve_bit[2] = XM_BIT;
    
      triVec.push_back(tri1);
      
      // "downward" pointing
      Tri tri2;
      // right handed tet..
      tri2.ve_ip[0] = (nx-1)*ny + j;
      tri2.ve_ip[1] = j+1;
      tri2.ve_ip[2] = (nx-1)*ny + (j+1);
      
      tri2.ve_bit[0] = 0;
      tri2.ve_bit[1] = XP_BIT;
      tri2.ve_bit[2] = 0;
      
      triVec.push_back(tri2);
    }
    
    // y boundaries..
    for (int i = 0; i < (nx-1); ++i) {
      // "upward" pointing
      Tri tri1;
      // right handed tet..
      tri1.ve_ip[0] = i*ny + (ny-1);
      tri1.ve_ip[1] = (i+1)*ny + (ny-1);
      tri1.ve_ip[2] = (i+1)*ny;
    
      tri1.ve_bit[0] = 0;
      tri1.ve_bit[1] = 0;
      tri1.ve_bit[2] = YP_BIT;
      
      triVec.push_back(tri1);
      
      // "downward" pointing
      Tri tri2;
      // right handed tet..
      tri2.ve_ip[0] = i*ny;
      tri2.ve_ip[1] = i*ny + (ny-1);
      tri2.ve_ip[2] = (i+1)*ny;;
      
      tri2.ve_bit[0] = 0;
      tri2.ve_bit[1] = YM_BIT;
      tri2.ve_bit[2] = 0;
      
      triVec.push_back(tri2);
    }
  
    // and finally the cornors...
    {
      // xm, ym cornor
      Tri tri1;
      tri1.ve_ip[0] = 0;
      tri1.ve_ip[1] = ny*(nx-1) + (ny-1);
      tri1.ve_ip[2] = (ny-1);
      
      tri1.ve_bit[0] = 0;
      tri1.ve_bit[1] = XM_BIT|YM_BIT;
      tri1.ve_bit[2] = YM_BIT;
      
      triVec.push_back(tri1);
      
      Tri tri2;
      tri2.ve_ip[0] = 0;
      tri2.ve_ip[1] = ny*(nx-1);
      tri2.ve_ip[2] = ny*(nx-1) + (ny-1);
      
      tri2.ve_bit[0] = 0;
      tri2.ve_bit[1] = XM_BIT;
      tri2.ve_bit[2] = XM_BIT|YM_BIT;
      
      triVec.push_back(tri2);    
    }
    
    
    // build neighbors.. (we only do this once, so efficiency doesn't matter)
    FOR_IT {
      FOR_IVE {
	bool found = false;
	for (int it_nbr = 0; it_nbr < triVec.size(); ++it_nbr) {
	  if (it_nbr != it) {
	    for (int ive_nbr = 0; ive_nbr < 3; ++ive_nbr) {
	      if (triVec[it].ve_ip[(ive+1)%3] == triVec[it_nbr].ve_ip[ive_nbr]) {
		if (triVec[it].ve_ip[(ive+2)%3] == triVec[it_nbr].ve_ip[(ive_nbr+1)%3]) {
		  found = true;
		  triVec[it].ve_itnbr[ive] = it_nbr;
		  triVec[it].ve_ivenbr[ive] = (ive_nbr+2)%3;	      
		  break;
		}
		else if (triVec[it].ve_ip[(ive+2)%3] == triVec[it_nbr].ve_ip[(ive_nbr+2)%3]) {
		  found = true;
		  triVec[it].ve_itnbr[ive] = it_nbr;
		  triVec[it].ve_ivenbr[ive] = (ive_nbr+1)%3;	      
		  break;
		} 
	      }
	    }
	    if (found == true) break; 
	  }
	}						
	assert(found);
      }
    }
    
  }
  
  void flipEdges() {
    
    bool done = false;
    int iter = 0;
    int nflip = 0;
    while (!done) {
      done = true;
      assert(iter < 100);
      FOR_IT {
	FOR_IVE {
	  
	  int it_nbr  = triVec[it].ve_itnbr[ive];
	  int ive_nbr = triVec[it].ve_ivenbr[ive];
	  
	  assert(it != it_nbr);
	  
	  triVec[it].updateVeX(xp);
	  triVec[it_nbr].updateVeX(xp);
	  
	  /*
	    cout << "it verticies: ";
	    FOR_I3 cout << triVec[it].ve_x[i][0] << " " << triVec[it].ve_x[i][1] << " " ;
	    cout << endl;
	    FOR_I3 cout << triVec[it].ve_ip[i] << " ";
	    cout << endl;
	    FOR_I3 cout << triVec[it].ve_bit[i] << " ";
	    cout << endl;
	    cout << "it_nbr verticies: ";
	    FOR_I3 cout << triVec[it_nbr].ve_x[i][0] << " " << triVec[it_nbr].ve_x[i][1] << " ";
	    cout << endl;
	    FOR_I3 cout << triVec[it_nbr].ve_ip[i] << " ";
	    cout << endl;
	    FOR_I3 cout << triVec[it_nbr].ve_bit[i] << " ";
	    cout << endl;
	  */
	  
	  int offset;
	  // transform bits so they are equal on verticies they share
	  if ( triVec[it].ve_ip[(ive+1)%3] == triVec[it_nbr].ve_ip[(ive_nbr+1)%3]) {
	    offset = 0;
	  }
	  else {
	    assert(triVec[it].ve_ip[(ive+1)%3] == triVec[it_nbr].ve_ip[(ive_nbr+2)%3]);
	    offset = 1;
	  }
	  
	  assert(triVec[it].ve_ip[(ive+1)%3] == triVec[it_nbr].ve_ip[(ive_nbr+offset+1)%3]);
	  assert(triVec[it].ve_ip[(ive+2)%3] == triVec[it_nbr].ve_ip[(ive_nbr-offset+2)%3]);
	
	  if ( (triVec[it].ve_bit[(ive+1)%3] & XM_BIT) != (triVec[it_nbr].ve_bit[(ive_nbr+offset+1)%3] & XM_BIT)) {
	    if (triVec[it].ve_bit[(ive+1)%3] & XM_BIT) {
	      FOR_I3 triVec[it].toggleXPBit(i);
	    }
	    else {
	      FOR_I3 triVec[it_nbr].toggleXPBit(i);
	    }
	  }
	  else if ( (triVec[it].ve_bit[(ive+1)%3] & XP_BIT) != (triVec[it_nbr].ve_bit[(ive_nbr+offset+1)%3] & XP_BIT)) {
	    if (triVec[it].ve_bit[(ive+1)%3] & XP_BIT) {
	      FOR_I3 triVec[it].toggleXMBit(i);
	    }
	    else {
	      FOR_I3 triVec[it_nbr].toggleXMBit(i);
	    }
	  }
	  if ( (triVec[it].ve_bit[(ive+1)%3] & YM_BIT) != (triVec[it_nbr].ve_bit[(ive_nbr+offset+1)%3] & YM_BIT)) {
	    if (triVec[it].ve_bit[(ive+1)%3] & YM_BIT) {
	      FOR_I3 triVec[it].toggleYPBit(i);
	    }
	    else {
	      FOR_I3 triVec[it_nbr].toggleYPBit(i);
	    }
	  }
	  else if ( (triVec[it].ve_bit[(ive+1)%3] & YP_BIT) != (triVec[it_nbr].ve_bit[(ive_nbr+offset+1)%3] & YP_BIT)) {
	    if (triVec[it].ve_bit[(ive+1)%3] & YP_BIT) {
	    FOR_I3 triVec[it].toggleYMBit(i);
	    }
	    else {
	      FOR_I3 triVec[it_nbr].toggleYMBit(i);
	    }
	  }
	  
	  if ( (triVec[it].ve_bit[(ive+2)%3] & XM_BIT) != (triVec[it_nbr].ve_bit[(ive_nbr-offset+2)%3] & XM_BIT)) {
	    if (triVec[it].ve_bit[(ive+2)%3] & XM_BIT) {
	      FOR_I3 triVec[it].toggleXPBit(i);
	    }
	    else {
	      FOR_I3 triVec[it_nbr].toggleXPBit(i);
	    }
	  }
	  else if ( (triVec[it].ve_bit[(ive+2)%3] & XP_BIT) != (triVec[it_nbr].ve_bit[(ive_nbr-offset+2)%3] & XP_BIT)) {
	    if (triVec[it].ve_bit[(ive+2)%3] & XP_BIT) {
	    FOR_I3 triVec[it].toggleXMBit(i);
	    }
	    else {
	      FOR_I3 triVec[it_nbr].toggleXMBit(i);
	    }
	  }
	  if ( (triVec[it].ve_bit[(ive+2)%3] & YM_BIT) != (triVec[it_nbr].ve_bit[(ive_nbr-offset+2)%3] & YM_BIT)) {
	    if (triVec[it].ve_bit[(ive+2)%3] & YM_BIT) {
	      FOR_I3 triVec[it].toggleYPBit(i);
	    }
	    else {
	      FOR_I3 triVec[it_nbr].toggleYPBit(i);
	    }
	  }
	  else if ( (triVec[it].ve_bit[(ive+2)%3] & YP_BIT) != (triVec[it_nbr].ve_bit[(ive_nbr-offset+2)%3] & YP_BIT)) {
	    if (triVec[it].ve_bit[(ive+2)%3] & YP_BIT) {
	      FOR_I3 triVec[it].toggleYMBit(i);
	    }
	    else {
	      FOR_I3 triVec[it_nbr].toggleYMBit(i);
	    }
	  }
	  
	  triVec[it].updateVeX(xp);
	  triVec[it_nbr].updateVeX(xp);
	  
	  double xc[2];
	  {
	    double bx = triVec[it].ve_x[1][0] - triVec[it].ve_x[0][0];
	    double by = triVec[it].ve_x[1][1] - triVec[it].ve_x[0][1];
	    double cx = triVec[it].ve_x[2][0] - triVec[it].ve_x[0][0];
	    double cy = triVec[it].ve_x[2][1] - triVec[it].ve_x[0][1];
	    double d = 2.0*(bx*cy - by*cx);
	    xc[0] = (cy*(bx*bx+by*by) - by*(cx*cx+cy*cy))/d + triVec[it].ve_x[0][0];
	    xc[1] = (bx*(cx*cx+cy*cy) - cx*(bx*bx+by*by))/d + triVec[it].ve_x[0][1];
	  }
	  
	  
	  double r_it = sqrt((triVec[it].ve_x[ive][0] - xc[0])*(triVec[it].ve_x[ive][0] - xc[0]) + 
			     (triVec[it].ve_x[ive][1] - xc[1])*(triVec[it].ve_x[ive][1] - xc[1])); 
	  double r_it_nbr = sqrt((triVec[it_nbr].ve_x[ive_nbr][0] - xc[0])*(triVec[it_nbr].ve_x[ive_nbr][0] - xc[0]) + 
				 (triVec[it_nbr].ve_x[ive_nbr][1] - xc[1])*(triVec[it_nbr].ve_x[ive_nbr][1] - xc[1])); 
	  
	  if ( (r_it - r_it_nbr) > 1e-9*r_it ) {
	    
	    // need temporary tri to store old info...
	    Tri tmp = triVec[it];
	    Tri tmp_nbr = triVec[it_nbr];
	    
	    triVec[it].ve_ip[(ive+1)%3] = tmp_nbr.ve_ip[ive_nbr];
	    triVec[it].ve_bit[(ive+1)%3] = tmp_nbr.ve_bit[ive_nbr];
	    
	    triVec[it].ve_itnbr[ive] =  tmp_nbr.ve_itnbr[(ive_nbr+offset+1)%3];
	    triVec[it].ve_ivenbr[ive] =  tmp_nbr.ve_ivenbr[(ive_nbr+offset+1)%3];
	    triVec[triVec[it].ve_itnbr[ive]].ve_itnbr[triVec[it].ve_ivenbr[ive]] = it;
	    triVec[triVec[it].ve_itnbr[ive]].ve_ivenbr[triVec[it].ve_ivenbr[ive]] = ive;
	    
	    triVec[it].ve_itnbr[(ive+2)%3] =  it_nbr;
	    triVec[it].ve_ivenbr[(ive+2)%3] =  (ive_nbr+offset+1)%3;  
	    
	    triVec[it_nbr].ve_ip[(ive_nbr-offset+2)%3] = tmp.ve_ip[ive];
	    triVec[it_nbr].ve_bit[(ive_nbr-offset+2)%3] = tmp.ve_bit[ive];
	    
	    triVec[it_nbr].ve_itnbr[ive_nbr] =  tmp.ve_itnbr[(ive+2)%3];
	    triVec[it_nbr].ve_ivenbr[ive_nbr] =  tmp.ve_ivenbr[(ive+2)%3];
	    triVec[triVec[it_nbr].ve_itnbr[ive_nbr]].ve_itnbr[triVec[it_nbr].ve_ivenbr[ive_nbr]] = it_nbr;
	    triVec[triVec[it_nbr].ve_itnbr[ive_nbr]].ve_ivenbr[triVec[it_nbr].ve_ivenbr[ive_nbr]] = ive_nbr;
	    
	    triVec[it_nbr].ve_itnbr[(ive_nbr+offset+1)%3] =  it;
	    triVec[it_nbr].ve_ivenbr[(ive_nbr+offset+1)%3] =  (ive+2)%3;
	    
	    done = false;
	    ++nflip;
	    
	  }
	  
	  triVec[it].updateVeX(xp);
	  triVec[it_nbr].updateVeX(xp);
	  
	  /*
	    cout << "it verticies: ";
	    FOR_I3 cout << triVec[it].ve_x[i][0] << " " << triVec[it].ve_x[i][1] << " " ;
	    cout << endl;
	    FOR_I3 cout << triVec[it].ve_ip[i] << " ";
	    cout << endl;
	    FOR_I3 cout << triVec[it].ve_bit[i] << " ";
	    cout << endl;
	    cout << "nbrs: ";
	    FOR_I3 cout << triVec[triVec[it].ve_itnbr[0]].ve_ip[i] << " " ;
	    cout << endl;
	    cout << triVec[triVec[it].ve_itnbr[0]].ve_ip[triVec[it].ve_ivenbr[0]] << endl;
	    cout << "nbrs: ";
	    FOR_I3 cout << triVec[triVec[it].ve_itnbr[1]].ve_ip[i] << " " ;
	    cout << endl;
	    cout << triVec[triVec[it].ve_itnbr[1]].ve_ip[triVec[it].ve_ivenbr[1]] << endl;
	    cout << "nbrs: ";
	    FOR_I3 cout << triVec[triVec[it].ve_itnbr[2]].ve_ip[i] << " " ;
	    cout << endl;
	    cout << triVec[triVec[it].ve_itnbr[2]].ve_ip[triVec[it].ve_ivenbr[2]] << endl;
	    cout << "it_nbr verticies: ";
	    FOR_I3 cout << triVec[it_nbr].ve_x[i][0] << " " << triVec[it_nbr].ve_x[i][1] << " ";
	    cout << endl;
	    FOR_I3 cout << triVec[it_nbr].ve_ip[i] << " ";
	    cout << endl;
	    FOR_I3 cout << triVec[it_nbr].ve_bit[i] << " ";
	    cout << endl;	
	    cout << "nbrs: ";
	    FOR_I3 cout << triVec[triVec[it_nbr].ve_itnbr[0]].ve_ip[i] << " " ;
	    cout << endl;	
	    cout << triVec[triVec[it_nbr].ve_itnbr[0]].ve_ip[triVec[it_nbr].ve_ivenbr[0]] << endl;
	    cout << "nbrs: ";
	    FOR_I3 cout << triVec[triVec[it_nbr].ve_itnbr[1]].ve_ip[i] << " " ;
	    cout << endl;
	    cout << triVec[triVec[it_nbr].ve_itnbr[1]].ve_ip[triVec[it_nbr].ve_ivenbr[1]] << endl;
	    cout << "nbrs: ";
	    FOR_I3 cout << triVec[triVec[it_nbr].ve_itnbr[2]].ve_ip[i] << " " ;
	    cout << endl;
	    cout << triVec[triVec[it_nbr].ve_itnbr[2]].ve_ip[triVec[it_nbr].ve_ivenbr[2]] << endl;
	  */
	  
	  triVec[it].resetBits();
	  triVec[it_nbr].resetBits();
	  
	  
	}
      }
      ++iter;
    }
    
    
  }

    
  inline double evalExp(const double _x[2], const double _y[2], const double _lam[2]) { 
    return(exp(-beta*((_x[0] - _y[0])*(_x[0] - _y[0]) + (_x[1] - _y[1])*(_x[1] - _y[1])) + _lam[0]*(_x[0] - _y[0]) + _lam[1]*(_x[1] - _y[1]) ));
  }
  /*
  void updateLambda() {
    
    // first quadrature points...    
    #pragma omp parallel for
    FOR_IT {
      triVec[it].updateVeX(xp);
      
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
	
	double g[2];
	double H[2][2];
	double denom;
	double lambda[2] = {triVec[it].lambda[iq][0], triVec[it].lambda[iq][1]};
	int iter = 0;
	bool done = false;
	while(1) {
	  FOR_I2 g[i] = 0.0;
	  FOR_I2 FOR_J2 H[i][j] = 0.0;
	  denom = 0.0;
	  
	  for (int ivar = -1; ivar <= 1; ++ivar) {
	    for (int jvar = -1; jvar <= 1; ++jvar) {
	      const int this_icell = icell + ivar;
	      const int this_jcell = jcell + jvar;
	      for (typename list<Point>::iterator il = cells[this_icell + ncells*this_jcell].begin(); il != cells[this_icell + ncells*this_jcell].end(); ++il) {
		const double val = evalExp(xquad, il->x, lambda);
		denom += val;
		FOR_I2 g[i] += (xquad[i] - il->x[i])*val;
		FOR_I2 FOR_J2 H[i][j] += (xquad[i] - il->x[i])*(xquad[j] - il->x[j])*val;
	      }
	    }
	  }
	  
	  FOR_I2 g[i] /= denom;
	  FOR_I2 FOR_J2 H[i][j] = H[i][j]/denom - g[i]*g[j];
	  
	  if (done) {
	    triVec[it].denom[iq] = denom;
	    FOR_I2 triVec[it].g[iq][i] = g[i];
	    FOR_I2 FOR_J2 triVec[it].H[iq][i][j] = H[i][j];
	    FOR_I2 triVec[it].lambda[iq][i] = lambda[i];
	    break;
	  }
	  
	  const double discrim = H[0][0]*H[1][1] - H[0][1]*H[1][0];
	  const double dl[2] = {-(g[0]*H[1][1]-g[1]*H[1][0])/discrim, -(-g[0]*H[0][1] + g[1]*H[0][0])/discrim};
	  
	  double st = 1.0;
	  double lambda_st[2];
	  FOR_I2 lambda_st[i] = lambda[i] + st*dl[i];
	  double denom_st = 0.0;
	  
	  for (int ivar = -1; ivar <= 1; ++ivar) {
	    for (int jvar = -1; jvar <= 1; ++jvar) {
	      const int this_icell = icell + ivar;
	      const int this_jcell = jcell + jvar;
	      for (typename list<Point>::iterator il = cells[this_icell + ncells*this_jcell].begin(); il != cells[this_icell + ncells*this_jcell].end(); ++il) {
		denom_st += evalExp(xquad, il->x, lambda_st);
	      }
	    }
	  }
	  
	  while(log(denom_st) > log(denom) + alpha*st*(g[0]*dl[0] + g[1]*dl[1])) {
	    st *= eta;
	    denom_st = 0.0;
	    FOR_I2 lambda_st[i] = lambda[i] + st*dl[i];
	    
	    for (int ivar = -1; ivar <= 1; ++ivar) {
	      for (int jvar = -1; jvar <= 1; ++jvar) {
		const int this_icell = icell + ivar;
		const int this_jcell = jcell + jvar;
		for (typename list<Point>::iterator il = cells[this_icell + ncells*this_jcell].begin(); il != cells[this_icell + ncells*this_jcell].end(); ++il) {
		  denom_st += evalExp(xquad, il->x, lambda_st);
		}
	      }
	    }
	  }
	  
	  FOR_I2 lambda[i] += st*dl[i];
	  ++iter;
	  
	  assert( iter < 50);
	  
	  // should write a NaN checker...
	  FOR_I2 assert( lambda[i] == lambda[i] );
	  
	  if (fabs(g[0]*dl[0] + g[1]*dl[1]) < 1e-12) done = true;
	}
      }
      
    }
    
    
  }
  */
  
  void updateLambda() {
        
    // first quadrature points...    
    #pragma omp parallel for
    FOR_IT {
      triVec[it].updateVeX(xp);
      
      FOR_IQ {
	double xquad[2] = {0.0, 0.0};
	FOR_I2 FOR_IVE xquad[i] += triVec[it].ve_x[ive][i]*quad.ve_wgt[iq][ive];
	
	calcLambda(triVec[it].lambda[iq], triVec[it].denom[iq], triVec[it].g[iq], triVec[it].H[iq], xquad);
      }
    }
    
    // then particles... (i.e. verticies)
    FOR_IP {
      double g[2];
      double H[2][2];
      calcLambda(lambdap[ip], denomp[ip], g, H, xp[ip]);
    }
  
  }
    
  void calcLambda(double this_lambda[2], double &this_denom, double this_g[2], double this_H[2][2], const double this_x[2]) {

    // make sure lambda has something in it!
    
    int icell, jcell;
    if (this_x[0] < 0.0 && this_x[1] >= 0.0) {
      icell = 1;
      jcell = 1;
    }
    else if (this_x[0] >= 0.0 && this_x[1] >= 0.0) {
      icell = 2;
      jcell = 1;
    }
    else if (this_x[0] < 0.0 && this_x[1] < 0.0) {
      icell = 1;
      jcell = 2;
    }
    else {
      assert(this_x[0] >= 0.0 && this_x[1] < 0.0);
      icell = 2;
      jcell = 2;
    }
    
    int iter = 0;
    bool done = false;
    while(1) {
      FOR_I2 this_g[i] = 0.0;
      FOR_I2 FOR_J2 this_H[i][j] = 0.0;
      this_denom = 0.0;
      
      for (int ivar = -1; ivar <= 1; ++ivar) {
	for (int jvar = -1; jvar <= 1; ++jvar) {
	  const int this_icell = icell + ivar;
	  const int this_jcell = jcell + jvar;
	  for (typename list<Point>::iterator il = cells[this_icell + ncells*this_jcell].begin(); il != cells[this_icell + ncells*this_jcell].end(); ++il) {
	    const double val = evalExp(this_x, il->x, this_lambda);
	    this_denom += val;
	    FOR_I2 this_g[i] += (this_x[i] - il->x[i])*val;
	    FOR_I2 FOR_J2 this_H[i][j] += (this_x[i] - il->x[i])*(this_x[j] - il->x[j])*val;
	  }
	}
      }
      
      FOR_I2 this_g[i] /= this_denom;
      FOR_I2 FOR_J2 this_H[i][j] = this_H[i][j]/this_denom - this_g[i]*this_g[j];
      
      if (done) break; 
	  
      const double discrim = this_H[0][0]*this_H[1][1] - this_H[0][1]*this_H[1][0];
      const double dl[2] = {-(this_g[0]*this_H[1][1]-this_g[1]*this_H[1][0])/discrim, -(-this_g[0]*this_H[0][1] + this_g[1]*this_H[0][0])/discrim};
      
      double st = 1.0;
      double lambda_st[2];
      FOR_I2 lambda_st[i] = this_lambda[i] + st*dl[i];
      double denom_st = 0.0;
      
      for (int ivar = -1; ivar <= 1; ++ivar) {
	for (int jvar = -1; jvar <= 1; ++jvar) {
	  const int this_icell = icell + ivar;
	  const int this_jcell = jcell + jvar;
	  for (typename list<Point>::iterator il = cells[this_icell + ncells*this_jcell].begin(); il != cells[this_icell + ncells*this_jcell].end(); ++il) {
	    denom_st += evalExp(this_x, il->x, lambda_st);
	  }
	}
      }
      
      while(log(denom_st) > log(this_denom) + alpha*st*(this_g[0]*dl[0] + this_g[1]*dl[1])) {
	st *= eta;
	denom_st = 0.0;
	FOR_I2 lambda_st[i] = this_lambda[i] + st*dl[i];
	
	for (int ivar = -1; ivar <= 1; ++ivar) {
	  for (int jvar = -1; jvar <= 1; ++jvar) {
	    const int this_icell = icell + ivar;
	    const int this_jcell = jcell + jvar;
	    for (typename list<Point>::iterator il = cells[this_icell + ncells*this_jcell].begin(); il != cells[this_icell + ncells*this_jcell].end(); ++il) {
	      denom_st += evalExp(this_x, il->x, lambda_st);
	    }
	  }
	}
      }
      
      FOR_I2 this_lambda[i] += st*dl[i];
      ++iter;
      
      assert( iter < 50);
      
      // should write a NaN checker...
      FOR_I2 assert( this_lambda[i] == this_lambda[i] );
      
      if (fabs(this_g[0]*dl[0] + this_g[1]*dl[1]) < 1e-12) done = true;
    }
    
  }
  
  void updateVol() {
        
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
		
	const double lambda[2] = {triVec[it].lambda[iq][0],triVec[it].lambda[iq][1]};
	const double weight = quad.p_wgt[iq]*triVec[it].area/triVec[it].denom[iq];
	
        #pragma omp parallel for collapse(2)
	for (int ivar = -1; ivar <= 1; ++ivar) {
	  for (int jvar = -1; jvar <= 1; ++jvar) {
	    const int this_icell = icell + ivar;
	    const int this_jcell = jcell + jvar;
	    for (typename list<Point>::iterator il = cells[this_icell + ncells*this_jcell].begin(); il != cells[this_icell + ncells*this_jcell].end(); ++il) {
	      il->tmp += weight*evalExp(xquad, il->x, lambda);
	    }
	  }
	}
	
      }
    }
        
    FOR_IP volp[ip] = 0.0;

    for (int i = 0; i < ncells; ++i) {
      for (int j = 0; j < ncells; ++j) {
	for (typename list<Point>::iterator il = cells[i + ncells*j].begin(); il != cells[i + ncells*j].end(); ++il) {
	  volp[il->ip] += il->tmp;
	}
      }
    }

    
    double vol_sum = 0.0;
    FOR_IP vol_sum += volp[ip];
    assert( fabs(vol_sum - (XMAX-XMIN)*(YMAX-YMIN)) < 1e-9); 
    //cout << "vol_sum: " << vol_sum << endl;
    
  }

  void setDt() {
    dt = 0.02;

  }
  
  void repairTris() {
    
    // flip points that have crossed periodic boundaries..
    FOR_IP {
      // store flipping info in ip_flag...
      ip_flag[ip] = 0;
      if (xp[ip][0] > XMAX) {
	xp[ip][0] -= (XMAX - XMIN);
	ip_flag[ip] |= XP_BIT;
      }
      else if (xp[ip][0] < XMIN) {
	xp[ip][0] += (XMAX - XMIN);
	ip_flag[ip] |= XM_BIT;
      }
      
      if (xp[ip][1] > YMAX) {
	xp[ip][1] -= (YMAX - YMIN);
	ip_flag[ip] |= YP_BIT;
      }
      else if (xp[ip][1] < YMIN) {
	xp[ip][1] += (YMAX - YMIN);
	ip_flag[ip] |= YM_BIT;
      }
    }
   
    FOR_IT {
      FOR_IVE {
	if      (ip_flag[triVec[it].ve_ip[ive]] & XM_BIT) triVec[it].toggleXMBit(ive);
	else if (ip_flag[triVec[it].ve_ip[ive]] & XP_BIT) triVec[it].toggleXPBit(ive);
	if      (ip_flag[triVec[it].ve_ip[ive]] & YM_BIT) triVec[it].toggleYMBit(ive);
	else if (ip_flag[triVec[it].ve_ip[ive]] & YP_BIT) triVec[it].toggleYPBit(ive);
      }
    }
    
    
    FOR_IT triVec[it].resetBits();
  
    flipEdges();
    
  }
  
  virtual void updateU(double (*_u)[2]) {
    
  }

  void updateRho() {
    
    FOR_IP rhop[ip] = mp[ip]/volp[ip];
    
  }

  void calcRhs(double (*rhs_gp)[2]) {
    
    FOR_IP FOR_I2 rhs_gp[ip][i] = 1.0;
    
  }
  
  void run() {
    cout << "run()" << endl;
    
    temporalHook();
    
    char filename[32];
    sprintf(filename,"tri.%08d.dat",step);
    writeTecplot(filename);

    double (*dg)[2] = new double[np][2];
    
    double (*x0)[2] = new double[np][2];
    double (*dx0)[2] = new double[np][2];
    double (*dx1)[2] = new double[np][2];
    double (*dx2)[2] = new double[np][2];

    bool done = false;
    while(!done) {
      
      ++step;
      setDt();
      time += dt;

      if (step%check_interval == 0) {
	cout <<
	  "\n----------------------------------------------------------\n" <<
	  " starting step: " << step << " time: " << time << " dt: " << dt <<
	  "\n----------------------------------------------------------" << endl;
      }
      /*
      if (step == 1) calcRhs(dg);  // this is only the case for Euler eqns. - viscous term relies on velocities and therefore requires two RHS calcs per timestep

      FOR_IP FOR_I2 gp[ip][i] += 0.5*dt*dg[ip][i];
      updateU();
      FOR_IP FOR_I2 xp[ip][i] += dt*up[ip][i];
      updateAfterAdvect();
      calcRhs(dg);
      FOR_IP FOR_I2 gp[ip][i] += 0.5*dt*dg[ip][i];
      */

      FOR_IP FOR_I2 x0[ip][i] = xp[ip][i];
      updateU(dx0);
      FOR_IP FOR_I2 xp[ip][i] = x0[ip][i] + dt*dx0[ip][i];
      //updateAfterAdvect();
      updateU(dx1);
      FOR_IP FOR_I2 xp[ip][i] = x0[ip][i] + dt*(dx0[ip][i] + dx1[ip][i])/4.0;
      //updateAfterAdvect();
      updateU(dx2);
      FOR_IP FOR_I2 xp[ip][i] = x0[ip][i] + dt*(dx0[ip][i] + dx1[ip][i] + 4.0*dx2[ip][i])/6.0;
      

      if (step%100 == 0) {
	updateAfterAdvect();    
	temporalHook();
      }
      
      if (step%write_interval == 0) {
	char filename[32];
	sprintf(filename,"tri.%08d.dat",step);
	writeTecplot(filename);
      }

      
      if (step >= nsteps) done = true;
      
    }
    
    delete[] dg;
    
  }
  
  
  void writeTecplot(char * filename) {
    
    cout << "writing " << filename << endl;
    
    FILE * fp = fopen(filename,"w");
    
    int tri_count = 0;
    FOR_IT {
      if ((triVec[it].ve_bit[0] == 0) && (triVec[it].ve_bit[1] == 0) && (triVec[it].ve_bit[2] == 0))
      ++tri_count;
    }
    
    fprintf(fp,"TITLE = \"%s\"\n",filename);
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"VOL\"\n");
    fprintf(fp,"\"M\"\n");
    fprintf(fp,"\"RHO\"\n");
    fprintf(fp,"\"RHO_ERROR\"\n");
    fprintf(fp,"\"U-X\"\n");
    fprintf(fp,"\"U_Y\"\n");
    
    fprintf(fp,"ZONE T=\"%s\"\n",filename);
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", np, tri_count);
    
    // back out pointwise quantity 
    double * my_rho = new double[np]; 
    
    #pragma omp parallel for 
    FOR_IP {
      my_rho[ip] = 0.0;
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
      
      const double coeff = 1.0/denomp[ip];
            	
      for (int ivar = -1; ivar <= 1; ++ivar) {
	for (int jvar = -1; jvar <= 1; ++jvar) {
	  const int this_icell = icell + ivar;
	  const int this_jcell = jcell + jvar;
	  for (typename list<Point>::iterator il = cells[this_icell + ncells*this_jcell].begin(); il != cells[this_icell + ncells*this_jcell].end(); ++il) {
	    my_rho[ip] += coeff*rhop[il->ip]*evalExp(xp[ip], il->x, lambdap[ip]);
	  }
	}
      }
    }

    FOR_IP  {      
      fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n", xp[ip][0], xp[ip][1], volp[ip], mp[ip], my_rho[ip], rho_error[ip], up[ip][0], up[ip][1]);
    }

    delete[] my_rho;
    
    FOR_IT {
    if ((triVec[it].ve_bit[0] == 0) && (triVec[it].ve_bit[1] == 0) && (triVec[it].ve_bit[2] == 0))
      fprintf(fp,"%d %d %d\n",triVec[it].ve_ip[0]+1,triVec[it].ve_ip[1]+1,triVec[it].ve_ip[2]+1);
    }
    
    fclose(fp);
    
  }
  
};
