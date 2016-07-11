#include "cJSON.h"
#define EPSZ 8.85418781762039e-12
#define MUZ 1.25663706143592e-06
#define CC 299792458
#define PI 3.14159265358979
class Grid
{
public:
	double ez, ex, ey, hz, hx, hy;
	unsigned char medium;
public:
	Grid()
	{
		ez = 0;
		ex = 0;
		ey = 0;
		hz = 0;
		hx = 0;
		hy = 0;
		medium = 0;
	}
};

class PmlGrid : public Grid
{
public:
	double psi_Ezy, psi_Eyz, psi_Exz, psi_Ezx, psi_Eyx, psi_Exy;
	double psi_Hzy, psi_Hyz, psi_Hxz, psi_Hzx, psi_Hyx, psi_Hxy;
public:
	PmlGrid()
	{
		psi_Ezy = 0;
		psi_Eyz = 0;
		psi_Exz = 0;
		psi_Ezx = 0;
		psi_Eyx = 0;
		psi_Exy = 0;

		psi_Hzy = 0;
		psi_Hyz = 0;
		psi_Hxz = 0;
		psi_Hzx = 0;
		psi_Hyx = 0;
		psi_Hxy = 0;
	}
};

class FdtdManager
{
private:
	int m_Nx, m_Ny, m_Nz;
	int m_Nxtot, m_Nytot, m_Nztot;
	int m_Npml;
	int m_Nmedium;
	int m_T;
	double m_dx, m_dy, m_dz, m_dt;
	double m_Epsr[256], m_Sigma[256];
	double *m_source;
	Grid**** grids;
	int source_x;
	int source_y;
	int source_z;
	cJSON * root;
private:
	void cross(double * x, double * y, double * r);
	void vec_sub(double * x, double * y, double * r);
	int is_between(double x, double a, double b);
	double dot(double * x, double * y);
	int testCone(double x0, double y0, double z0,
		double x1, double y1, double z1,
		double r, double x, double y, double z);
	float testCylinder(double x0, double y0, double z0,
		double x1, double y1, double z1,
		double l, double r,
		double x, double y, double z);
	int testCuboid(double x0, double y0, double z0,
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double x3, double y3, double z3,
		double x, double y, double z);
	void copy_x(Grid ****, int source, int destination, int y_min, int y_max, int z_min, int z_max);
	void copy_y(Grid ****, int source, int destination, int x_min, int x_max, int z_min, int z_max);
	void copy_z(Grid ****, int source, int destination, int x_min, int x_max, int y_min, int y_max);
	char* load_file(char const* path);

	void xx(int i, int j, int k, double& eps, double& sig);
	void yy(int i, int j, int k, double& eps, double& sig);
	void zz(int i, int j, int k, double& eps, double& sig);
	double CAx(int i, int j, int k);
	double CAy(int i, int j, int k);
	double CAz(int i, int j, int k);
	double CBx(int i, int j, int k);
	double CBy(int i, int j, int k);
	double CBz(int i, int j, int k);
	inline double CPx(int i, int j, int k) { return 1.0; }
	inline double CPy(int i, int j, int k) { return 1.0; }
	inline double CPz(int i, int j, int k) { return 1.0; }
	inline double CQx(int i, int j, int k) { return m_dt / MUZ; }
	inline double CQy(int i, int j, int k) { return m_dt / MUZ; }
	inline double CQz(int i, int j, int k) { return m_dt / MUZ; }
	void calcExAC(int i, int j, int k, double& A, double& C);
	void calcEyAC(int i, int j, int k, double& A, double& C);
	void calcEzAC(int i, int j, int k, double& A, double& C);
	void calcHxAC(int i, int j, int k, double& A, double& C);
	void calcHyAC(int i, int j, int k, double& A, double& C);
	void calcHzAC(int i, int j, int k, double& A, double& C);
	void updateEField();
	void updateHField();
	void updatePmlEField();
	void updatePmlHField();
public:
	void init(char * filename);
	void loadJson(char* filename);
	void beginFdtd();
};