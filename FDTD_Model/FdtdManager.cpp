#include <cstdio>
#include <cmath>
#include <cstring>
#include <omp.h>
#include "cJSON.h"
#include "FdtdManager.h"
#define testBall(x0, y0, z0, r, x, y, z) ((pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)) < pow(r, 2))
#define get_Item cJSON_GetObjectItem
#define get_Size cJSON_GetArraySize
#define get_ArrItem cJSON_GetArrayItem

///////////////public method/////////////////////
/////////////////////////////////////////////////
void FdtdManager::init(char* filename)
{
	// FILE * output = fopen(filename, "w+");
	char* charList = load_file(filename);
	root = cJSON_Parse(charList);
	// delete[] charList;
	charList = NULL;
	double x_len = get_Item(root, "x_len")->valuedouble;
	double y_len = get_Item(root, "y_len")->valuedouble;
	double z_len = get_Item(root, "z_len")->valuedouble;
	double dx = get_Item(root, "dx")->valuedouble;
	double dy = get_Item(root, "dy")->valuedouble;
	double dz = get_Item(root, "dz")->valuedouble;

	double source_x_double = get_Item(root, "source_x")->valuedouble;
	double source_y_double = get_Item(root, "source_y")->valuedouble;
	double source_z_double = get_Item(root, "source_z")->valuedouble;

//	cJSON * obj = get_Item(root, "snapshot"
	int default_material_type = get_Item(root, "default_material_type")->valueint;
	///////////////////////////////////////////////
	m_Npml = get_Item(root, "pml_layers")->valueint + 1;
	////////////////////////////////////////////////
	//grids
	m_Nx = (int)x_len / dx;
	m_Ny = (int)y_len / dy;
	m_Nz = (int)z_len / dz;

	source_x = round(source_x_double / dx) + m_Npml;
	source_y = round(source_y_double / dy) + m_Npml;
	source_z = round(source_z_double / dz) + m_Npml;

	m_dx = dx;
	m_dy = dy;
	m_dz = dz;

	m_Nxtot = 2 * m_Npml + m_Nx;
	m_Nytot = 2 * m_Npml + m_Ny;
	m_Nztot = 2 * m_Npml + m_Nz;

	m_dt = 0.99 / CC / sqrt(1.0 / (m_dx*m_dx) + 1.0 / (m_dy*m_dy) + 1.0 / (m_dz*m_dz));
	printf("time step: %ef\n", m_dt);
	double f = get_Item(root, "f")->valuedouble;
	double timewindow = get_Item(root, "timewindow")->valuedouble;
	m_T = (int)(timewindow / m_dt) + 1;
	m_source = new double[m_T];
	double ks = 2 * PI*PI*f*f;
	double X = 1 / f;

	for (int i = 0; i < m_T; ++i) {
		double t = m_dt*i;
		m_source[i] = 2 * ks*sqrt(exp(1 / 2 / ks))*exp(-ks*(t - X)*(t - X))*(t - X);
	}
	
	///////源归一化//////
	double min = m_source[0];
	double max = m_source[0];
	for (int i = 0; i < m_T; ++i)
	{
		if (m_source[i] > max)
		{
			max = m_source[i];
		}
		if (m_source[i] < min)
		{
			min = m_source[i];
		}
	}
	for (int i = 0; i < m_T; ++i)
	{
		m_source[i] = 1000000 * (m_source[i] / (max - min));
	}
	////////////////////

	int xx = m_Nx + 2 * m_Npml;
	int yy = m_Ny + 2 * m_Npml;
	int zz = m_Nz + 2 * m_Npml;

	/////////set property//////

	///////////////////////////

	printf("%d\n%d\n%d\n", xx, yy, zz);

	////////////init as one dim array////////////
	//Node * Grid = (Node *) calloc(xx * yy * zz, sizeof(Node));
	///////////////multi dynamic array////////////
	//Node *** Grid = calloc(xx, sizeof(Node **));
	grids = new Grid***[xx];
	for (int i = 0; i < xx; i++)
	{
		grids[i] = new Grid**[yy];
		for (int j = 0; j < yy; j++)
		{
			grids[i][j] = new Grid*[zz];
		}
	}

	int x_max = m_Nx + m_Npml;
	int y_max = m_Ny + m_Npml;
	int z_max = m_Nz + m_Npml;



	for (int i = 0; i < xx; i++)
	{
		for (int j = 0; j < yy; j++)
		{
			for (int k = 0; k < zz; k++)
			{
				if (
					(i < x_max) && (i >= m_Npml) &&
					(j < y_max) && (j >= m_Npml) &&
					(k < z_max) && (k >= m_Npml)
					)
				{
					grids[i][j][k] = new Grid;
				}
				else
				{
					grids[i][j][k] = new PmlGrid;
				}
			}
		}
	}

	printf("memory allocated!\n");

	cJSON * model = get_Item(root, "model");
	//  #pragma omp parallel for 
	for (int i = 0; i < m_Nx; i++)
	{
		for (int j = 0; j < m_Ny; j++)
		{
			for (int k = 0; k < m_Nz; k++)
			{
				grids[i + m_Npml][j + m_Npml][k + m_Npml]->medium = (char)get_Item(root, "default_material_type")->valueint;
				for (int m = 0; m < get_Size(model); m++)
				{
					cJSON * obj = get_ArrItem(model, m);
					char * type = get_Item(obj, "type")->valuestring;
					char material_type = (char)get_Item(obj, "material_type")->valueint;
					cJSON * args_json = get_Item(obj, "args");
					////////////////////cuboid//////////////////////////
					if (strcmp(type, "cuboid") == 0)
					{
						//            printf("cuboid\n"); 
						double args[12];
						for (int index = 0; index < 12; index++)
						{
							args[index] = get_ArrItem(args_json, index)->valuedouble;
							//              printf("%f\n", args[index]);
						}

						if (testCuboid(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10], args[11], dx * (i + 0.5), dy * (j + 0.5), dz * (k + 0.5)))
						{
							//              printf("cuboid\n");
							grids[i + m_Npml][j + m_Npml][k + m_Npml]->medium = material_type;
						}
					}
					//////////////////ball//////////////////////////
					if (strcmp(type, "ball") == 0)
					{
						//            printf("ball\n");
						double args[4];
						for (int index = 0; index < 4; index++)
						{
							args[index] = get_ArrItem(args_json, index)->valuedouble;
						}

						if (testBall(args[0], args[1], args[2], args[3], dx * (i + 0.5), dy * (j + 0.5), dz * (k + 0.5)))
						{
							//              printf("ball\n");
							grids[i + m_Npml][j + m_Npml][k + m_Npml]->medium = material_type;
						}
					}
					/////////////////////cylinder///////////////////
					if (strcmp(type, "cylinder") == 0)
					{
						//            printf("cylinder\n");
						double args[8];
						for (int index = 0; index < 8; index++)
						{
							args[index] = get_ArrItem(args_json, index)->valuedouble;
						}

						if (testCylinder(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], dx * (i + 0.5), dy * (j + 0.5), dz * (k + 0.5)))
						{
							//             printf("cylinder\n");
							grids[i + m_Npml][j + m_Npml][k + m_Npml]->medium = material_type;
							//              printf("%d\n", Grid[Grid_Index].material_type);
						}
					}
					////////////////////cone////////////////////////
					if (strcmp(type, "cone") == 0)
					{
						double args[7];
						for (int index = 0; index < 7; index++)
						{
							args[index] = get_ArrItem(args_json, index)->valuedouble;
						}

						if (testCone(args[0], args[1], args[2], args[3], args[4], args[5], args[6], dx * (i + 0.5), dy * (j + 0.5), dz * (k + 0.5)))
						{
							grids[i + m_Npml][j + m_Npml][k + m_Npml]->medium = material_type;
						}
					}

				}
			}
		}
	}

	////////set pml media////////////
	//up
	for (int i = m_Npml + m_Nz; i < m_Npml * 2 + m_Nz; i++)
		copy_z(grids,
			m_Npml + m_Nz - 1,
			i,
			m_Npml,
			m_Npml + m_Nx,
			m_Npml,
			m_Npml + m_Ny);
	//down
	for (int i = 0; i < m_Npml; i++)
		copy_z(grids, m_Npml, i, m_Npml, m_Npml + m_Nx, m_Npml, m_Npml + m_Ny);
	//left
	for (int i = 0; i < m_Npml; i++)
		copy_y(grids, m_Npml, i, m_Npml, m_Npml + m_Nx, m_Npml, m_Npml + m_Nz);
	//right
	for (int i = m_Npml + m_Ny; i < m_Npml * 2 + m_Ny; i++)
		copy_y(grids, m_Npml + m_Ny - 1, i, m_Npml, m_Npml + m_Nx, m_Npml, m_Npml + m_Nz);
	//back
	for (int i = 0; i < m_Npml; i++)
		copy_x(grids, m_Npml, i, m_Npml, m_Npml + m_Ny, m_Npml, m_Npml + m_Nz);
	//      copy_x(grids, m_Npml, i, m_Npml, m_Npml + m_Ny, m_Npml, m_Npml + m_Nz);
	//front
	for (int i = m_Npml + m_Nx; i < m_Npml * 2 + m_Nx; i++)
		copy_x(grids, m_Npml + m_Nx - 1, i, m_Npml, m_Npml + m_Ny, m_Npml, m_Npml + m_Nz);
	//edge
	//4 edges down
	for (int i = 0; i < m_Npml; i++)
	{
		copy_z(grids, m_Npml, i,
			m_Npml, m_Npml + m_Nx,
			0, m_Npml
			);
		copy_z(grids, m_Npml, i,
			m_Npml, m_Npml + m_Nx,
			m_Npml + m_Ny, m_Npml * 2 + m_Ny
			);
		copy_z(grids, m_Npml, i,
			0, m_Npml,
			m_Npml, m_Npml + m_Ny
			);
		copy_z(grids, m_Npml, i,
			m_Npml + m_Nx, m_Npml * 2 + m_Nx,
			m_Npml, m_Npml + m_Ny
			);
	}
	//4 edges up
	for (int i = m_Npml + m_Nz; i < m_Npml * 2 + m_Nz; i++)
	{
		copy_z(grids, m_Npml + m_Nz - 1, i,
			m_Npml, m_Npml + m_Nx,
			0, m_Npml
			);
		copy_z(grids, m_Npml + m_Nz - 1, i,
			m_Npml, m_Npml + m_Nx,
			m_Npml + m_Ny, m_Npml * 2 + m_Ny
			);
		copy_z(grids, m_Npml + m_Nz - 1, i,
			0, m_Npml,
			m_Npml, m_Npml + m_Ny
			);
		copy_z(grids, m_Npml + m_Nz - 1, i,
			m_Npml + m_Nx, m_Npml * 2 + m_Nx,
			m_Npml, m_Npml + m_Ny
			);
	}
	//2 edges back
	for (int i = 0; i < m_Npml; i++)
	{
		copy_x(grids, m_Npml, i,
			0, m_Npml,
			m_Npml, m_Npml + m_Nz
			);
		copy_x(grids, m_Npml, i,
			m_Npml + m_Ny, 2 * m_Npml + m_Ny,
			m_Npml, m_Npml + m_Nz
			);
	}
	//2 edges front
	for (int i = m_Npml + m_Nx; i < m_Npml * 2 + m_Nx; i++)
	{
		copy_x(grids, m_Npml + m_Nx - 1, i,
			0, m_Npml,
			m_Npml, m_Npml + m_Nz
			);
		copy_x(grids, m_Npml + m_Nx - 1, i,
			m_Npml + m_Ny, 2 * m_Npml + m_Ny,
			m_Npml, m_Npml + m_Nz
			);

	}
	//corners
	//4 corners down
	for (int i = 0; i < m_Npml; i++)
	{
		copy_z(grids, m_Npml, i,
			0, m_Npml,
			0, m_Npml
			);
		copy_z(grids, m_Npml, i,
			0, m_Npml,
			m_Npml + m_Ny, m_Npml * 2 + m_Ny
			);
		copy_z(grids, m_Npml, i,
			m_Npml + m_Nx, m_Npml * 2 + m_Nx,
			0, m_Npml
			);
		copy_z(grids, m_Npml, i,
			m_Npml + m_Nx, m_Npml * 2 + m_Nx,
			m_Npml + m_Ny, m_Npml * 2 + m_Ny
			);
	}
	//4 corners up
	for (int i = m_Npml + m_Nz; i < m_Npml * 2 + m_Nz; i++)
	{
		copy_z(grids, m_Npml + m_Nz - 1, i,
			0, m_Npml,
			0, m_Npml
			);
		copy_z(grids, m_Npml + m_Nz - 1, i,
			0, m_Npml,
			m_Npml + m_Ny, m_Npml * 2 + m_Ny
			);
		copy_z(grids, m_Npml + m_Nz - 1, i,
			m_Npml + m_Nx, m_Npml * 2 + m_Nx,
			0, m_Npml
			);
		copy_z(grids, m_Npml + m_Nz - 1, i,
			m_Npml + m_Nx, m_Npml * 2 + m_Nx,
			m_Npml + m_Ny, m_Npml * 2 + m_Ny
			);
	}

	//set medium para
	cJSON * material_para = get_Item(root, "material_para");
	m_Nmedium = get_Size(material_para);
	for (int m = 0; m < get_Size(material_para); m++)
	{
		cJSON * material = get_ArrItem(material_para, m);
		m_Epsr[m] = (float)get_Item(material, "eps")->valuedouble;
		m_Sigma[m] = (float)get_Item(material, "sig")->valuedouble;
		;
	}
	//set npml
	m_Npml = m_Npml - 1;
	//check results
	FILE * result = fopen("result.txt", "w+");
	for (int i = 0; i < xx; i++)
	{
		for (int j = 0; j < yy; j++)
		{
			for (int k = 0; k < zz; k++)
			{
				fprintf(result, "%d\n", grids[i][j][k]->medium);
			}
		}
	}

	//      copy_x(grids, m_Npml + m_Nx - 1， i, m_Npml + m_Ny, m_Npml, m_Npml + m_Nz);
	//write to file
	// if(Grid[i][j][k].material_type)
	//  {
	//  fprintf(output, "%d\n", Grid[i][j][k].material_type);
	// }
	// else
	// {
	//   fprintf(output, "%d\n", 0);
	// }
}

// void FdtdManager::buildModel()
// {
// 	printf("Building Model\n");
// }
// ////////////////private method///////////////////
/////////////////////////////////////////////////

void FdtdManager::cross(double * x, double * y, double * r)
{
	r[0] = -x[2] * y[1] + x[1] * y[2];
	r[1] = x[2] * y[0] - x[0] * y[2];
	r[2] = -x[1] * y[0] + x[0] * y[1];
}

void FdtdManager::vec_sub(double * x, double * y, double * r)
{
	r[0] = x[0] - y[0];
	r[1] = x[1] - y[1];
	r[2] = x[2] - y[2];
}

int FdtdManager::is_between(double x, double a, double b)
{
	if (x > a && x < b) return 1;
	if (x < a && x > b) return 1;
	return 0;
}

double FdtdManager::dot(double * x, double * y)
{
	return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

int FdtdManager::testCone(double x0, double y0, double z0, double x1, double y1, double z1, double r, double x, double y, double z)
{
	double apex[3] = { x0, y0, z0 }; //顶点
	double center[3] = { x1, y1, z1 }; //圆心
	double point[3] = { x, y, z };  //测试点
	double dir[3];
	vec_sub(apex, center, dir);
	double dir_len = sqrt(pow(dir[0], 2) + pow(dir[1], 2) + pow(dir[2], 2));
	dir[0] /= dir_len;
	dir[1] /= dir_len;
	dir[2] /= dir_len;
	double apex2point[3];
	vec_sub(apex, point, apex2point);
	double cone_distance = dot(apex2point, dir);
	double axis_len = sqrt(
		pow(x0 - x1, 2) +
		pow(y0 - y1, 2) +
		pow(z0 - z1, 2)
		);
	if (cone_distance < 0 || cone_distance > axis_len)
	{
		return 0;
	}

	double radium_intersect = (r * cone_distance) / axis_len; //截面圆半径
	double apex2point_len_squared =
		pow(apex2point[0], 2) +
		pow(apex2point[1], 2) +
		pow(apex2point[2], 2);

	double point2axis_distance = sqrt(apex2point_len_squared - pow(cone_distance, 2)); //点到轴距离
	if (point2axis_distance > radium_intersect)
	{
		return 0;
	}
	return 1;
}
int FdtdManager::testCuboid(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double x, double y, double z)
{
	double p1[3];
	p1[0] = x0;
	p1[1] = y0;
	p1[2] = z0;
	double p2[3];
	p2[0] = x1;
	p2[1] = y1;
	p2[2] = z1;
	double p3[3];
	p3[0] = x2;
	p3[1] = y2;
	p3[2] = z2;
	double p4[3];
	p4[0] = x3;
	p4[1] = y3;
	p4[2] = z3;
	double px[3];
	px[0] = x;
	px[1] = y;
	px[2] = z;

	double u1[3], u2[3], u3[3], u[3], v[3], w[3];
	vec_sub(p1, p2, u1);
	vec_sub(p1, p3, u2);
	vec_sub(p1, p4, u3);
	cross(u1, u2, u);
	if (!is_between(dot(px, u), dot(p1, u), dot(p4, u)))
	{
		return 0;
	}
	cross(u1, u3, v);
	if (!is_between(dot(px, v), dot(p1, v), dot(p3, v)))
	{
		return 0;
	}
	cross(u2, u3, w);
	if (!is_between(dot(px, w), dot(p1, w), dot(p2, w)))
	{
		return 0;
	}
	return 1;
}


float FdtdManager::testCylinder(double x0, double y0, double z0, double x1, double y1, double z1, double l, double r, double x, double y, double z)
{
	double dx, dy, dz;
	double dot, dsq;
	double pdx, pdy, pdz;
	dx = x1 - x0;
	dy = y1 - y0;
	dz = z1 - z0;
	pdx = x - x0;
	pdy = y - y0;
	pdz = z - z0;
	dot = pdx * dx + pdy * dy + pdz * dz;
	if (dot < 0.0f || dot > l)
	{
		return(0.0);
	}
	else
	{
		dsq = (pdx*pdx + pdy*pdy + pdz*pdz) - dot*dot / l;
		if (dsq > r)
		{
			return(0.0);
		}
		else
		{
			return(1.0);
		}
	}
}


//to load .json file as char[]
char* FdtdManager::load_file(char const* path)
{
	char* buffer = 0;
	long length;
	FILE * f = fopen(path, "rb"); //was "rb"

	if (f)
	{
		fseek(f, 0, SEEK_END);
		length = ftell(f);
		fseek(f, 0, SEEK_SET);
		//buffer = (char*)malloc((length + 1)*sizeof(char));
		buffer = new char[length + 1];
		if (buffer)
		{
			fread(buffer, sizeof(char), length, f);
		}
		fclose(f);
		buffer[length + 1] = '\0';
		;
	}
	else
	{
		printf("unable to open %s\n", path);
	}
	return buffer;
}

		// for (int i = 0; i < length; i++) {
		//     printf("buffer[%d] == %c\n", i, buffer[i]);
		// }
		//printf("buffer = %s\n", buffer);



		void FdtdManager::copy_x(Grid **** grids, int source, int destination, int y_min, int y_max, int z_min, int z_max)
	{
		for (int j = y_min; j < y_max; j++)
		{
			for (int k = z_min; k < z_max; k++)
			{
				grids[destination][j][k]->medium = grids[source][j][k]->medium;
			}
		}
	}

	void FdtdManager::copy_y(Grid **** grids, int source, int destination, int x_min, int x_max, int z_min, int z_max)
	{
		for (int i = x_min; i < x_max; i++)
		{
			for (int k = z_min; k < z_max; k++)
			{
				grids[i][destination][k]->medium = grids[i][source][k]->medium;
			}
		}
	}

	void FdtdManager::copy_z(Grid **** grids, int source, int destination, int x_min, int x_max, int y_min, int y_max)
	{
		for (int i = x_min; i < x_max; i++)
		{
			for (int j = y_min; j < y_max; j++)
			{
				grids[i][j][destination]->medium = grids[i][j][source]->medium;
			}
		}
	}


	void FdtdManager::loadJson(char* filename)
	{
		//test();
		m_Nx = 80;
		m_Ny = 80;
		m_Nz = 80;
		m_Npml = 20;
		m_Nxtot = 2 * m_Npml + m_Nx + 2;
		m_Nytot = 2 * m_Npml + m_Ny + 2;
		m_Nztot = 2 * m_Npml + m_Nz + 2;
		m_Nmedium = 1;
		m_Epsr[0] = 1;
		m_Sigma[0] = 0;
		m_dx = 0.1;
		m_dy = 0.1;
		m_dz = 0.1;
		m_dt = 0.99 / CC / sqrt(1.0 / (m_dx*m_dx) + 1.0 / (m_dy*m_dy) + 1.0 / (m_dz*m_dz));
		printf("dt :%ef\n", m_dt);
		double f = 150e6;
		double timewindow = 60e-9;
		m_T = (int)(timewindow / m_dt) + 1;
		m_source = new double[m_T];
		double ks = 2 * PI*PI*f*f;
		double X = 1 / f;
		for (int i = 0; i < m_T; ++i) {
			double t = m_dt*i;
			m_source[i] = 2 * ks*sqrt(exp(1 / 2 / ks))*exp(-ks*(t - X)*(t - X))*(t - X);
		}
		grids = new Grid***[m_Nxtot];
		for (int i = 0; i < m_Nxtot; ++i) {
			grids[i] = new Grid**[m_Nytot];
			for (int j = 0; j < m_Nytot; ++j) {
				grids[i][j] = new Grid*[m_Nztot];
				for (int k = 0; k < m_Nztot; ++k) {
					if (i <= m_Npml || j <= m_Npml || k <= m_Npml || i > m_Npml + m_Nx || j > m_Npml + m_Ny || k > m_Npml + m_Nz) {
						grids[i][j][k] = new PmlGrid();
						//	printf("%d ", grids[i][j][k]->medium);
					}
					else {
						grids[i][j][k] = new Grid();
						//	printf("%d ", grids[i][j][k]->medium);
					}
				}
			}
		}
		printf("loading %s\n", filename);
		printf("Total Step is : %d\n", m_T);
	}


	void FdtdManager::xx(int i, int j, int k, double& eps, double& sig)
	{
		Grid* grid = grids[i][j][k];
		if (k == 1 || j == 1) {
			eps = EPSZ*m_Epsr[grid->medium];
			sig = m_Sigma[grid->medium];
		}
		else {
			eps = EPSZ*(m_Epsr[grid->medium] + m_Epsr[grids[i][j - 1][k]->medium]
				+ m_Epsr[grids[i][j][k - 1]->medium] + m_Epsr[grids[i][j - 1][k - 1]->medium]) / 4.0;
			sig = (m_Sigma[grid->medium] + m_Sigma[grids[i][j - 1][k]->medium]
				+ m_Sigma[grids[i][j][k - 1]->medium] + m_Sigma[grids[i][j - 1][k - 1]->medium]) / 4.0;
		}

	}

	void FdtdManager::yy(int i, int j, int k, double& eps, double& sig)
	{
		Grid* grid = grids[i][j][k];
		if (k == 1 || i == 1) {
			eps = EPSZ*m_Epsr[grid->medium];
			sig = m_Sigma[grid->medium];
		}
		else {
			eps = EPSZ*(m_Epsr[grid->medium] + m_Epsr[grids[i - 1][j][k]->medium]
				+ m_Epsr[grids[i][j][k - 1]->medium] + m_Epsr[grids[i - 1][j][k - 1]->medium]) / 4.0;
			sig = (m_Sigma[grid->medium] + m_Sigma[grids[i - 1][j][k]->medium]
				+ m_Sigma[grids[i][j][k - 1]->medium] + m_Sigma[grids[i - 1][j][k - 1]->medium]) / 4.0;
		}
	}

	void FdtdManager::zz(int i, int j, int k, double& eps, double& sig)
	{
		Grid* grid = grids[i][j][k];
		if (i == 1 || j == 1) {
			eps = EPSZ*m_Epsr[grid->medium];
			sig = m_Sigma[grid->medium];
		}
		else {
			eps = EPSZ*(m_Epsr[grid->medium] + m_Epsr[grids[i - 1][j][k]->medium]
				+ m_Epsr[grids[i][j - 1][k]->medium] + m_Epsr[grids[i - 1][j - 1][k]->medium]) / 4.0;
			sig = (m_Sigma[grid->medium] + m_Sigma[grids[i - 1][j][k]->medium]
				+ m_Sigma[grids[i][j - 1][k]->medium] + m_Sigma[grids[i - 1][j - 1][k]->medium]) / 4.0;
		}
	}

	double FdtdManager::CAx(int i, int j, int k)
	{
		double sig, eps;
		xx(i, j, k, eps, sig);
		double temp = sig*m_dt / eps / 2.0;
		//printf("sig = %ef eps = %ef CAx:%ef ", sig, eps, (1.0-temp)/(1.0+temp));
		return (1.0 - temp) / (1.0 + temp);
	}

	double FdtdManager::CBx(int i, int j, int k)
	{
		double sig, eps;
		xx(i, j, k, eps, sig);
		double temp = m_dt / eps;
		//printf("sig = %ef eps = %ef CBx:%ef ", sig, eps, temp);
		return temp / (1.0 + temp*sig / 2.0);
	}

	double FdtdManager::CAy(int i, int j, int k)
	{
		double sig, eps;
		yy(i, j, k, eps, sig);
		double temp = sig*m_dt / eps / 2.0;
		//printf("sig = %ef eps = %ef CAy:%ef ", sig, eps, (1.0-temp)/(1.0+temp));
		return (1.0 - temp) / (1.0 + temp);
	}

	double FdtdManager::CBy(int i, int j, int k)
	{
		double sig, eps;
		yy(i, j, k, eps, sig);
		double temp = m_dt / eps;
		//printf("sig = %ef eps = %ef CBy:%ef ", sig, eps, temp);
		return temp / (1.0 + temp*sig / 2.0);
	}

	double FdtdManager::CAz(int i, int j, int k)
	{
		double sig, eps;
		zz(i, j, k, eps, sig);
		double temp = sig*m_dt / eps / 2.0;
		//printf("sig = %ef epsr = %ef CAz:%ef ", sig, eps, (1.0-temp)/(1.0+temp));
		return (1.0 - temp) / (1.0 + temp);
	}

	double FdtdManager::CBz(int i, int j, int k)
	{
		double sig, eps;
		zz(i, j, k, eps, sig);
		double temp = m_dt / eps;
		//printf("sig = %ef eps = %ef CBz:%ef ", sig, eps,temp);
		return temp / (1.0 + temp*sig / 2.0);
	}

	void FdtdManager::updateEField()
	{
		//Grid**** grids = this->grids;
#pragma omp parallel for
		for (int i = m_Npml + 1; i <= m_Npml + m_Nx; ++i) {
			for (int j = m_Npml + 1; j <= m_Npml + m_Ny; ++j) {
				for (int k = m_Npml + 1; k <= m_Npml + m_Nz; ++k) {
					grids[i][j][k]->ex = CAx(i, j, k)*grids[i][j][k]->ex + CBx(i, j, k)*((grids[i][j][k]->hz - grids[i][j - 1][k]->hz) / m_dy - (grids[i][j][k]->hy - grids[i][j][k - 1]->hy) / m_dz);
					grids[i][j][k]->ey = CAy(i, j, k)*grids[i][j][k]->ey + CBy(i, j, k)*((grids[i][j][k]->hx - grids[i][j][k - 1]->hx) / m_dz - (grids[i][j][k]->hz - grids[i - 1][j][k]->hz) / m_dx);
					grids[i][j][k]->ez = CAz(i, j, k)*grids[i][j][k]->ez + CBz(i, j, k)*((grids[i][j][k]->hy - grids[i - 1][j][k]->hy) / m_dx - (grids[i][j][k]->hx - grids[i][j - 1][k]->hx) / m_dy);
				}
			}
		}
	}
	void FdtdManager::updateHField()
	{
		//Grid**** grids = this->grids;
#pragma omp parallel for
		for (int i = m_Npml + 1; i <= m_Npml + m_Nx; ++i) {
			for (int j = m_Npml + 1; j <= m_Npml + m_Ny; ++j) {
				for (int k = m_Npml + 1; k <= m_Npml + m_Nz; ++k) {
					grids[i][j][k]->hx = CPx(i, j, k)*grids[i][j][k]->hx - CQx(i, j, k)*((grids[i][j + 1][k]->ez - grids[i][j][k]->ez) / m_dy - (grids[i][j][k + 1]->ey - grids[i][j][k]->ey) / m_dz);
					grids[i][j][k]->hy = CPy(i, j, k)*grids[i][j][k]->hy - CQy(i, j, k)*((grids[i][j][k + 1]->ex - grids[i][j][k]->ex) / m_dz - (grids[i + 1][j][k]->ez - grids[i][j][k]->ez) / m_dx);
					grids[i][j][k]->hz = CPz(i, j, k)*grids[i][j][k]->hz - CQz(i, j, k)*((grids[i + 1][j][k]->ey - grids[i][j][k]->ey) / m_dx - (grids[i][j + 1][k]->ex - grids[i][j][k]->ex) / m_dy);
				}
			}
		}
	}

	void FdtdManager::calcExAC(int i, int j, int k, double& A, double& C)
	{
		double sigmax = 1.0 / 30 / PI / m_dx / sqrt(m_Epsr[grids[i][j][k]->medium]);
		double sig = 0;
		if (i <= m_Npml) {
			sig = sigmax*pow((double)(m_Npml + 1 - i) / m_Npml, 4);
		}
		else if (i > m_Npml + m_Nx) {
			sig = sigmax*pow((double)(i - m_Npml - m_Nx - 1) / m_Npml, 4);
		}
		double temp = sig / EPSZ;
		A = exp(-temp*m_dt);
		C = A - 1.0;
	}
	void FdtdManager::calcEyAC(int i, int j, int k, double& A, double& C)
	{
		double sigmax = 1.0 / 30 / PI / m_dy / sqrt(m_Epsr[grids[i][j][k]->medium]);
		double sig = 0;
		if (j <= m_Npml) {
			sig = sigmax*pow((double)(m_Npml + 1 - j) / m_Npml, 4);
		}
		else if (j > m_Npml + m_Ny) {
			sig = sigmax*pow((double)(j - m_Npml - m_Ny - 1) / m_Npml, 4);
		}
		double temp = sig / EPSZ;
		A = exp(-temp*m_dt);
		C = A - 1.0;
	}
	void FdtdManager::calcEzAC(int i, int j, int k, double& A, double& C)
	{
		double sigmax = 1.0 / 30 / PI / m_dz / sqrt(m_Epsr[grids[i][j][k]->medium]);
		double sig = 0;
		if (k <= m_Npml) {
			sig = sigmax*pow((double)(m_Npml + 1 - k) / m_Npml, 4);
		}
		else if (k > m_Npml + m_Nz) {
			sig = sigmax*pow((double)(k - m_Npml - m_Nz - 1) / m_Npml, 4);
		}
		double temp = sig / EPSZ;
		A = exp(-temp*m_dt);
		C = A - 1.0;
	}
	void FdtdManager::calcHxAC(int i, int j, int k, double& A, double& C)
	{
		double sigmax = 1.0 / 30 / PI / m_dx / sqrt(m_Epsr[grids[i][j][k]->medium]);
		double sig = 0;
		if (i <= m_Npml) {
			sig = sigmax*pow((double)(m_Npml + 0.5 - i) / m_Npml, 4);
		}
		else if (i > m_Npml + m_Nx) {
			sig = sigmax*pow((double)(i - m_Npml - m_Nx - 0.5) / m_Npml, 4);
		}
		double temp = sig / EPSZ;
		A = exp(-temp*m_dt);
		C = A - 1.0;
	}
	void FdtdManager::calcHyAC(int i, int j, int k, double& A, double& C)
	{
		double sigmax = 1.0 / 30 / PI / m_dy / sqrt(m_Epsr[grids[i][j][k]->medium]);
		double sig = 0;
		if (j <= m_Npml) {
			sig = sigmax*pow((double)(m_Npml + 0.5 - j) / m_Npml, 4);
		}
		else if (j > m_Npml + m_Ny) {
			sig = sigmax*pow((double)(j - m_Npml - m_Ny - 0.5) / m_Npml, 4);
		}
		double temp = sig / EPSZ;
		A = exp(-temp*m_dt);
		C = A - 1.0;
	}
	void FdtdManager::calcHzAC(int i, int j, int k, double& A, double& C)
	{
		double sigmax = 1.0 / 30 / PI / m_dz / sqrt(m_Epsr[grids[i][j][k]->medium]);
		double sig = 0;
		if (k <= m_Npml) {
			sig = sigmax*pow((double)(m_Npml + 0.5 - k) / m_Npml, 4);
		}
		else if (k > m_Npml + m_Nz) {
			sig = sigmax*pow((double)(k - m_Npml - m_Nz - 0.5) / m_Npml, 4);
		}
		double temp = sig / EPSZ;
		A = exp(-temp*m_dt);
		C = A - 1.0;
	}
	void FdtdManager::updatePmlEField()
	{
#pragma omp parallel for
		for (int i = 1; i <= 2 * m_Npml + m_Nx; ++i) {
			for (int j = 1; j <= 2 * m_Npml + m_Ny; ++j) {
				for (int k = 1; k <= 2 * m_Npml + m_Nz; ++k) {
					if (i > m_Npml && i <= m_Npml + m_Nx && j > m_Npml && j <= m_Npml + m_Ny && k > m_Npml && k <= m_Npml + m_Nz) {
						continue;
					}
					double Ax, Ay, Az, Cx, Cy, Cz;
					calcExAC(i, j, k, Ax, Cx);
					calcEyAC(i, j, k, Ay, Cy);
					calcEzAC(i, j, k, Az, Cz);
					((PmlGrid*)grids[i][j][k])->psi_Ezy = Ay*((PmlGrid*)grids[i][j][k])->psi_Ezy + Cy*((grids[i][j][k]->hz - grids[i][j - 1][k]->hz) / m_dy);
					((PmlGrid*)grids[i][j][k])->psi_Eyz = Az*((PmlGrid*)grids[i][j][k])->psi_Eyz + Cz*((grids[i][j][k]->hy - grids[i][j][k - 1]->hy) / m_dz);
					grids[i][j][k]->ex = CAx(i, j, k)*grids[i][j][k]->ex + CBx(i, j, k)*((grids[i][j][k]->hz - grids[i][j - 1][k]->hz) / m_dy - (grids[i][j][k]->hy - grids[i][j][k - 1]->hy) / m_dz + ((PmlGrid*)grids[i][j][k])->psi_Ezy - ((PmlGrid*)grids[i][j][k])->psi_Eyz);
					((PmlGrid*)grids[i][j][k])->psi_Exz = Az*((PmlGrid*)grids[i][j][k])->psi_Exz + Cz*((grids[i][j][k]->hx - grids[i][j][k - 1]->hx) / m_dz);
					((PmlGrid*)grids[i][j][k])->psi_Ezx = Ax*((PmlGrid*)grids[i][j][k])->psi_Ezx + Cx*((grids[i][j][k]->hz - grids[i - 1][j][k]->hz) / m_dx);
					grids[i][j][k]->ey = CAy(i, j, k)*grids[i][j][k]->ey + CBy(i, j, k)*((grids[i][j][k]->hx - grids[i][j][k - 1]->hx) / m_dz - (grids[i][j][k]->hz - grids[i - 1][j][k]->hz) / m_dx + ((PmlGrid*)grids[i][j][k])->psi_Exz - ((PmlGrid*)grids[i][j][k])->psi_Ezx);
					((PmlGrid*)grids[i][j][k])->psi_Eyx = Ax*((PmlGrid*)grids[i][j][k])->psi_Eyx + Cx*((grids[i][j][k]->hy - grids[i - 1][j][k]->hy) / m_dx);
					((PmlGrid*)grids[i][j][k])->psi_Exy = Ay*((PmlGrid*)grids[i][j][k])->psi_Exy + Cy*((grids[i][j][k]->hx - grids[i][j - 1][k]->hx) / m_dy);
					grids[i][j][k]->ez = CAz(i, j, k)*grids[i][j][k]->ez + CBz(i, j, k)*((grids[i][j][k]->hy - grids[i - 1][j][k]->hy) / m_dx - (grids[i][j][k]->hx - grids[i][j - 1][k]->hx) / m_dy + ((PmlGrid*)grids[i][j][k])->psi_Eyx - ((PmlGrid*)grids[i][j][k])->psi_Exy);
				}
			}
		}
	}

	void FdtdManager::updatePmlHField()
	{
#pragma omp parallel for
		for (int i = 1; i <= 2 * m_Npml + m_Nx; ++i) {
			for (int j = 1; j <= 2 * m_Npml + m_Ny; ++j) {
				for (int k = 1; k <= 2 * m_Npml + m_Nz; ++k) {
					if (i > m_Npml && i <= m_Npml + m_Nx && j > m_Npml && j <= m_Npml + m_Ny && k > m_Npml && k <= m_Npml + m_Nz) {
						continue;
					}
					double Ax, Ay, Az, Cx, Cy, Cz;
					calcHxAC(i, j, k, Ax, Cx);
					calcHyAC(i, j, k, Ay, Cy);
					calcHzAC(i, j, k, Az, Cz);
					((PmlGrid*)grids[i][j][k])->psi_Hzy = Ay*((PmlGrid*)grids[i][j][k])->psi_Hzy + Cy*((grids[i][j + 1][k]->ez - grids[i][j][k]->ez) / m_dy);
					((PmlGrid*)grids[i][j][k])->psi_Hyz = Az*((PmlGrid*)grids[i][j][k])->psi_Hyz + Cz*((grids[i][j][k + 1]->ey - grids[i][j][k]->ey) / m_dz);
					grids[i][j][k]->hx = CPx(i, j, k)*grids[i][j][k]->hx - CQx(i, j, k)*((grids[i][j + 1][k]->ez - grids[i][j][k]->ez) / m_dy - (grids[i][j][k + 1]->ey - grids[i][j][k]->ey) / m_dz + ((PmlGrid*)grids[i][j][k])->psi_Hzy - ((PmlGrid*)grids[i][j][k])->psi_Hyz);
					((PmlGrid*)grids[i][j][k])->psi_Hxz = Az*((PmlGrid*)grids[i][j][k])->psi_Hxz + Cz*((grids[i][j][k + 1]->ex - grids[i][j][k]->ex) / m_dz);
					((PmlGrid*)grids[i][j][k])->psi_Hzx = Ax*((PmlGrid*)grids[i][j][k])->psi_Hzx + Cx*((grids[i + 1][j][k]->ez - grids[i][j][k]->ez) / m_dx);
					grids[i][j][k]->hy = CPy(i, j, k)*grids[i][j][k]->hy - CQy(i, j, k)*((grids[i][j][k + 1]->ex - grids[i][j][k]->ex) / m_dz - (grids[i + 1][j][k]->ez - grids[i][j][k]->ez) / m_dx + ((PmlGrid*)grids[i][j][k])->psi_Hxz - ((PmlGrid*)grids[i][j][k])->psi_Hzx);
					((PmlGrid*)grids[i][j][k])->psi_Hyx = Ax*((PmlGrid*)grids[i][j][k])->psi_Hyx + Cx*((grids[i + 1][j][k]->ey - grids[i][j][k]->ey) / m_dx);
					((PmlGrid*)grids[i][j][k])->psi_Hxy = Ay*((PmlGrid*)grids[i][j][k])->psi_Hxy + Cy*((grids[i][j + 1][k]->ex - grids[i][j][k]->ex) / m_dy);
					grids[i][j][k]->hz = CPz(i, j, k)*grids[i][j][k]->hz - CQz(i, j, k)*((grids[i + 1][j][k]->ey - grids[i][j][k]->ey) / m_dx - (grids[i][j + 1][k]->ex - grids[i][j][k]->ex) / m_dy + ((PmlGrid*)grids[i][j][k])->psi_Hyx - ((PmlGrid*)grids[i][j][k])->psi_Hxy);
				}
			}
		}
	}


	void FdtdManager::beginFdtd()
	{
		char filename[64];
		//printf("OMP: %d\n",omp_get_nested());
		//omp_set_nested(1);
		//printf("OMP: %d\n",omp_get_nested());
		printf("Begin Calculate!\n");
		for (int t = 0; t < m_T; ++t) {
			printf("\rStep: %6d", t + 1);
			updateEField();
			updatePmlEField();
			//grids[m_Nxtot / 2][m_Nytot / 2][m_Nztot / 2]->ez = grids[m_Nxtot / 2][m_Nytot / 2][m_Nztot / 2]->ez + m_source[t];
			grids[source_x][source_y][source_z]->ez = grids[source_x][source_y][source_z]->ez + m_source[t];
			updateHField();
			updatePmlHField();

			sprintf(filename, "Exy%d.txt", t + 1);
			if ((t + 1) % (get_Item(root,"snapshot_step")->valueint) == 0 && get_Item(root, "snapshot_xy_flag")->valueint)
			{
				FILE* fid = fopen(filename, "w");
				int z = round(get_Item(root, "snapshot_xy_position")->valuedouble / m_dz) + m_Npml + 1;
				for (int i = 0; i < m_Nxtot; ++i)
				{
					for (int j = 0; j < m_Nytot; ++j)
					{
						fprintf(fid, "%ef ", grids[i][j][z]->ez);
					}
					fprintf(fid, "\n");
				}
				fclose(fid);
			}

			if ((t + 1) % (get_Item(root, "snapshot_step")->valueint) == 0 && get_Item(root, "snapshot_xz_flag")->valueint)
			{
				sprintf(filename, "Exz%d.txt", t + 1);
				FILE* fid = fopen(filename, "w");
				int y = round(get_Item(root, "snapshot_xz_position")->valuedouble / m_dy) + m_Npml + 1;
				for (int i = 0; i < m_Nxtot; ++i)
				{
					for (int k = 0; k < m_Nztot; ++k)
					{
						fprintf(fid, "%ef ", grids[i][y][k]->ez);
					}
					fprintf(fid, "\n");
				}
				fclose(fid);
			}

			//if ((t + 1) % 5 == 0)
			//{
			//	sprintf(filename, "Exz%d.txt", t + 1);
			//	FILE* fid = fopen(filename, "w");
			//	for (int i = 0; i < m_Nxtot; ++i) {
			//		for (int k = 0; k < m_Nztot; ++k) {
			//			fprintf(fid, "%ef ", grids[i][121][k]->ez);
			//		}
			//		fprintf(fid, "\n");
			//	}
			//	fclose(fid);
			//}
		}
		printf("\n");
		/*
		FILE* fid = fopen("Ax.txt", "w");
		for (int i = 0; i < m_Nxtot; ++i) {
		for (int k = 0; k < m_Nztot; ++k) {
		double a, c;
		calcHxAC(i,m_Nytot/2,k, a, c);
		fprintf(fid, "%ef\t", a);
		}
		fprintf(fid, "\n");
		}
		fclose(fid);*/
	}

