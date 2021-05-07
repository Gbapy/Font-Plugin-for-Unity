#pragma once

#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>

using namespace std;

#define		MARGIN	30
#define		PI		3.1415926535897932384626433832795
#define		bwThreshold	30
typedef		unsigned char	UCHAR;
typedef		uint32_t		UINT;
typedef		unsigned long	DWORD;

typedef struct _VECTOR2I_
{
	int val[2];

	_VECTOR2I_() {
		val[0] = 0; val[1] = 0;
	}

	_VECTOR2I_(int i, int j) {
		val[0] = i; val[1] = j;
	}
}Vec2i, *PVec2i;

typedef struct _VECTOR3D_
{
	_VECTOR3D_() {
	
	}

	_VECTOR3D_(double x, double y, double z) {
		val[0] = x; val[1] = y; val[2] = z;
	}
	double	val[3];
}VECTOR3D, *PVECTOR3D;

typedef struct _VERTEX_
{
	_VERTEX_() {
		val[0] = 0; val[1] = 0; val[2] = 0; val[3] = 1;
		nor = VECTOR3D(0, 0, 0);
		ntriCnt = 0;
		optIndex = -1;
	}

	_VERTEX_(double x, double y, double z) {
		val[0] = x; val[1] = y; val[2] = z; val[3] = 1;
		nor = VECTOR3D(0, 0, 0);
		ntriCnt = 0;
		optIndex = -1;
	}
	double		val[4];
	double		uv[2];
	VECTOR3D	nor;
	int			ntriCnt;
	int			optIndex;

	void inverseDirection() {
		this->val[0] = -this->val[0];
		this->val[1] = -this->val[1];
		this->val[2] = -this->val[2];
	}
}VERTEX, *PVERTEX;

typedef struct _TRIANGLE_
{
	_TRIANGLE_() {
		idx[0] = -1; idx[1] = -1; idx[2] = -1;
	}

	_TRIANGLE_(int idx1, int idx2, int idx3) {
		idx[0] = idx1; idx[1] = idx2; idx[2] = idx3;
	}
	int		idx[3];
	double	param[4];
	void flip() {
		int tmp = idx[1]; idx[1] = idx[2]; idx[2] = tmp;
	}
}TRIANGLE, *PTRIANGLE;

typedef struct _SUBNODE_
{
	vector<VERTEX>		vertex;
	vector<TRIANGLE>	index;
}SUBNODE, *PSUBNODE;

typedef struct _MATRIX_ 
{
	double val[4][4];
}MATRIX, *PMATRIX;

enum FBXENUM
{
	SUCCESS_IMPORT,
	INVALID_FBX_FILE,
	UNSURPPORTED_FBX_VERSION,
	UNSURPPORTED_PROPERTY_TYPE,
	INVALID_VERTEX_PAIR,
	INVALID_POLYGON_PAIR,
	COMPRESSED_PROPERTY_FOUND
};

extern Vec2i	neighbor[8];
extern SUBNODE	gem;
extern SUBNODE	prong;
extern int RESIZED_WIDTH;
extern int RESIZED_HEIGHT;
extern int ROI;

extern double	*depth;
extern int		*indexMap1;
extern int		*indexMap2;
extern VERTEX	*vertMap;
extern unsigned char	*renderImg;

extern double	roi_Depth;
extern int		gridInv;
extern double	maxDepth;
extern double	minDepth;
extern double	x_Rotate;
extern double	y_Rotate;
extern double	z_Rotate;
extern double	zoomFactor;
extern double	normalFactor;

extern MATRIX s_Mat;
extern MATRIX t_Mat;
extern MATRIX rX_Mat;
extern MATRIX rY_Mat;
extern MATRIX rZ_Mat;
extern MATRIX TMat;

extern vector<SUBNODE>	tNodes;
extern vector<SUBNODE>	eNodes;

extern VECTOR3D tDirV;

void normalizeDepth();
void doProjection(VERTEX *a);
VERTEX doCross(VERTEX *a, VERTEX *b);
VECTOR3D doCross(VECTOR3D *a, VECTOR3D *b);
void doNormalize(VERTEX *v);
void doNormalize(VECTOR3D *v);
double doMagnitude(VECTOR3D v);
double doMagnitude(VERTEX v);
void doProjection(VECTOR3D *a);
VERTEX doMinus(VERTEX *a, VERTEX *b);
VECTOR3D doMinus(VECTOR3D *a, VECTOR3D *b);
void OptimizeMesh(vector<SUBNODE> *model, bool optFlag, float optRate);
void OptimizeMesh(SUBNODE *node, bool optFlag, float optRate);
void GetBound(VERTEX *mx, VERTEX *mn, SUBNODE *node);
void GetBound(VERTEX *mx, VERTEX *mn, vector<SUBNODE> *nodes);
void Move2Center(SUBNODE *node);
void DoScale(SUBNODE *node, double scaleFactor);
void Move2Center(vector<SUBNODE> *nodes);
void DoScale(vector<SUBNODE> *nodes, double scaleFactor);
void cloneNodes(SUBNODE *src, SUBNODE *dst);
void cloneNodes(vector<SUBNODE> *src, vector<SUBNODE> *dst);
void Move2Position(SUBNODE *node, double x, double y, double z);
void Move2Position(vector<SUBNODE> *nodes, double x, double y, double z);
void doSetVertexNormal(vector<SUBNODE> *model);
void xRotateMatrix(double alpha);
void yRotateMatrix(double alpha);
void zRotateMatrix(double alpha);
void doSubnodeTransform(PSUBNODE node, PMATRIX mat);
void doTransform(vector<SUBNODE> *nodes, MATRIX *mat);

void InitRenderBuffer(int width, int height);
void ReleaseRenderBuffer();
void generateBitmapImage(unsigned char* image, int height, int width, char* imageFileName);