#ifndef __triangle_H__
#define __MODEL_H__

#include <vector>
#include "geometry.h"
#include "tgaimage.h"
#include "Face.h"

class Triangle {
private:
public:
	
	float x1,x2,x3;
	float y1,y2,y3;
	float z1,z2,z3;
	Vec3f s1,s2,s3;
	Face *face ;//=NULL;
	
	TGAColor full,edge1,edge2,edge3;
	
	Triangle(int x1,int y1,int x2,int y2,int x3,int y3);
	Triangle(int x1,int y1,int x2,int y2,int x3,int y3,TGAColor full);
	Triangle(int x1,int y1,int x2,int y2,int x3,int y3,TGAColor edge1,TGAColor edge2,TGAColor edge3);
	Triangle(int x1,int y1,int x2,int y2,int x3,int y3,TGAColor full,TGAColor edge1,TGAColor edge2,TGAColor edge3);
	Triangle(Vec2i s1,Vec2i s2,Vec2i s3,TGAColor full);
	Triangle(Vec3f ss1,Vec3f ss2,Vec3f ss3,TGAColor full);
	Triangle(Vec3f ss1,Vec3f ss2,Vec3f ss3,Face full);
	~Triangle();
	
	bool neg();
	
	void setZ(float z1,float z2,float z3);
	
	
};

#endif //__MODEL_H__
