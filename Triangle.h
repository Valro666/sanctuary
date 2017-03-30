#ifndef __triangle_H__
#define __MODEL_H__

#include <vector>
#include "geometry.h"
#include "tgaimage.h"

class Triangle {
private:
public:
	
	int x1,x2,x3;
	int y1,y2,y3;
	
	TGAColor full,edge1,edge2,edge3;
	
	Triangle(int x1,int y1,int x2,int y2,int x3,int y3);
	Triangle(int x1,int y1,int x2,int y2,int x3,int y3,TGAColor full);
	Triangle(int x1,int y1,int x2,int y2,int x3,int y3,TGAColor edge1,TGAColor edge2,TGAColor edge3);
	Triangle(int x1,int y1,int x2,int y2,int x3,int y3,TGAColor full,TGAColor edge1,TGAColor edge2,TGAColor edge3);
	~Triangle();
	
	
};

#endif //__MODEL_H__
