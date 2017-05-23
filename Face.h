#ifndef __Face_H__
#define __Face_H__

#include <vector>
#include "geometry.h"
#include "tgaimage.h"

class Face {
private:
public:

	float s1,s2,s3;
	float n1,n2,n3;
	float t1,t2,t3;

	//Face(int i,char *filename);
	Face(float x,float y,float z,float xt,float yt,float zt,float xv,float yv,float zv);


};

#endif //__Face_H__
