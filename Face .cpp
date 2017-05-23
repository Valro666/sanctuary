#include <vector>
#include "geometry.h"
#include "tgaimage.h"
#include "Face.h"


	Face::Face(float v,float vt,float vn,float v2,float vt2,float vn2,float v3,float vt3,float vn3){

	s1=v;
	s2=v2;
	s3=v3;
	n1=vn;
	n2=vn2;
	n3=vn3;
	t1=vt;
	t2=vt2;
	t3=vt3;

	}



