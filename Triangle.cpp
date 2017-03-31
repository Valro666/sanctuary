#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <math.h>
#include "Triangle.h"

Triangle::	~Triangle(){	
}

Triangle::	Triangle(int x1,int y1,int x2,int y2,int x3,int y3){
	
	this->x1=x1;
	this->x2=x2;
	this->x3=x3;
	
	this->y1=y1;
	this->y2=y2;
	this->y3=y3;
	
	
}

Triangle::	Triangle(int x1,int y1,int x2,int y2,int x3,int y3, TGAColor full){
	
	this->x1=x1;
	this->x2=x2;
	this->x3=x3;
	
	this->y1=y1;
	this->y2=y2;
	this->y3=y3;
	
	this->full=full;
	this->edge1=full;
	this->edge2=full;
	this->edge3=full;
	
	
}

Triangle::	Triangle(int x1,int y1,int x2,int y2,int x3,int y3, TGAColor edge1, TGAColor edge2, TGAColor edge3){
	
	this->x1=x1;
	this->x2=x2;
	this->x3=x3;
	
	this->y1=y1;
	this->y2=y2;
	this->y3=y3;
	
	this->full=edge1;
	this->edge1=edge1;
	this->edge2=edge2;
	this->edge3=edge3;
	
	
}

Triangle::	Triangle(int x1,int y1,int x2,int y2,int x3,int y3,TGAColor full, TGAColor edge1, TGAColor edge2, TGAColor edge3){
	
	this->x1=x1;
	this->x2=x2;
	this->x3=x3;
	
	this->y1=y1;
	this->y2=y2;
	this->y3=y3;
	
	this->full=full;
	this->edge1=edge1;
	this->edge2=edge2;
	this->edge3=edge3;
	
	
}

Triangle::	Triangle(Vec2i s1,Vec2i s2,Vec2i s3,TGAColor full){
	x1=s1.x;
	x2=s1.x;
	x3=s3.x;
	
	y1=s1.y;
	y2=s2.y;
	y3=s3.y;
	
	this->full=full;
	this->edge1=full;
	this->edge2=full;
	this->edge3=full;
	
}

