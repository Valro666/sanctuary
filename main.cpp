
#include <vector>
#include <cmath>
#include <string>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
//#include "geometry3.h"
#include "triangle.h"
#include <stdio.h>   
#include <vector>   
#include <stdlib.h>    
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <complex>
#include <limits>

const TGAColor white = TGAColor(255, 255, 255, 255);

const TGAColor rose = TGAColor(255, 0, 255, 255);

const TGAColor bleu = TGAColor(0, 0, 255, 255);

const TGAColor cyan = TGAColor(0, 100, 100, 100);

 float Zcam = 10;
 
float za=0, zb=0 , zc=0;

const int width  = 800;
const int height = 800;

//typedef vec<2,  int> Vec2f;

TGAColor randC(){
	int a = rand() * 255 ;
	int b = rand() * 255 ;
	int c = rand() * 255;
	
	return TGAColor(a, b, c, 255);
}

TGAColor randC(int encity){
	int a = rand() * 255 ;
	int b = rand() * 255 ;
	int c = rand() * 255;
	
	return TGAColor(a*encity, b*encity, c*encity, 255);
}

void line(TGAImage &img){

    for(int i = 0 ; i < 100;i++){
        img.set(i,i,white);
    }
}
//attention a la taille entre chaque pixel
void line(TGAImage &img,TGAColor c,int x,int y ,int xx,int yy){

    for (float t=0.; t<1.; t+=.01) {
            int x0 = x*(1.-t) + xx*t;
            int y0 = y*(1.-t) + yy*t;
            img.set(x0, y0, c);
        }
}

void line2(TGAImage &image,TGAColor c,int x0,int y0 ,int x1,int y1){

    for (int x=x0; x<=x1; x++) {
            float t = (x-x0)/(float)(x1-x0);
            int y = y0*(1.-t) + y1*t;
            image.set(x, y, c);
        }
}

void line3(TGAImage &image, TGAColor color, int x0, int y0, int x1, int y1) { 
	int tmp=-1;
    bool steep = false; 
    if (std::abs(x0-x1)<std::abs(y0-y1)) { // if the line is steep, we transpose the image 
        std::swap(x0, y0); 
        std::swap(x1, y1); 
        steep = true; 
    } 
    if (x0>x1) { // make it left-to-right 
        std::swap(x0, x1); 
        std::swap(y0, y1); 
    } 
    for (int x=x0; x<x1; x++) { 
    if(x<0){
    	break;
	}
        float t = (x-x0)/(float)(x1-x0); 
        float y = y0*(1.-t) + y1*t;
		
        if (steep) { 
            image.set(y, x, color); // if transposed, de-transpose 
        } else { 
            image.set(x, y, color); 
        }
		tmp = (int)y;  
    } 
}

void fond(TGAImage &image){

    for(int i = 0 ; i < image.get_width();i++){
        for(int z = 0 ; z < image.get_height();z++){
             image.set(i,z,rose);
        }
    }

}

void drawModel(Model mod,TGAImage img){
	
	
}
 bool bary(Triangle tri,int x , int y){
	
	
	double denominator= ((tri.y2 - tri.y3)*(tri.x1 - tri.x3) + (tri.x3 - tri.x2)*(tri.y1 - tri.y3));
 	double a = ((tri.y2 - tri.y3)*(x - tri.x3) + (tri.x3 - tri.x2)*(y - tri.y3)) / denominator;
 	double b = ((tri.y3 - tri.y1)*(x - tri.x3) + (tri.x1 - tri.x3)*(y - tri.y3)) / denominator;
 	double c = 1 - a - b;
	
	if(a>0)
		if(a<1)
			if(b>0)
				if(b<1)
					if(c>0)
						if(c<1) return true;
						
	return false;
	
}

bool bary2(Triangle tri,int x , int y){
	
	
	float denominator= ((tri.y2 - tri.y3)*(tri.x1 - tri.x3) + (tri.x3 - tri.x2)*(tri.y1 - tri.y3));
 	float a = ((tri.y2 - tri.y3)*(x - tri.x3) + (tri.x3 - tri.x2)*(y - tri.y3)) / denominator;
 	float b = ((tri.y3 - tri.y1)*(x - tri.x3) + (tri.x1 - tri.x3)*(y - tri.y3)) / denominator;
 	float c = 1 - a - b;
	za = a;
	zb = b;
	zc = c;
	
//	printf("%f %f %f\n",za,zb,zc);
	if(a>=0.)
		if(a<=1.)
			if(b>=0)
				if(b<=1.)
					if(c>=0)
						if(c<=1.) return true;
						
	return false;
	
}

void triangle2(Triangle tri,TGAImage &image){
	
	int x0=tri.x1 ,x1=tri.x2 ,x2=tri.x3 ;
	int y0=tri.y1 ,y1=tri.y2 ,y2=tri.y3 ;
	TGAColor color = tri.full;	
	
	line3(image,tri.edge1,x0,y0,x1,y1);
	line3(image,tri.edge2,x0,y0,x2,y2);
	line3(image,tri.edge3,x2,y2,x1,y1);
//*/
	int xmin=x0 ,xmax=x2;
	int ymin=y0 ,ymax=y2;
	
	
	if(xmin>x1) xmin = x1;
	if(xmin>x2) xmin = x2;
	if(ymin>y1) ymin = y1;
	if(ymin>y2) ymin = y2;
	
	if(xmax<x0) xmax = x0;
	if(xmax<x1) xmax = x1;
	if(ymax<y0) ymax = y0;
	if(ymax<y1) ymax = y1;
	
	for(int x = xmin;x < xmax;x++)
		for(int y = ymin;y<ymax;y++)
			if(bary(tri,x,y)) image.set(x,y,color);
			
	
}

void triangleBuff(Triangle tri,  float *buffer,TGAImage &image){
	
	int x0=tri.x1 ,x1=tri.x2 ,x2=tri.x3 ;
	int y0=tri.y1 ,y1=tri.y2 ,y2=tri.y3 ;
	TGAColor color = tri.full;	
	
//	line3(image,tri.edge1,x0,y0,x1,y1);
//	line3(image,tri.edge2,x0,y0,x2,y2);
//	line3(image,tri.edge3,x2,y2,x1,y1);
	int xmin=x0 ,xmax=x2;
	int ymin=y0 ,ymax=y2;
	
	if(xmin>x1) xmin = x1;
	if(xmin>x2) xmin = x2;
	if(ymin>y1) ymin = y1;
	if(ymin>y2) ymin = y2;
	
	if(xmax<x0) xmax = x0;
	if(xmax<x1) xmax = x1;
	if(ymax<y0) ymax = y0;
	if(ymax<y1) ymax = y1;
	
	 int tp = false ;
	
	for(int x = xmin;x < xmax;x++)
		for(int y = ymin;y<ymax;y++){
			za=0, zb=0 ,zc=0 ;
			if(bary2(tri,x,y)){
				float z = tri.z1*za+tri.z2*zb+tri.z3*zc;
				float mod =(1.-(z/Zcam));
				//float mod =1;
				z=z/(1-z/Zcam);
				
				float xx , yy;
				
				
				
				xx = x/mod;
				yy = y/mod;
				
				
				//printf("%f %f %f %d %f %f %f %f\n",z,Zcam,mod,x,xx,za,zb,zc);
				
				if(int(xx+yy*width)<=(width*height))
				if (buffer[int(xx+yy*width)]<=z) {
                buffer[int(xx+yy*width)] = z;
                image.set(xx,yy, color);
            	}
			}
		}
}
/*
Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
    Vec3f s[2];
    for (int i=2; i--; ) {
        s[i][0] = C[i]-A[i];
        s[i][1] = B[i]-A[i];
        s[i][2] = A[i]-P[i];
    }
    Vec3f u = cross(s[0], s[1]);
    return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    //if (std::abs(u[2])>1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
    //    return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    //return Vec3f(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}
*/
void triangle(Triangle tri,TGAImage &image){
	
	int x0=tri.x1 ,x1=tri.x2 ,x2=tri.x3 ;
	int y0=tri.y1 ,y1=tri.y2 ,y2=tri.y3 ;
	TGAColor color = tri.full;	
	int tmp=-1;
    for (int x=x0; x<=x1; x++) { 
    if(x<0){
    	break;
	}
        float t = (x-x0)/(float)(x1-x0); 
        float y = y0*(1.-t) + y1*t; 
    	line3(image,color,x,y,x2,y2);
		tmp = (int)y; 
    }
	line3(image,white,x0,y0,x1,y1);
	line3(image,white,x0,y0,x2,y2);
	line3(image,white,x2,y2,x1,y1);
}

//double buffer[width][height];

int main(int argc, char** argv) {
	
	
	
	Model *model = NULL;
	
	model = new Model("objet/diablo.obj");
	
	if(argv[1]!=NULL){
	Zcam = atoi(argv[1]);
	}
	
	srand(time(NULL));
	
    TGAImage image(800, 800, TGAImage::RGB); // 130 200 100 300
    /*
    Triangle *tri= new Triangle(10,10,130,200,100,300,randC(),randC(),randC(),randC());   
	Triangle *tri2= new Triangle(10,10,150,150,350,250,randC(),randC(),randC(),randC());
	Triangle *tri3= new Triangle(450,550,600,600,350,250,randC(),randC(),randC(),randC());
	
    triangle2(*tri,image);
	triangle2(*tri2,image);
	triangle2(*tri3,image);
	
	/*
	line2(image,rose,10, 10, 150, 150);
	line2(image,rose,150, 150, 350, 250);
	line2(image,rose,10, 10, 350, 250);
	
    line3(image,bleu,600, 400, 130, 200);
    line3(image,bleu,130, 200, 100, 300);
    line3(image,bleu,600, 400, 100, 300);
    
//*/
	
	
	Vec3f light_dir(0,0,-1);
    //std::vector<int> buffer(width*height,-std::numeric_limits<int>::max());
	//printf("fin buff \n");
	float *buffer = new float[width*height];
    for (int i=width*height; i--; buffer[i] = -std::numeric_limits<float>::max());
	
	for (int i=0; i<model->nfaces(); i++) { 
    	std::vector<int> face = model->face(i); 
    	Vec3f screen_coords[3]; 
    	Vec3f world_coords[3]; 
        
   		for (int j=0; j<3; j++) { 
       		Vec3f v = model->vert(face[j]); 
       		screen_coords[j] = Vec3f((v.x+1.)*width/2., (v.y+1.)*height/2.,v.z); 
       		world_coords[j]  = v; 
		}
    
		Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]); 
    	n.normalize(); 
    	float intensity = n*light_dir; 
   
   		if (intensity>0) { /*
        	Triangle *tri = new Triangle(screen_coords[0].x,screen_coords[0].y, 
			screen_coords[1].x,screen_coords[1].y, 
			screen_coords[2].x,screen_coords[2].y,TGAColor(intensity*255, intensity*255, intensity*255, 255));
			tri->setZ(screen_coords[0].z,screen_coords[1].z,screen_coords[2].z);//*/
			
			
			Triangle *tri = new Triangle(screen_coords[0],screen_coords[1],screen_coords[2],
			TGAColor(intensity*255, intensity*255, intensity*255, 255));//*/
			
			
			triangleBuff(*tri,buffer, image); 
    	} 
	}
	
	image.flip_vertically();
    image.write_tga_file("dessin.tga");
    return 0;
}
