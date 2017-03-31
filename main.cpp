
#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "triangle.h"
#include <stdio.h>   
#include <stdlib.h>    
#include <time.h>

const TGAColor white = TGAColor(255, 255, 255, 255);

const TGAColor rose = TGAColor(255, 0, 255, 255);

const TGAColor bleu = TGAColor(0, 0, 255, 255);

const TGAColor cyan = TGAColor(0, 100, 100, 100);

const int width  = 800;
const int height = 800;

TGAColor randC(){
	int a = rand() * 255 ;
	int b = rand() * 255 ;
	int c = rand() * 255;
	
	return TGAColor(0, a, b, c);
}

void line(TGAImage &img){

    for(int i = 0 ; i < 100;i++){
        img.set(i,i,white);
       // img.set(i,i+1,white);
       // img.set(i,i+2,white);
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
    for (int x=x0; x<=x1; x++) { 
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

void triangle(Triangle tri,TGAImage &image){
	
	
	
	int x0=tri.x1 ,x1=tri.x2 ,x2=tri.x3 ;
	int y0=tri.y1 ,y1=tri.y2 ,y2=tri.y3 ;
	TGAColor color = tri.full;	
	
	line3(image,tri.edge1,x0,y0,x1,y1);
	
	line3(image,tri.edge2,x0,y0,x2,y2);
	
	line3(image,tri.edge3,x2,y2,x1,y1);
	
	int tmp=-1;
     
    for (int x=x0; x<=x1; x++) { 
        float t = (x-x0)/(float)(x1-x0); 
        float y = y0*(1.-t) + y1*t; 
    		line3(image,color,x,y,x2,y2);
		tmp = (int)y;  
    }
}

int main(int argc, char** argv) {
	
	srand(time(NULL));
	
    TGAImage image(800, 600, TGAImage::RGB); // 130 200 100 300
    /*
    Triangle *tri= new Triangle(10,10,130,200,100,300,randC(),randC(),randC(),randC());   
	Triangle *tri2= new Triangle(10,10,150,150,350,250,randC(),randC(),randC(),randC());
	Triangle *tri3= new Triangle(450,550,600,600,350,250,randC(),randC(),randC(),randC());
	
    triangle(*tri,image);
	triangle(*tri2,image);
	triangle(*tri3,image);
	
	/*
	line2(image,rose,10, 10, 150, 150);
	line2(image,rose,150, 150, 350, 250);
	line2(image,rose,10, 10, 350, 250);
	
    line3(image,bleu,600, 400, 130, 200);
    line3(image,bleu,130, 200, 100, 300);
    line3(image,bleu,600, 400, 100, 300);
    
//*/
	
	Model *model = NULL;
	
	 model = new Model("objet/african_head.obj");
	for (int i=0; i<model->nfaces(); i++) { 
    std::vector<int> face = model->face(i); 
    for (int j=0; j<3; j++) { 
        Vec3f v0 = model->vert(face[j]); 
        Vec3f v1 = model->vert(face[(j+1)%3]); 
        int x0 = (v0.x+1.)*width/2.; 
        int y0 = (v0.y+1.)*height/2.; 
        int x1 = (v1.x+1.)*width/2.; 
        int y1 = (v1.y+1.)*height/2.; 
        line( image, white,x0, y0, x1, y1); 
    }
    
    
	for (int i=0; i<model->nfaces(); i++) { 
	
    std::vector<int> face = model->face(i); 
    Vec2i screen_coords[3]; 
    for (int j=0; j<3; j++) { 
    
        Vec3f world_coords = model->vert(face[j]); 
        screen_coords[j] = Vec2i((world_coords.x+1.)*width/2., (world_coords.y+1.)*height/2.); 
    } 
    
    Triangle *tri = new Triangle(screen_coords[0].x,screen_coords[0].y, 
	screen_coords[1].x,screen_coords[1].y, 
	screen_coords[2].x,screen_coords[2].y,randC());
    triangle(*tri, image); 
}
	 
}
	//*/
	
	image.flip_vertically();
    image.write_tga_file("dessin.tga");
    return 0;
}
