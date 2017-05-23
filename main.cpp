
#include <vector>
#include <cmath>
#include <string>
#include "tgaimage.h"
#include "model.h"
#include "geometry3.h"
//#include "geometry.h"
//#include "geometry3.h"
#include "triangle.h"
#include <stdio.h>
#include "Face.h"
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

std::vector<Vec3f> sommet;
std::vector<Vec3f> texture;
std::vector<Vec3f> vecteur;
std::vector<Face> modele;

Vec3f eye(0.25,0,1); // x y distance : position oeil  [0,1]
Vec3f center(0,0,0); // x y z : coordonnée focus de l oeil [0,1]
Vec3f light_src(0,0,-1); // x y z : [0,1]
Vec3f varying_intensity;

Matrix ModelView  ;
Matrix Projection;
Matrix ViewPort  ;

TGAImage diffu ;

void load_texture(std::string filename, const char *suffix, TGAImage &img) {
    std::string texfile(filename);
    size_t dot = texfile.find_last_of(".");
    if (dot!=std::string::npos) {
        texfile = texfile.substr(0,dot) + std::string(suffix);
        std::cerr << "texture file " << texfile << " loading " << (img.read_tga_file(texfile.c_str()) ? "ok" : "failed") << std::endl;
        img.flip_vertically();
    }
}

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

}

void triangleBuff(Triangle tri,  float *buffer,TGAImage &image){

	int x0=tri.x1 ,x1=tri.x2 ,x2=tri.x3 ;
	int y0=tri.y1 ,y1=tri.y2 ,y2=tri.y3 ;
	TGAColor color = tri.full;

	//line3(image,tri.edge1,x0,y0,x1,y1);
	//line3(image,tri.edge2,x0,y0,x2,y2);
	//line3(image,tri.edge3,x2,y2,x1,y1);
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
				//float mod =(1.-(z/Zcam));
				float mod =1;
				//z=z/(1-z/Zcam);
				float xx , yy;
				xx = x/mod;
				yy = y/mod;
				//printf("%f %f %f %d %f %f %f %f\n",z,Zcam,mod,x,xx,za,zb,zc);
				if(int(xx+yy*width)<=(width*height))
				if (buffer[int(xx+yy*width)]<=z) {
					TGAColor couleur = TGAColor(0, 100, 100, 100);
					buffer[int(xx+yy*width)] = z;
					image.set(xx,yy, color);
            	}
			}
		}
}

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

Matrix viewport(int x, int y, int w, int h) {
    Matrix m = Matrix::identity(4);
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = Zcam/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = Zcam/2.f;
    return m;
}

Matrix lookat(Vec3f eye, Vec3f center, Vec3f up) {
    Vec3f z = (eye-center).normalize();
    Vec3f x = (up^z).normalize();
    Vec3f y = (z^x).normalize();
    Matrix res = Matrix::identity(4);
    for (int i=0; i<3; i++) {
        res[0][i] = x[i];
        res[1][i] = y[i];
        res[2][i] = z[i];
        res[i][3] = -center[i];
    }
    return res;
}

Vec3f normale(Triangle tri){

	
    float X = (tri.y2-tri.y1)*(tri.z3-tri.z1) - (tri.z2-tri.z1)*(tri.y3-tri.y1);
    float Y = (tri.z2-tri.z1)*(tri.x3-tri.x1) - (tri.x2-tri.x1)*(tri.z3-tri.z1);
    float Z = (tri.x2-tri.x1)*(tri.y3-tri.y1) - (tri.y2-tri.y1)*(tri.x3-tri.x1);

    float norme= X*X+Y*Y+Z*Z;

    norme = sqrt((norme));

    X=X/norme;
    Y=Y/norme;
    Z=Z/norme;
	
	//Vec3f n = (vecteur[tri.face->n1]+vecteur[tri.face->n2]+vecteur[tri.face->n3]);
	//n.normalize();

	Vec3f ff(X,Y,Z);
	ff.normalize();
	return ff;
//return n;
}

float norme(Vec3f vec){

    float norme= vec.x*vec.x+vec.y*vec.y+vec.z*vec.z;

    norme = sqrt((norme));

return norme;
}

Vec3f cross(Vec3f u,Vec3f v){
	float x= (u.y*v.z)-(u.z*v.y);
	float y= (u.z*v.x)-(u.x*v.z);
	float z= (u.x*v.y)-(u.y*v.x);


return Vec3f(x,y,z);
}

float scalaire(Vec3f vec,Vec3f vec2){

    float norme= vec.x*vec2.x+vec.y*vec2.y+vec.z*vec2.z;

return norme;
}




void aff(Vec3f pts,TGAImage &image){
	Vec3f truc = pts;
	int xmin = pts.x-50;
	int xmax = pts.x+50;
	int ymin = pts.y-50;
	int ymax = pts.y+50;

	if(xmin<0){
		xmin = 0;
	}
	if(xmax>width){
		xmax = width;
	}

	if(ymax>height){
		ymax = height;
	}

	if(ymin<0){
		ymin = 0;
	}

	for(int x = xmin ; x < xmax ; x++){
		for(int y = ymin ; x < ymax ; y++){
			image.set(x,y, rose);
			printf("%f %f %d %d\n",x,y,width,height);
		}
	}
}

Vec3f mult(Vec3f p , float i){
	float x = p.x*i;
	float y = p.y*i;
	float z = p.z*i;
	return Vec3f(x,y,z);
}

Vec3f normale2(Triangle tri){
	Vec3f f1 = mult(vecteur[tri.face->n1],za);
	Vec3f f2 = mult(vecteur[tri.face->n2],zb);
	Vec3f f3 = mult(vecteur[tri.face->n3],zc);

	Vec3f tot = f1+f2+f3;
	return tot;
}

Vec3f reflection(Triangle tri,Vec3f light){
    float X = 0;
    float Y = 0;
    float Z = 0;

    Vec3f norm ;
    norm= normale(tri);

    float val= scalaire(normale(tri),light);
     val *=2;

     norm.x = norm.x * val;
     norm.y = norm.y * val;
     norm.z = norm.z * val;

     norm = norm - light;

return norm;
}

void triangleBuffIntense(Triangle tri,  float *buffer,TGAImage &image){

	float intensity;

	float x0=tri.x1 ,x1=tri.x2 ,x2=tri.x3 ;
	float y0=tri.y1 ,y1=tri.y2 ,y2=tri.y3 ;

	TGAColor color = tri.full;

	float xmin=x0 ,xmax=x2;
	float ymin=y0 ,ymax=y2;

	if(xmin>x1) xmin = x1;
	if(xmin>x2) xmin = x2;
	if(ymin>y1) ymin = y1;
	if(ymin>y2) ymin = y2;

	if(xmax<x0) xmax = x0;
	if(xmax<x1) xmax = x1;
	if(ymax<y0) ymax = y0;
	if(ymax<y1) ymax = y1;

    int tp = false ;

	for(int x = xmin;x < xmax;x++){
            int g = -1;
            int d = 0;
		for(int y = ymin;y<ymax;y++){
            if(bary2(tri,x,y)){
                if(g==-1){
                    g=y;
                }
                d = y;
                }
            }

				Vec3f p1 = tri.s1;
				Vec3f p2 = tri.s2;
				Vec3f p3 = tri.s3;
				Vec3f tmp = tri.s3;

				if(p1.x<p2.x){
					tmp = p1;
					p1 = p2;
					p2 = tmp;
				}

				if(p1.x<p3.x){
					tmp = p1;
					p1 = p3;
					p3 = tmp;
				}

        for(int y = ymin;y<ymax;y++){

			za=0, zb=0 ,zc=0 ;
			if(bary2(tri,x,y)){

				float z = tri.z1*za+tri.z2*zb+tri.z3*zc;
				float mod =1;
				float xx , yy;
				xx = x/mod;
				yy = y/mod;
				float tx = x0*za + x1*zb + x2*zc;
				float ty = x0*za + x1*zb + x2*zc;

				//coordonnée source lumiere
				Vec3f source = light_src;
				source.normalize();//= Vec3f(ViewPort*Projection*ModelView*Matrix(light_src));
				// normale triangle
				Vec3f N = normale(tri);
				//Vec3f N = normale2(tri);
				N.normalize();
				// reflection lumiere
				Vec3f ref = reflection(tri,source);
				ref.normalize();

                //k
				float renvoie=0.5;
				float diffuse=1;
				float specule=0;
				float brillance=10;
				//i
                float ir=10;
				float id=125;
				float is=10;

				float R = renvoie*ir;
				float D = diffuse*id*scalaire(source,N);
				float S = specule*is*pow(scalaire(ref,eye),brillance);

                intensity = R+D+S;

				if (buffer[int(xx+yy*width)]<=z) {

					Vec3f tt = texture[(tri.face->t1)];
					Vec3f tt2 = texture[(tri.face->t2)];
					Vec3f tt3 = texture[(tri.face->t3)];
					int tx= tt.x-1+tt2.x-1;
					tx/=2;
					int ty= tt.y-1+tt2.y-1;
					ty/=2;
					
					TGAColor couleur = TGAColor(255*intensity, 255*intensity, 255*intensity,255);
					
					buffer[int(xx+yy*width)] = z;
					//printf("%f\n",tri.s1.x);

					// recuperation de la texture mais ca marche pas 
					int xture = za*texture[tri.face->t1-1].x+ zb*texture[tri.face->t2-1].x+ zc*texture[tri.face->t3-1].x;
					int yture = za*texture[tri.face->t1-1].y+ zb*texture[tri.face->t2-1].y+ zc*texture[tri.face->t3-1].y;
					TGAColor color = diffu.get(xture*1.25,height-(yture*0.50));
					float rouge = color.r;
					float vert = color.g;
					float bleu = color.b;
					rouge*=intensity;
					vert*=intensity;
					bleu*=intensity;
					
					//printf("%f %c %c\n",rouge,vert,bleu);
					image.set(xx,yy, couleur);
					//image.set(xx,yy, color);
					//image.set(xx,yy, TGAColor(rouge,vert,bleu,255));
            	}
			}
		}
	}
}


int nombreFace(const char *filename){

	int res = 0 ;
	std::ifstream in;
	in.open (filename, std::ifstream::in);
	if (in.fail()) return res;
	std::string line;
	while (!in.eof()) {
		std::getline(in, line);
        std::istringstream iss(line.c_str());
		if (!line.compare(0, 2, "v ")){
			line = line.substr(2);
			int space = line.find_first_of(" ");
			float x = atof(line.substr(0,space).c_str());
			line = line.substr(space+1);
			space = line.find_first_of(" ");
			float y = atof(line.substr(0,space).c_str());
			line = line.substr(space+1);
			space = line.find_first_of(" ");
			float z = atof(line.substr(0,space).c_str());
			//printf("%f %f %f\n",x,y,z);
			Vec3f truc = Vec3f(x,y,z);
			sommet.push_back(truc);
		}
		if (!line.compare(0, 3, "vt ")){
			//printf("%s\n",line.c_str());
			line = line.substr(4);
			int space = line.find_first_of(" ");
			float x = atof(line.substr(0,space).c_str());
			line = line.substr(space+1);
			space = line.find_first_of(" ");
			x = x * width;
			float y = atof(line.substr(0,space).c_str());
			line = line.substr(space+1);
			space = line.find_first_of(" ");
			y = y * height;
			float z = atof(line.substr(0,space).c_str());
			//printf("%f %f %f\n",x,y,z);
			Vec3f truc = Vec3f(x,y,z);
			texture.push_back(truc);
		}
		if (!line.compare(0, 3, "vn ")){
			//printf("%s\n",line.c_str());
			line = line.substr(4);
			int space = line.find_first_of(" ");
			float x = atof(line.substr(0,space).c_str());
			line = line.substr(space+1);
			space = line.find_first_of(" ");
			float y = atof(line.substr(0,space).c_str());
			line = line.substr(space+1);
			space = line.find_first_of(" ");
			float z = atof(line.substr(0,space).c_str());
			//printf("%f %f %f\n",x,y,z);
			Vec3f truc = Vec3f(x,y,z);
			vecteur.push_back(truc);
		}

		if (!line.compare(0, 2, "f ")){

			line = line.substr(2);
			int space = line.find_first_of("/");
			float x = atof(line.substr(0,space).c_str());
			line = line.substr(space+1);
			space = line.find_first_of("/");
			float y = atof(line.substr(0,space+1).c_str());
			line = line.substr(space+1);
			space = line.find_first_of(" ");
			float z = atof(line.substr(0,space+1).c_str());
			line = line.substr(space+1);

			space = line.find_first_of("/");
			float xt = atof(line.substr(0,space).c_str());
			line = line.substr(space+1);
			space = line.find_first_of("/");
			float yt = atof(line.substr(0,space+1).c_str());
			line = line.substr(space+1);
			space = line.find_first_of(" ");
			float zt = atof(line.substr(0,space+1).c_str());
			line = line.substr(space+1);

			space = line.find_first_of("/");
			float xv = atof(line.substr(0,space).c_str());
			line = line.substr(space+1);
			space = line.find_first_of("/");
			float yv = atof(line.substr(0,space+1).c_str());
			line = line.substr(space+1);
			float zv = atof(line.c_str());

			modele.push_back(Face(x,y,z,xt,yt,zt,xv,yv,zv));
			//line = line.substr(line.length()-1);
			//printf("%s\n%f %f %f %f %f %f %f %f %f\n",line.c_str(),x,y,z,xt,yt,zt,xv,yv,zv);
			res++;
		}
    }
	//printf("%d\n",res);
	return res;
}


int main(int argc, char** argv) {

	diffu.read_tga_file("objet/african_head_diffuse.tga");

	int nb = nombreFace("objet/african_head.obj");
	//int nb = nombreFace("objet/diablo.obj");

	ModelView  = lookat(eye, center, Vec3f(0,1,0));
	Matrix Projection = Matrix::identity(4);
	Matrix ViewPort   = viewport(width/8, height/8, width*3/4, height*3/4);

	srand(time(NULL));

    TGAImage image(width, height, TGAImage::RGB);
	
	float *buffer = new float[width*height];
    for (int i=width*height; i--; buffer[i] = -std::numeric_limits<float>::max());

	for(int i=0; i<nb; i++){

		Face f = modele[i];
		Vec3f s1=sommet[f.s1-1];
		s1 = Vec3f(ViewPort*Projection*ModelView*Matrix(s1));
		Vec3f s2=sommet[f.s2-1];
		s2 = Vec3f(ViewPort*Projection*ModelView*Matrix(s2));
		Vec3f s3=sommet[f.s3-1];		
		s3 = Vec3f(ViewPort*Projection*ModelView*Matrix(s3));

		Triangle *tri = new Triangle(s1,s2,s3,f);//
		triangleBuffIntense(*tri,buffer, image);

	}

	image.flip_vertically();
    image.write_tga_file("dessin.tga");
    return 0;
}
