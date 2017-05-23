
#include <vector>
#include <cmath>
#include <string>
#include "tgaimage.h"
#include "model.h"
#include "geometry3.h"
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

Vec3f eye(0.5,0.5,1); // x y distance : position oeil  [0,1]
Vec3f center(0,0,0); // x y z : coordonnée focus de l oeil [0,1]
Vec3f light_src(0,0,-0.5); // x y z : [0,1]
Vec3f varying_intensity;

Matrix ModelView  ;
Matrix Projection;
Matrix ViewPort  ;

Model *model = NULL;

TGAImage diffu ;

//typedef vec<2,  int> Vec2f;
/*
virtual Vec4f vertex(int iface, int nthvert) {
        varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        return Viewport*Projection*ModelView*gl_Vertex; // transform it to screen coordinates
    }*/



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
/*
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
			if(bary(tri,x,y)) image.set(x,y,color);//*/


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
	
	Vec3f n = (vecteur[tri.face->n1]+vecteur[tri.face->n2]+vecteur[tri.face->n3]);
	n.normalize();

return Vec3f(X,Y,Z);
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

Vec3f reflection(Triangle tri,Vec3f light){
    float X = 0;
    float Y = 0;
    float Z = 0;

    Vec3f norm ;
    norm= normale(tri);
	//Vec3f norm = (vecteur[tri.face->n1]+vecteur[tri.face->n2]+vecteur[tri.face->n3])*-1;

    float val= scalaire(normale(tri),light);
     val *=2;

     norm.x = norm.x * val;
     norm.y = norm.y * val;
     norm.z = norm.z * val;

     norm = norm - light;

return norm;
}


void aff(Vec3f pts,TGAImage &image){
	Vec3f truc = pts;//(ViewPort*Projection*ModelView*Matrix(pts));
	//truc.normalize();
	//printf("%f %f %f \n",pts.x,pts.y,pts.z);
	//printf("%f %f %f \n",truc.x,truc.y,truc.z);
	

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
void triangleBuffIntense(Triangle tri,  float *buffer,TGAImage &image){
	//if (tri.y1==tri.y2 && tri.y1==tri.y3) return;

	//tri.s1 = Vec3f(ViewPort*Projection*ModelView*Matrix(tri.s1));
	//tri.s2 = Vec3f(ViewPort*Projection*ModelView*Matrix(tri.s2));
	//tri.s3 = Vec3f(ViewPort*Projection*ModelView*Matrix(tri.s3));

	//Vec3f s1 = Vec3f(ViewPort*Projection*ModelView*Matrix(tri.s1));
	//Vec3f s2 = Vec3f(ViewPort*Projection*ModelView*Matrix(tri.s2));
	//Vec3f s3 = Vec3f(ViewPort*Projection*ModelView*Matrix(tri.s3));

	//printf("%f %f %f\n",tri.s1.x,tri.s2.x,tri.s3.x);

	float intensity;
	//int x0=tri.s1.x ,x1=tri.s2.x ,x2=tri.s3.x ;
	//int y0=tri.s1.y ,y1=tri.s1.y ,y2=tri.s1.y ;

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

				//Vec3f norme = (vecteur[tri.face->n1]+vecteur[tri.face->n2]+vecteur[tri.face->n3]);
				Vec3f norme ;
				norme = normale(tri);
				norme.normalize();
				Vec3f light_dst(xx,yy,z);
				// Vec3f light_dst = tri.s1;
				/* 
				Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]); 
				n.normalize(); 
				float intensity = n*light_dir;
				*/
				//Vec3f light = cross(light_src,light_dst);
				Vec3f light = light_src;
				light.normalize();
				Vec3f reflec ;
				reflec = reflection(tri,light);

				//coordonnée source lumiere
				Vec3f source = light_src;
				source.normalize();//= Vec3f(ViewPort*Projection*ModelView*Matrix(light_src));
				// normale triangle
				Vec3f N = normale(tri);
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
				//float D = diffuse*id*scalaire(light_dir,norme);
				float D = diffuse*id*scalaire(source,N);
				float S = specule*is*pow(scalaire(ref,eye),brillance);
				//float S = specule*is*pow(brillance,scalaire(reflec,eye));

                intensity = R+D+S;
				//intensity = R;

				//intensity = 1;

				//if(int(xx+yy*width)<(width*height))
				if (buffer[int(xx+yy*width)]<=z) {

					Vec3f tt = texture[(tri.face->t1)];
					Vec3f tt2 = texture[(tri.face->t2)];
					Vec3f tt3 = texture[(tri.face->t3)];
					int tx= tt.x-1+tt2.x-1;
					tx/=2;
					int ty= tt.y-1+tt2.y-1;
					ty/=2;

					//printf("%f %f\n",tx,ty);

					//int r = diffu.get(tx,ty).r*intensity;
					//int v = diffu.get(tx,ty).g*intensity;
					//int b = diffu.get(tx,ty).b*intensity;
					//Vec3f truc = ViewPort*Projection*ModelView*Matrix(light_src);
					//TGAColor couleur = TGAColor(r,v,b,255);
					//TGAColor couleur = diffu.get(truc.x,truc.y);
					//printf("%f %f %f\n",r,v,b);
					
					TGAColor couleur = TGAColor(255*intensity, 255*intensity, 255*intensity,255);
					
					//TGAColor couleur = TGAColor(255*z, 255*z, 255*z,255);
					buffer[int(xx+yy*width)] = z;
					image.set(xx,yy, couleur);
					//image.set(truc.x,truc.y, rose);
					//image.set(xx,yy, color);
            	}
			}
		}
	}
}
/*
Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        gl_Vertex = Viewport*Projection*ModelView*gl_Vertex;     // transform it to screen coordinates
        varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity
        return gl_Vertex;
    }
//*/

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

	//model = new Model("objet/diablo.obj");
	//model = new Model("objet/african_head.obj");
	//model = new Model("objet/y-wing.obj"); // ne marche pas
	//model = new Model("objet/imperial.obj"); // tres long et ne sais pas s il marche
	//load_texture("objet/african_head","_diffuse.tga",diffu);

	diffu.read_tga_file("objet/african_head_diffuse.tga");
	//diffu.get(0,0);

	//Face *face;
	//Model::diffuse(Vec2f(1,1));

	int nb = nombreFace("objet/african_head.obj");
	//int nb = nombreFace("objet/diablo.obj");
	//int nb = nombreFace("objet/y-wing.obj");
	//int nb = nombreFace("objet/imperial.obj");

	ModelView  = lookat(eye, center, Vec3f(0,1,0));
	Matrix Projection = Matrix::identity(4);
	Matrix ViewPort   = viewport(width/8, height/8, width*3/4, height*3/4);

	//diffu = new TGAImage("objet/african_head_diffuse.tga");

	srand(time(NULL));

    TGAImage image(width, height, TGAImage::RGB);
    //std::vector<int> buffer(width*height,-std::numeric_limits<int>::max());
	//printf("fin buff \n");
	
	float *buffer = new float[width*height];
    for (int i=width*height; i--; buffer[i] = -std::numeric_limits<float>::max());

	//Vec3f poule= Vec3f(ViewPort*Projection*ModelView*Matrix(light_src));
	//printf("%f",poule.x);

	for(int i=0; i<nb; i++){
		//float intensity = n*light_dir;
		//float intensity = 1;
		Face f = modele[i];
		Vec3f s1=sommet[f.s1-1];
		//printf("%f %f %f\n",s1.x,s1.y,s1.z);
		s1 = Vec3f(ViewPort*Projection*ModelView*Matrix(s1));
		//printf("%f %f %f\n",s1.x,s1.y,s1.z);
		Vec3f s2=sommet[f.s2-1];
		s2 = Vec3f(ViewPort*Projection*ModelView*Matrix(s2));

		Vec3f s3=sommet[f.s3-1];		
		s3 = Vec3f(ViewPort*Projection*ModelView*Matrix(s3));

		Triangle *tri = new Triangle(s1,s2,s3,f);//
		triangleBuffIntense(*tri,buffer, image);
		//triangleBuff(*tri,buffer, image);
		//triangle2(*tri,image);
	}

	//aff(light_src,image);

	image.flip_vertically();
    image.write_tga_file("dessin.tga");
    return 0;
}
