#include "tgaimage.h"
const TGAColor white = TGAColor(255, 255, 255, 255);

const TGAColor rose = TGAColor(255, 0, 255, 255);

const TGAColor bleu = TGAColor(0, 0, 255, 255);

void line(TGAImage &img){

    for(int i = 0 ; i < 100;i++){
        img.set(i,i,white);
       // img.set(i,i+1,white);
       // img.set(i,i+2,white);
    }



}

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

void fond(TGAImage &image){

    for(int i = 0 ; i < image.get_width();i++){
        for(int z = 0 ; z < image.get_height();z++){
             image.set(i,z,rose);
        }
    }

}

int main(int argc, char** argv) {
    TGAImage image(800, 600, TGAImage::RGB);
    fond(image);
    line(image,white,130,210,420,100);
    line2(image,bleu,0,130,100,100);
    line2(image,bleu,400, 0, 800, 600);


    image.write_tga_file("dessin.tga");
    return 0;
}
