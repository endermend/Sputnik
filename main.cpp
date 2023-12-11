#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#define EARTH_R 6376778.19784

using namespace std;

typedef struct _Point{
    double x,y;
}Point,*pPoint;

double distance (double a, double b){
    if (a > b){
        if(a > M_PI + b){
            return 2*M_PI+b-a;
        }
        return a-b;
    }
    if(b > M_PI + a){
        return 2*M_PI+a-b;
    }
    return b-a;
}

void calibrate(pPoint a,pPoint b){
    if (a->x > b->x){
        if(a->x > M_PI + b->x){
            b->x += 2*M_PI;
        }
    }
    else{
        if(b->x > M_PI + a->x){
            a->x += 2*M_PI;
        }
    }

    if (a->y > b->y){
        if(a->y > M_PI + b->y){
            b->y += 2*M_PI;
        }
    }
    else{
        if(b->y > M_PI + a->y){
            a->y += 2*M_PI;
        }
    }
}

double dist(double x, double y, double z){
    return sqrt(x*x+y*y+z*z);
}

double toangleX (double x, double y, double radius){
    return acos(x/radius)*(y>=0?1:-1);
}

double toangleY (double z, double x, double circle_rad, double sphere_rad){
    return  acos(circle_rad/sphere_rad*(x>=0?1:-1))*(z>=0?1:-1);
}

//https://www.desmos.com/calculator/bxkkezeyoy
double calc_radius(double dsq, double alpha){
    double cos = std::cos(alpha);
    double sin = std::sin(alpha);
    double l = (sqrt(dsq)*cos-sqrt(EARTH_R*EARTH_R-dsq*sin*sin))*sin;
    return std::asin(l/EARTH_R);
}

bool does_circle_contain_point(Point center,Point point,int radius){
    return (center.x-point.x)*(center.x-point.x)+(center.y-point.y)*(center.y-point.y)<=radius*radius;
}
//Idea: https://www.geogebra.org/calculator/zubjv67k
//Realization:https://www.desmos.com/calculator/xvaas2yyes
double len_sqr(Point center, Point point, Point vec, double radius){
    double dy = point.y-center.y;
    //tangens
    double t = vec.x==0?(radius*radius*radius):vec.y/vec.x;//if vec == 0 -> tangents = infinity, so i use radius^3
    //a,b,c: ax^2-bx+c=0
    double a = 2*(t*t+1);//multiply by 2 to save some time for calculation
    double b = 2*(t*t*point.x-t*dy+center.x);
    double c = (center.x*center.x-radius*radius+t*t*point.x*point.x-2*t*dy*point.x+dy*dy);
    double D = b*b-2*a*c;
    //if vector goes to right calculate right intersection, else - left
    double dx_coord = (b+(vec.x >=0?sqrt(D):-sqrt(D)))/a - point.x;
    double dy_coord = t*dx_coord;
    //return the square of len to save some time for calculation
    return dx_coord*dx_coord+dy_coord*dy_coord;
}

typedef struct _ImportantPoint{
    Point coordinates;
    int sputnik_id;
    double relax_time;
    bool connected;
}ImportantPoint,*pImportantPoint;

typedef struct _Sputnik{
    Point coordinates;
    double alpha;
    double heigth;

    double radius(){
        return calc_radius(heigth*heigth,alpha);
    }
}Sputnik,*pSputnik;

int main()
{
    ifstream cells("cells.txt");
    cells.tie(0);
    ofstream out("output.txt");
    size_t S,P,T;
    //A number of sputniks,points and requests
    //cin>>S>>P>>T;

    cells>>P;
    //A list of coordinates of points
    ImportantPoint *points = (ImportantPoint*)calloc(P,sizeof(ImportantPoint));
    for(size_t i = 0; i < P; i++){
        //TODO AHTUNG!
        double x,y,z;
        cells>>x>>y>>z;
        double heigth = sqrt(x*x+y*y+z*z);
        double radius = sqrt(x*x+y*y);
        points[i].coordinates.x = toangleX(x,y,radius);
        points[i].coordinates.y = toangleY(z,x,radius,heigth);
        //cout<<points[i].coordinates.x<<" "<<points[i].coordinates.y<<endl;
    }

    cells.close();
    //A list of alpha, height and coordinates of the sputnik
    ifstream sats("satellites.txt");
    Sputnik *sputniks = (Sputnik*)calloc(S,sizeof(Sputnik));
    sats.tie(0);
    sats>>S>>T;
    //Read the first data
    double last_t = 0;
    //cin>>last_t;
    for(size_t i = 0; i < S; i++){
        double x,y,z;
        int alpha;
        sats>>x>>y>>z>>alpha;
        sputniks[i].alpha = alpha;
        sputniks[i].heigth=dist(x,y,z);
        double radius = sqrt(x*x+y*y);
        //cout<<radius<<" "<<sputniks[i].heigth<<endl;
        sputniks[i].coordinates.x = toangleX(x,y,radius);
        sputniks[i].coordinates.y = toangleY(z,x,radius,sputniks[i].heigth);
    }

    //Time, Position x, Position y, Height above earth
    for(size_t t = 0; t < T-1; t++){
        cout<<t<<endl;
        double current_time = last_t + 1;
        //cin>>current_time;
        double delta_time = current_time-last_t;
        for(size_t i = 0; i < S; i++){
            double x,y,z;
            sats>>x>>y>>z>>sputniks[i].alpha;

            double next_height=dist(x,y,z);
            double radiusXY = sqrt(sputniks[i].coordinates.x*sputniks[i].coordinates.x+sputniks[i].coordinates.y*sputniks[i].coordinates.y);
            Point spherical = {toangleX(sputniks[i].coordinates.x,sputniks[i].coordinates.y,radiusXY),toangleY(z,x,radiusXY,sputniks[i].heigth)};
            double delta_heigth = sputniks[i].heigth/next_height;


            double speedX = distance(spherical.x,sputniks[i].coordinates.x)*delta_heigth/(delta_time);
            double speedY = distance(spherical.y,sputniks[i].coordinates.y)*delta_heigth/(delta_time);

            //cout<<sputniks[i].coordinates.x<<" "<<sputniks[i].coordinates.y<<endl;
            sputniks[i].heigth = next_height;
            double radius = sputniks[i].radius();
            for(size_t j = 0; j < P; j++){
                //if point is already connecte - skip this point
                if(points[j].connected){
                    continue;
                }
                calibrate(&sputniks[i].coordinates,&points[j].coordinates);
                //if circle don't contain a point - skip this point

                if(!does_circle_contain_point(sputniks[i].coordinates,points[j].coordinates,radius)){
                    continue;
                }
                /*
                cout<<sputniks[i].coordinates.x<<" "<<sputniks[i].coordinates.y<<" "
                <<points[j].coordinates.x<<" "<<points[j].coordinates.y<<" "<<radius<<endl;
                */
                //Calculate the length of section from the point to the intersection of opposite vector of sputnik's speed from the point and circle
                double len = len_sqr(sputniks[i].coordinates,points[j].coordinates,{-speedX,-speedY},radius);
                //time is ratio of length and speed
                double time = sqrt(len/(speedX*speedX+speedY*speedY));

                //search for maximum time
                if(time>points[j].relax_time){
                    points[j].relax_time = time;
                    points[j].sputnik_id = i+1;
                }
            }
            //Update sputnik's data
            last_t = current_time;
            sputniks[i].coordinates = spherical;
        }
        for(size_t i = 0; i < P; i++){
            if(points[i].sputnik_id > 0){
                points[i].connected = true;
            }
            out<<points[i].sputnik_id<<" ";

            points[i].relax_time -= delta_time;
            if(points[i].connected && points[i].relax_time <= 0){
                points[i].connected = 0;
                points[i].sputnik_id = 0;
                points[i].relax_time = 0;
            }
        }
        out<<endl;
    }
    for(size_t i = 0; i < P; i++){
        if(points[i].connected){
            goto print;
        }
        for(int j = 0; j < S; j++){
            calibrate(&sputniks[j].coordinates,&points[i].coordinates);
            if(does_circle_contain_point(sputniks[j].coordinates,points[i].coordinates,sputniks[j].radius())){
                points[i].sputnik_id = j+1;
                break;
            }
        }
        print:
            out<< points[i].sputnik_id<<" ";
    }
    return 0;
}

