#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

typedef struct _Attributes{
    float vector_tangens;
    int positive;
}Attributes;

typedef enum _Mode{
    Seeking,Connected
}Mode;

typedef struct _Point{
    float x,y;
}Point,*pPoint;

bool does_circle_contain_point(Point center,Point point,int radius){
    return (center.x-point.x)*(center.x-point.x)+(center.y-point.y)*(center.y-point.y)<=radius*radius;
}

//https://www.desmos.com/calculator/xvaas2yyes
float len_sqr(Point center, Point point, float radius,Attributes attributes){
    float dy = point.y-center.y;
    float t = attributes.vector_tangens;
    float a = (t*t+1);
    float b = (2*t*t*point.x-2*t*dy+2*center.x);
    float c = (center.x*center.x-radius*radius+t*t*point.x*point.x-2*t*dy*point.x+dy*dy);
    float D = b*b-4*a*c;
    float x_coord = 0;
    if((t>=0)^(attributes.positive)){
        x_coord = (b-sqrt(D))/2/a;
    }
    else{
         x_coord = (b+sqrt(D))/2/a;
    }
    float dy_coord = t*(x_coord-point.x);
    return (x_coord-point.x)*(x_coord-point.x)+dy_coord*dy_coord;
}
typedef struct _ImportantPoint{
    Point coordinates;
    int sputnik_id;
    float relax_time;
    Mode mode;
}ImportantPoint,*pImportantPoint;

typedef struct _Sputnik{
    Point coordinates;
    float tan_alpha;
    float heigth;
}Sputnik,*pSputnik;

int main()
{
    size_t S,P,T;//A number of sputniks,points and requests
    cin>>S>>P>>T;

    //A list of coordinates of points
    vector <ImportantPoint> points(P);
    for(size_t i = 0; i < P; i++){
        cin>>points[i].coordinates.x
        >>points[i].coordinates.y;
    }

    //A list of alpha, height and coordinates of the sputnik
    vector <Sputnik> sputniks(S);
    for(size_t i = 0; i < S; i++){
        float next_alpha;
        cin>>next_alpha;
        sputniks[i].tan_alpha = tan(next_alpha);
    }

    //Read the first data
    float last_t;
    cin>>last_t;
    for(size_t i = 0; i < S; i++){
        cin>>sputniks[i].coordinates.x
        >>sputniks[i].coordinates.y
        >>sputniks[i].heigth;
    }
    for(size_t t = 0; t < T-1; t++){
        float current_time;
        cin>>current_time;
        float delta_time = current_time-last_t;
        for(size_t i = 0; i < S; i++){
            float next_x, next_y, next_height;
            cin>>next_x>>next_y>>next_height;
            float delta_heigth = sputniks[i].heigth/next_height;
            float speedX = (next_x-sputniks[i].coordinates.x)*delta_heigth/(delta_time);
            float speedY = (next_y-sputniks[i].coordinates.y)*delta_heigth/(delta_time);
            sputniks[i].coordinates.x = next_x;
            sputniks[i].coordinates.y = next_y;
            float radius = sputniks[i].heigth*sputniks[i].tan_alpha;
            for(size_t j = 0; j < P; j++){
                if(points[j].mode == Connected
                   || !does_circle_contain_point(sputniks[i].coordinates,points[j].coordinates,radius)){
                    continue;
                }
                float len = len_sqr(sputniks[i].coordinates,points[j].coordinates,radius,{speedX==0?(radius*radius):speedY/speedX,speedY>=0});
                float time = sqrt(len/(speedX*speedX+speedY*speedY));
                if(time>points[j].relax_time){
                    points[j].relax_time = time;
                    points[j].sputnik_id = i+1;
                }
            }
        }
        for(size_t i = 0; i < P; i++){
            if(points[i].sputnik_id > 0){
                points[i].mode = Connected;
            }
            cout<<points[i].sputnik_id<<" ";
        }
        cout<<endl;
        for(size_t i = 0; i < P; i++){
            points[i].relax_time -= delta_time;
            if(points[i].mode == Connected && points[i].relax_time <= 0){
                points[i].mode = Seeking;
                points[i].sputnik_id = 0;
                points[i].relax_time = 0;
            }
        }
    }

    return 0;
}

