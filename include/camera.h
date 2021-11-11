#ifndef _CAMERA_H_
#define _CAMERA_H_

#include "ray.h"

class camera{
public:
    explicit camera(){}

    // verticalFov:角度制
    explicit camera(float verticalFov, float aspect,vec3 lookFrom,vec3 lookAt,vec3 lookUp)
        :origin_(lookFrom){
        float theta = verticalFov * RTW_PI / 180;
        float half_height = tan(theta/2);
        float half_width = aspect * half_height;

        //camera的uvw坐标系,w相当于之前的z,且-z指向camera
        vec3 u, v, w;
        w = unit_vector(lookFrom - lookAt);
        u = unit_vector(cross(lookUp,w));
        v = cross(w,u);

        //lower_left_corner_ = vec3(-half_width, -half_height, -1.0);
        lower_left_corner_ =  origin_ - half_width*u -half_height*v -w;
        horizontal_ = 2 * half_width * u;
        vertical_   = 2 * half_height* v;
    }

    // for antiliasing,move the origin,simu multi-rays for same one pixel 
    ray get_ray(float u, float v){
        return ray(origin_,lower_left_corner_ + u*horizontal_ + v*vertical_ - origin_);
    }

private:
    vec3 origin_ = vec3(0.0, 0.0, 0.0);
    vec3 lower_left_corner_ = vec3(-2.0, -1.0 , -1.0);
    vec3 horizontal_ = vec3(4.0 , 0.0, 0.0);
    vec3 vertical_ = vec3(0.0,2.0,0.0);
};



#endif