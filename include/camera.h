#ifndef _CAMERA_H_
#define _CAMERA_H_

#include "ray.h"

class camera{
public:
    explicit camera(){}

    // for antiliasing,move the origin,simu multi-rays for same one pixel 
    ray get_ray(float u, float v){
        return ray(origin_,lower_left_corner_ + u*horizontal_ + v*vertical_ - origin_);
    }

private:
    vec3 origin_ = vec3(0,0,0);
    vec3 lower_left_corner_ = vec3(-2.0, -1.0 , -1.0);
    vec3 horizontal_ = vec3(4.0 , 0.0, 0.0);
    vec3 vertical_ = vec3(0.0,2.0,0.0);
};



#endif