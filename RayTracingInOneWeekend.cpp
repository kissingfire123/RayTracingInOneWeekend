#include <iostream>
#include <fstream>
#include <functional>
#include <io.h>
#include "include/vec3.h"
#include "include/ray.h"
#include "include/sphere.h"
#include "include/hitable_list.h"
#include "include/camera.h"

const static int g_Width = 800;
const static int g_Height = 400;
const static int g_MAX_FLOAT = 1000.0;
const static int g_RayNums = 100;

int Ch1OutputImage(){
    std::ofstream imageFile("Image01.ppm");
    
    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = g_Width -1; i >= 0; --i){
            float r = float(i) / float(g_Width);
            float g = float(j) / float(g_Height);
            float b = 0.3;
            int ir = int(255.99 * r);
            int ig = int(255.99 * g);
            int ib = int(255.99 * b);
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
    }
    imageFile.close();
    return 0;
}

//the vec3 class
int Ch2OutputImage(){
    std::ofstream imageFile("Image02.ppm");

    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = g_Width -1; i >= 0; --i){
            vec3 color(float(i) / float(g_Width),float(j) / float(g_Height),0.3);
            int ir = int(255.99 * color.r());
            int ig = int(255.99 * color.g());
            int ib = int(255.99 * color.b());
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
    }

    imageFile.close();
    return 0;
}

//rays,a simple camera,and background
// P = Origin + t * Direction
// coordinate: top:y+  , right: x+  , outPcScreen: z+
int Ch3SimpleCamImage(){
    auto getColor = [](const ray&r) -> vec3{
      vec3 unit_direction = unit_vector(r.direction()); 
      float t = 0.5 * (unit_direction.y() + 1.0);
      return (1.0 - t) * vec3(1.0,1.0,1.0) + t *vec3(0.5,0.7,1.0);  
    }; 

    std::ofstream imageFile("Image03.ppm");
    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    vec3 lower_left_corner_P(-2.0,-1.0,-1.0);
    vec3 horizontalDir(4.0,0.0,0.0);
    vec3 verticalDir(0.0,2.0,0.0);
    vec3 originP(0.0,0.0,0.0);

    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = g_Width -1; i >= 0; --i){
            float u = float(i) / float(g_Width);
            float v = float(j) / float(g_Height);
            ray r(originP,lower_left_corner_P + u*horizontalDir + v*verticalDir);
            
            vec3 color = getColor(r);
            int ir = int(255.99 * color.r());
            int ig = int(255.99 * color.g());
            int ib = int(255.99 * color.b());
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
    }
    imageFile.close();
    return 0; 
}

// 1. point_P = pointA + t * direction_B
// 2. point_P , circleCenter_C,vector_length(p-c) = radius_R
//    dot((A​​ + t*​B ​- ​C​),(​A​ + t*​B​ - ​C​)) = R*R
// 3. t*t*dot(B​ ​,​B​) + 2*t*dot(​B,A​-​C​) + dot(​A-C,A​-​C​) - R*R = 0
int Ch4AddSphere(){
    //resolve: either hit a sphere,
    auto hit_sphere = [](const vec3& center, float radius,const ray& r)-> bool{
        vec3 ac = r.origin() - center; // vector(A-C)
        float a = dot(r.direction() , r.direction());
        float b = 2.0 * dot(ac, r.direction());
        float c = dot(ac,ac) - radius*radius;
        float deltaDiscriminant = b * b - 4 * a * c;
        return deltaDiscriminant >= 0.0001f;
    };


    auto getColor = [&](const vec3& center, float radius,const ray&r) -> vec3{
      if(hit_sphere(center,radius,r)){  
          return vec3(1.0,0.0,0.0);
      }
      vec3 unit_direction = unit_vector(r.direction()); 
      float t = 0.5 * (unit_direction.y() + 1.0);
      return (1.0 - t) * vec3(1.0,1.0,1.0) + t *vec3(0.5,0.7,1.0);  
    }; 

    std::ofstream imageFile("Image04_add_sphere.ppm");
    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    vec3 lower_left_corner_P(-2.0,-1.0,-1.0);
    vec3 horizontalDir(4.0,0.0,0.0);
    vec3 verticalDir(0.0,2.0,0.0);
    vec3 originP(0.0,0.0,0.0);

    vec3 circleCenter (0.0,0.0,-1.0);
    float radius = 0.5;
    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = g_Width -1; i >= 0; --i){
            float u = float(i) / float(g_Width);
            float v = float(j) / float(g_Height);
            ray r(originP,lower_left_corner_P + u*horizontalDir + v*verticalDir);

            vec3 color = getColor(circleCenter,radius,r);
            int ir = int(255.99 * color.r());
            int ig = int(255.99 * color.g());
            int ib = int(255.99 * color.b());
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
    }
    imageFile.close();

    return 0;
}

// Ch5 Suface normals and multiple objects
int Ch5NormalsAndMultipleObj(){
    //resolve: either hit a sphere,get delta
    auto hit_sphere = [](const vec3& center, float radius,const ray& r)-> float{
        vec3 ac = r.origin() - center; // vector(A-C)
        float a = dot(r.direction() , r.direction());
        float b = 2.0 * dot(ac, r.direction());
        float c = dot(ac,ac) - radius*radius;
        float deltaDiscriminant = b * b - 4 * a * c;
        if(deltaDiscriminant <= 0.0001f){
            return -1.0;
        }
        else{
            return (-b - sqrt(deltaDiscriminant))/(2.0*a);
        }
    };

    auto getColor = [&](const vec3& center, float radius,const ray&r) -> vec3{
      float t =  hit_sphere(center,radius,r);
      if(t > 0.0001f){  
          vec3 N = unit_vector(r.at_Parameter(t) - center);
          return 0.5*vec3(N.x() + 1,N.y() + 1, N.z() + 1);
      }
      vec3 unit_direction = unit_vector(r.direction()); 
      t = 0.5 * (unit_direction.y() + 1.0);
      return (1.0 - t) * vec3(1.0,1.0,1.0) + t *vec3(0.5,0.7,1.0);  
    }; 


    std::ofstream imageFile("Image05_normals.ppm");
    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    vec3 lower_left_corner_P(-2.0,-1.0,-1.0);
    vec3 horizontalDir(4.0,0.0,0.0);
    vec3 verticalDir(0.0,2.0,0.0);
    vec3 originP(0.0,0.0,0.0);

    vec3 circleCenter (0.0,0.0,-1.0);
    float radius = 0.5;
    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = g_Width -1; i >= 0; --i){
            float u = float(i) / float(g_Width);
            float v = float(j) / float(g_Height);
            ray r(originP,lower_left_corner_P + u*horizontalDir + v*verticalDir);

            vec3 color = getColor(circleCenter,radius,r);
            int ir = int(255.99 * color.r());
            int ig = int(255.99 * color.g());
            int ib = int(255.99 * color.b());
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
    }
    imageFile.close();
    
    return 0;
}

//Ch5: multi-object and ray
int Ch5MultiObjHitableWith_tRange(){
    auto getColor = [&](const ray&r,hitable *world) -> vec3{
        hit_record reco;
        //根据光线击中的最近点，进行渲染着色
        if(world -> hit(r,0.0,g_MAX_FLOAT,reco)){
            return 0.5*vec3(reco.normal_.x() + 1,reco.normal_.y() + 1, reco.normal_.z() + 1);
        }
        else{
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5 * (unit_direction.y() + 1.0);
            return (1.0 - t) * vec3(1.0,1.0,1.0) + t *vec3(0.5,0.7,1.0);  
        } 
    };
    
    std::string imgFilePath("Image05_with_tRange.ppm");
    if(access(imgFilePath.c_str(),0) == 0){
        std::remove(imgFilePath.c_str());
    }
    std::ofstream imageFile(imgFilePath);
    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    vec3 lower_left_corner_P(-2.0,-1.0,-1.0);
    vec3 horizontalDir(4.0,0.0,0.0);
    vec3 verticalDir(0.0,2.0,0.0);
    vec3 originP(0.0,0.0,0.0);

// multi-object
    std::shared_ptr<hitable> list[2];
    list[0] = std::make_shared<sphere>(vec3(0,0,-1),0.5);
    list[1] = std::make_shared<sphere>(vec3(0,-60.5,-1),60); 
    std::shared_ptr<hitable> world = std::make_shared<hitable_list>(list,2);

    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = g_Width -1; i >= 0; --i){
            float u = float(i) / float(g_Width);
            float v = float(j) / float(g_Height);
            ray r(originP,lower_left_corner_P + u*horizontalDir + v*verticalDir);

            vec3 color = getColor(r,world.get());
            int ir = int(255.99 * color.r());
            int ig = int(255.99 * color.g());
            int ib = int(255.99 * color.b());
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
    }
    imageFile.close();
    std::cout<< imgFilePath;
    return 0;
}


//Ch6: Antialiasing 
int Ch6Antialiasing(){
    auto getColor = [&](const ray&r,hitable *world) -> vec3{
        hit_record reco;
        //根据光线击中的最近点，进行渲染着色
        if(world -> hit(r,0.0,g_MAX_FLOAT,reco)){
            return 0.5*vec3(reco.normal_.x() + 1,reco.normal_.y() + 1, reco.normal_.z() + 1);
        }
        else{
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5 * (unit_direction.y() + 1.0);
            return (1.0 - t) * vec3(1.0,1.0,1.0) + t *vec3(0.5,0.7,1.0);  
        } 
    };
    
    std::string imgFilePath("Image06_AntiAliasing.ppm");
    if(access(imgFilePath.c_str(),0) == 0){
        std::remove(imgFilePath.c_str());
    }
    std::ofstream imageFile(imgFilePath);
    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    vec3 lower_left_corner_P(-2.0,-1.0,-1.0);
    vec3 horizontalDir(4.0,0.0,0.0);
    vec3 verticalDir(0.0,2.0,0.0);
    vec3 originP(0.0,0.0,0.0);

// multi-object
    std::shared_ptr<hitable> list[2];
    list[0] = std::make_shared<sphere>(vec3(0,0,-1),0.5);
    list[1] = std::make_shared<sphere>(vec3(0,-60.5,-1),60); 
    std::shared_ptr<hitable> world = std::make_shared<hitable_list>(list,2);
    camera cam;//多条光线打向同一个pixel，模拟MSAA进行抗混叠
    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = g_Width -1; i >= 0; --i){
            vec3 color(0,0,0);
            for(int s = 0; s< g_RayNums ; ++s){
                float u = float(i + random_double())/ float(g_Width);
                float v = float(j + random_double())/ float(g_Height);
                ray r = cam.get_ray(u,v);
                color += getColor(r,world.get());
            }
            color /= float(g_RayNums);
            int ir = int(255.99 * color.r());
            int ig = int(255.99 * color.g());
            int ib = int(255.99 * color.b());
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
    }

    imageFile.close();
    std::cout<< imgFilePath;
    return 0;
}

int Ch7DiffuseMaterial(){
    auto random_in_unit_sphere = [&](){
        vec3 p;
        do{
            p = 2.0 * vec3(random_double(),random_double(),random_double()) - vec3(1,1,1);
        }while (p.length_squared() >= 1.0);
        return p;
    };
    using getColorFuncType = std::function<vec3(const ray&r,hitable *world)>;
    getColorFuncType getColor = [&](const ray&r,hitable *world) -> vec3{
        hit_record reco;
        //根据光线击中的最近点，进行渲染着色
        if(world -> hit(r,0.0,g_MAX_FLOAT,reco)){
            vec3 target = reco.p_ + reco.normal_ + random_in_unit_sphere();
            return 0.5* getColor(ray(reco.p_,target-reco.p_),world);
        }
        else{
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5 * (unit_direction.y() + 1.0);
            return (1.0 - t) * vec3(1.0,1.0,1.0) + t *vec3(0.5,0.7,1.0);  
        } 
    };
    

    std::string imgFilePath("Image07_DiffuseMaterial.ppm");
    if(access(imgFilePath.c_str(),0) == 0){
        std::remove(imgFilePath.c_str());
    }
    std::ofstream imageFile(imgFilePath);
    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    vec3 lower_left_corner_P(-2.0,-1.0,-1.0);
    vec3 horizontalDir(4.0,0.0,0.0);
    vec3 verticalDir(0.0,2.0,0.0);
    vec3 originP(0.0,0.0,0.0);

// multi-object
    std::shared_ptr<hitable> list[2];
    list[0] = std::make_shared<sphere>(vec3(0,0,-1),0.5);
    list[1] = std::make_shared<sphere>(vec3(0,-60.5,-1),60); 
    std::shared_ptr<hitable> world = std::make_shared<hitable_list>(list,2);
    camera cam;//多条光线打向同一个pixel，模拟MSAA进行抗混叠
    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = g_Width -1; i >= 0; --i){
            vec3 color(0,0,0);
            for(int s = 0; s< g_RayNums ; ++s){
                float u = float(i + random_double())/ float(g_Width);
                float v = float(j + random_double())/ float(g_Height);
                ray r = cam.get_ray(u,v);
                color += getColor(r,world.get());
            }
            color /= float(g_RayNums);
            int ir = int(255.99 * color.r());
            int ig = int(255.99 * color.g());
            int ib = int(255.99 * color.b());
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
    }
    imageFile.close();
    std::cout<< imgFilePath;   
    return 0;
}

int main(int argc,char* argv[]){
    std::cout << "OutputImage: \n" ;
    //Ch1OutputImage();
    //Ch2OutputImage();
    //Ch3SimpleCamImage();
    //Ch4AddSphere();
    //Ch5NormalsAndMultipleObj();
    //Ch5MultiObjHitableWith_tRange();
    //Ch6Antialiasing();
    Ch7DiffuseMaterial();
    return 0;
}