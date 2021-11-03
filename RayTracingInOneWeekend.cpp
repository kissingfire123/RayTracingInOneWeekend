#include <iostream>
#include <fstream>
#include <functional>
#include <io.h>
#include <chrono>

#include "progress.h"
#include "vec3.h"
#include "ray.h"
#include "sphere.h"
#include "hitable_list.h"
#include "camera.h"
#include "material.h"

//Declaration
int Ch1OutputImage(std::string imgFilePath);
int Ch2OutputImage(std::string imgFilePath);
int Ch3SimpleCamImage(std::string imgFilePath);
int Ch4AddSphere(std::string imgFilePath);
int Ch5NormalsAndMultipleObj(std::string imgFilePath);
int Ch5MultiObjHitableWith_tRange(std::string imgFilePath);
int Ch6Antialiasing(std::string imgFilePath);
int Ch7DiffuseMaterial(std::string imgFilePath);
int Ch8MaterialMetal(std::string imgFilePath);


const static int g_Width = 800;
const static int g_Height = 400;
const static int g_MAX_FLOAT = 1000;
const static int g_RayNums = 100;
const static int g_DepthThreshold = 50;


int main(int argc, char* argv[]) { 
#ifdef _DEBUG
    std::cout << "Now running model is Debug\n";  // debug enable
#else
    std::cout << "Now running model is Release\n";// more faster
#endif
    //以下每个Chx...函数都可以单独运行,互不影响,可以屏蔽其中任意几个单独运行其他的 

    Ch1OutputImage("Image01.ppm");
    Ch2OutputImage("Image02.ppm");
    Ch3SimpleCamImage("Image03.ppm");
    Ch4AddSphere("Image04_add_sphere.ppm");
    Ch5NormalsAndMultipleObj("Image05_normals.ppm");
    Ch5MultiObjHitableWith_tRange("Image05_with_tRange.ppm");
    Ch6Antialiasing("Image06_AntiAliasing.ppm");
    Ch7DiffuseMaterial("Image07_DiffuseMaterial.ppm");
    Ch8MaterialMetal("Image08_MetalMaterial.ppm");

    return 0;
}

int Ch1OutputImage(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath,g_Height);
    std::ofstream imageFile(imgFilePath);
    
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
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    return 0;
}

//the vec3 class
int Ch2OutputImage(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath, g_Height);
    std::ofstream imageFile(imgFilePath);

    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = g_Width -1; i >= 0; --i){
            vec3 color(float(i) / float(g_Width),float(j) / float(g_Height),0.3);
            int ir = int(255.99 * color.r());
            int ig = int(255.99 * color.g());
            int ib = int(255.99 * color.b());
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
        rtwProgress.Refresh(g_Height - j);
    }

    imageFile.close();
    return 0;
}

//rays,a simple camera,and background
// P = Origin + t * Direction
// coordinate: top:y+  , right: x+  , outPcScreen: z+
int Ch3SimpleCamImage(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath, g_Height);
    auto getColor = [](const ray&r) -> vec3{
      vec3 unit_direction = unit_vector(r.direction()); 
      float t = 0.5 * (unit_direction.y() + 1.0);
      return (1.0 - t) * vec3(1.0,1.0,1.0) + t *vec3(0.5,0.7,1.0);  
    }; 

    std::ofstream imageFile(imgFilePath);
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
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    return 0; 
}

// 1. point_P = pointA + t * direction_B
// 2. point_P , circleCenter_C,vector_length(p-c) = radius_R
//    dot((A​​ + t*​B ​- ​C​),(​A​ + t*​B​ - ​C​)) = R*R
// 3. t*t*dot(B​ ​,​B​) + 2*t*dot(​B,A​-​C​) + dot(​A-C,A​-​C​) - R*R = 0
int Ch4AddSphere(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath, g_Height);
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

    std::ofstream imageFile(imgFilePath);
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
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    return 0;
}

// Ch5 Suface normals and multiple objects
int Ch5NormalsAndMultipleObj(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath, g_Height);
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


    std::ofstream imageFile(imgFilePath);
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
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    return 0;
}

//Ch5: multi-object and ray
int Ch5MultiObjHitableWith_tRange(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath, g_Height);
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
    list[1] = std::make_shared<sphere>(vec3(0,-60.5,-1),60.0);
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
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    return 0;
}


//Ch6: Antialiasing 
int Ch6Antialiasing(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath, g_Height);
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
    list[1] = std::make_shared<sphere>(vec3(0,-60.5,-1),60.0);
    std::shared_ptr<hitable> world = std::make_shared<hitable_list>(list,2);
    camera cam;//Ch6: 多条光线打向同一个pixel，模拟MSAA进行抗混叠
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
        rtwProgress.Refresh(g_Height - j);
    }

    imageFile.close();
    return 0;
}

//Ch7: 模拟杂乱无章随机的漫反射
int Ch7DiffuseMaterial(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath, g_Height);
    // Ch7:取单位半径球内的任意一点vec3(x,y,z),用作表现反射的随机性
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
        //Ch6:根据光线击中的最近点，进行渲染着色
        if(world -> hit(r,0.001,g_MAX_FLOAT,reco)){//Ch7: 0.001f,表示去除靠近0的浮点值，避免浮点精度带来的毛刺
            vec3 target = reco.p_ + reco.normal_ + random_in_unit_sphere();//Ch7:p_+normal_得到hit-point的球心，再随机选个方向作为反射关系
            return 0.5* getColor(ray(reco.p_,target-reco.p_),world);//Ch7:递归调用，即多次反射，直到hit-miss
        }
        else{
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5 * (unit_direction.y() + 1.0);
            return (1.0 - t) * vec3(1.0,1.0,1.0) + t *vec3(0.5,0.7,1.0);  
        } 
    };
    

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
    list[1] = std::make_shared<sphere>(vec3(0,-60.5,-1),60.0);
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
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    return 0;
}


int Ch8MaterialMetal(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath, g_Height);
   //Ch8: modify getColor()
   using getColorFuncType = std::function<vec3(const ray&r,hitable *world,int depth)>;
    getColorFuncType getColor = [&](const ray&r,hitable *world,int depth) -> vec3{
        hit_record reco;
        //Ch6:根据光线击中的最近点，进行渲染着色
        if(world -> hit(r,0.001,g_MAX_FLOAT,reco)){//Ch7: 0.001f,表示去除靠近0的浮点值，避免浮点精度带来的毛刺
            ray scattered;
            vec3 attenuation;//Ch8 : 材料属性,反射率,吸光率
            if(depth < g_DepthThreshold && reco.mate_ptr->scatter(r,reco,attenuation,scattered)){
                return attenuation *getColor(scattered,world,depth + 1);//递归,继续反射
            }
            else{
                return vec3(0,0,0);
            }
        }
        else{
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5 * (unit_direction.y() + 1.0);
            return (1.0 - t) * vec3(1.0,1.0,1.0) + t *vec3(0.5,0.7,1.0);  
        } 
    };
    

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
    std::shared_ptr<hitable> list[4];
    list[0] = std::make_shared<sphere>(vec3(0,0,-1),       0.5,  std::make_shared<lambertian>(vec3(0.8,0.3,0.3)));
    list[1] = std::make_shared<sphere>(vec3(0,-100.5,-1),100.0, std::make_shared<lambertian>(vec3(0.8, 0.8, 0.3)));
    list[2] = std::make_shared<sphere>(vec3(1,0,-1),       0.4,  std::make_shared<metal>(vec3(0.8, 0.6, 0.2)));
    list[3] = std::make_shared<sphere>(vec3(-1,0,-1),      0.5,  std::make_shared<metal>(vec3(0.8,0.8,0.8))); 
    std::shared_ptr<hitable> world = std::make_shared<hitable_list>(list,4);
    camera cam;//多条光线打向同一个pixel，模拟MSAA进行抗混叠
    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = g_Width -1; i >= 0; --i){
            vec3 color(0,0,0);
            for(int s = 0; s< g_RayNums ; ++s){
                float u = float(i + random_double())/ float(g_Width);
                float v = float(j + random_double())/ float(g_Height);
                ray r = cam.get_ray(u,v);
                color += getColor(r,world.get(),0);
            }
            color /= float(g_RayNums);
            int ir = int(255.99 * color.r());
            int ig = int(255.99 * color.g());
            int ib = int(255.99 * color.b());
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    return 0;
}
