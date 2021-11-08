#include <iostream>
#include <fstream>
#include <functional>
#if    RTW_OS_WIN //通过CMake自行定义,区分系统类型
#include <io.h>
#elif  RTW_OS_MAC
#include <unistd.h>
#endif
#include <chrono>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"//vscode插件(PPM/PGM Viewer)可打开*.ppm,但不如*.bmp方便,都输出
#include "progress.h"


#include "vec3.h"
#include "ray.h"
#include "sphere.h"
#include "hitable_list.h"
#include "camera.h"
#include "material.h"

//Declaration
int Ch1_OutputImage(std::string imgFilePath);
int Ch2_OutputImage(std::string imgFilePath);
int Ch3_SimpleCamImage(std::string imgFilePath);
int Ch4_AddSphere(std::string imgFilePath);
int Ch5_NormalsAndMultipleObj(std::string imgFilePath);
int Ch5_MultiObjHitableWith_tRange(std::string imgFilePath);
int Ch6_Antialiasing(std::string imgFilePath);
int Ch7_DiffuseMaterial(std::string imgFilePath);
int Ch8_MaterialMetal(std::string imgFilePath);
int Ch9_Dielectrics(std::string imgFilePath);


const static int g_Width  = 800;//长宽比 2:1,和拟定的像素坐标比例一致
const static int g_Height = 400;
const static double g_MAX_TmFloat = 10000;//std::numeric_limits<double>::infinity();
const static int g_RayNums = 100;
const static int g_DepthThreshold = 50;


int main(int argc, char* argv[]) { 
#if  RTW_OS_WIN //have no idea for Xcode
#ifdef _DEBUG
    std::cout << "Now running model is Debug\n";  // debug enable
#else
    std::cout << "Now running model is Release\n";// more faster
#endif
#endif
    //以下每个Chx...函数都可以单独运行,互不影响,可以屏蔽其中任意几个单独运行其他的 

    Ch1_OutputImage("./Image01.ppm");
    Ch2_OutputImage("./Image02.ppm");
    Ch3_SimpleCamImage("./Image03.ppm");
    Ch4_AddSphere("./Image04_add_sphere.ppm");
    Ch5_NormalsAndMultipleObj("./Image05_normals.ppm");
    Ch5_MultiObjHitableWith_tRange("./Image05_with_tRange.ppm");
    Ch6_Antialiasing("./Image06_AntiAliasing.ppm");
    Ch7_DiffuseMaterial("./Image07_DiffuseMaterial.ppm");
    Ch8_MaterialMetal("./Image08_MetalMaterial.ppm");
    Ch9_Dielectrics("./Image09_DilectricsMaterial.ppm");
    return 0;
}

//Ch1: 对每个pixel逐次赋值Rgb,生成一张图片
int Ch1_OutputImage(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath,g_Height);
    std::ofstream imageFile(imgFilePath);
    std::vector<unsigned char> imgData;
    
    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = 0 ; i < g_Width; ++i){
            float r = float(i) / float(g_Width);
            float g = float(j) / float(g_Height);
            float b = 0.3;
            int ir = int(255.99 * r);  imgData.push_back(ir);
            int ig = int(255.99 * g);  imgData.push_back(ig);
            int ib = int(255.99 * b);  imgData.push_back(ib);
            imageFile << ir << " " << ig << " " << ib << "\n";
        }
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    imgFilePath.replace(imgFilePath.find(".ppm"),4,".bmp");
    stbi_write_bmp(imgFilePath.c_str(),g_Width,g_Height,3,imgData.data());
    return 0;
}

// Ch2: 在Ch1的基础上,加入vec3 class
int Ch2_OutputImage(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath, g_Height);
    std::ofstream imageFile(imgFilePath);
    std::vector<unsigned char> imgData;

    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = 0 ; i < g_Width; ++i){
            vec3 color(float(i) / float(g_Width),float(j) / float(g_Height),0.3);
            int ir = int(255.99 * color.r());   imgData.push_back(ir);
            int ig = int(255.99 * color.g());   imgData.push_back(ig);
            int ib = int(255.99 * color.b());   imgData.push_back(ib);
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
        rtwProgress.Refresh(g_Height - j);
    }

    imageFile.close();
    imgFilePath.replace(imgFilePath.find(".ppm"), 4, ".bmp");
    stbi_write_bmp(imgFilePath.c_str(), g_Width, g_Height, 3, imgData.data());
    return 0;
}

// Ch3： 根据向量原理，引入光线的概念，虚拟一台照相机观察光打在一个有限平面
//rays,a simple camera,and background
// P = Origin + t * Direction
// coordinate: top:y+  , right: x+  , outPcScreen: z+
int Ch3_SimpleCamImage(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath, g_Height);
    auto getColor = [](const ray&r) -> vec3{
      vec3 unit_direction = unit_vector(r.direction()); 
      float t = 0.5 * (unit_direction.y() + 1.0);
      return (1.0 - t) * vec3(1.0,1.0,1.0) + t *vec3(0.5,0.7,1.0);  
    }; 

    std::vector<unsigned char> imgData;
    std::ofstream imageFile(imgFilePath);
    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    vec3 lower_left_corner_P(-2.0,-1.0,-1.0);
    vec3 horizontalDir(4.0,0.0,0.0);
    vec3 verticalDir(0.0,2.0,0.0);
    vec3 originP(0.0,0.0,0.0);

    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = 0 ; i < g_Width; ++i){
            float u = float(i) / float(g_Width);
            float v = float(j) / float(g_Height);
            ray r(originP,lower_left_corner_P + u*horizontalDir + v*verticalDir);
            
            vec3 color = getColor(r);
            int ir = int(255.99 * color.r());  imgData.push_back(ir);
            int ig = int(255.99 * color.g());  imgData.push_back(ig);
            int ib = int(255.99 * color.b());  imgData.push_back(ib);
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    imgFilePath.replace(imgFilePath.find(".ppm"), 4, ".bmp");
    stbi_write_bmp(imgFilePath.c_str(), g_Width, g_Height, 3, imgData.data());
    return 0; 
}

// 1. point_P = pointA + t * direction_B
// 2. point_P , circleCenter_C,vector_length(p-c) = radius_R
//    dot((A​​ + t*​B ​- ​C​),(​A​ + t*​B​ - ​C​)) = R*R
// 3. t*t*dot(B​ ​,​B​) + 2*t*dot(​B,A​-​C​) + dot(​A-C,A​-​C​) - R*R = 0
int Ch4_AddSphere(std::string imgFilePath){
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

    std::vector<unsigned char> imgData;
    std::ofstream imageFile(imgFilePath);
    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    vec3 lower_left_corner_P(-2.0,-1.0,-1.0);
    vec3 horizontalDir(4.0,0.0,0.0);
    vec3 verticalDir(0.0,2.0,0.0);
    vec3 originP(0.0,0.0,0.0);

    vec3 circleCenter (0.0,0.0,-1.0);
    float radius = 0.5;
    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = 0 ; i < g_Width; ++i){
            float u = float(i) / float(g_Width);
            float v = float(j) / float(g_Height);
            ray r(originP,lower_left_corner_P + u*horizontalDir + v*verticalDir);

            vec3 color = getColor(circleCenter,radius,r);
            int ir = int(255.99 * color.r());  imgData.push_back(ir);
            int ig = int(255.99 * color.g());  imgData.push_back(ig);
            int ib = int(255.99 * color.b());  imgData.push_back(ib);
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    imgFilePath.replace(imgFilePath.find(".ppm"), 4, ".bmp");
    stbi_write_bmp(imgFilePath.c_str(), g_Width, g_Height, 3, imgData.data());
    return 0;
}

// Ch5 Suface normals and multiple objects
int Ch5_NormalsAndMultipleObj(std::string imgFilePath){
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

    std::vector<unsigned char> imgData;
    std::ofstream imageFile(imgFilePath);
    imageFile << "P3\n" << g_Width << " "  << g_Height << "\n255\n";
    vec3 lower_left_corner_P(-2.0,-1.0,-1.0);
    vec3 horizontalDir(4.0,0.0,0.0);
    vec3 verticalDir(0.0,2.0,0.0);
    vec3 originP(0.0,0.0,0.0);

    vec3 circleCenter (0.0,0.0,-1.0);
    float radius = 0.5;
    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = 0 ; i < g_Width; ++i){
            float u = float(i) / float(g_Width);
            float v = float(j) / float(g_Height);
            ray r(originP,lower_left_corner_P + u*horizontalDir + v*verticalDir);

            vec3 color = getColor(circleCenter,radius,r);
            int ir = int(255.99 * color.r());   imgData.push_back(ir);
            int ig = int(255.99 * color.g());   imgData.push_back(ig);
            int ib = int(255.99 * color.b());   imgData.push_back(ib);
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    imgFilePath.replace(imgFilePath.find(".ppm"), 4, ".bmp");
    stbi_write_bmp(imgFilePath.c_str(), g_Width, g_Height, 3, imgData.data());
    return 0;
}

//Ch5: multi-object and ray
int Ch5_MultiObjHitableWith_tRange(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath, g_Height);
    auto getColor = [&](const ray&r,hitable *world) -> vec3{
        hit_record reco;
        //根据光线击中的最近点，进行渲染着色
        if(world -> hit(r,0.0,g_MAX_TmFloat,reco)){
            return 0.5*vec3(reco.normal_.x() + 1,reco.normal_.y() + 1, reco.normal_.z() + 1);
        }
        else{
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5 * (unit_direction.y() + 1.0);
            return (1.0 - t) * vec3(1.0,1.0,1.0) + t *vec3(0.5,0.7,1.0);  
        } 
    };
    
    std::vector<unsigned char> imgData;
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
        for(int i = 0 ; i < g_Width; ++i){
            float u = float(i) / float(g_Width);
            float v = float(j) / float(g_Height);
            ray r(originP,lower_left_corner_P + u*horizontalDir + v*verticalDir);

            vec3 color = getColor(r,world.get());
            int ir = int(255.99 * color.r());   imgData.push_back(ir);
            int ig = int(255.99 * color.g());   imgData.push_back(ig);
            int ib = int(255.99 * color.b());   imgData.push_back(ib);
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    imgFilePath.replace(imgFilePath.find(".ppm"), 4, ".bmp");
    stbi_write_bmp(imgFilePath.c_str(), g_Width, g_Height, 3, imgData.data());
    return 0;
}


//Ch6: Antialiasing 
int Ch6_Antialiasing(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath, g_Height);
    auto getColor = [&](const ray&r,hitable *world) -> vec3{
        hit_record reco;
        //根据光线击中的最近点，进行渲染着色
        if(world -> hit(r,0.0,g_MAX_TmFloat,reco)){
            return 0.5*vec3(reco.normal_.x() + 1,reco.normal_.y() + 1, reco.normal_.z() + 1);
        }
        else{
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5 * (unit_direction.y() + 1.0);
            return (1.0 - t) * vec3(1.0,1.0,1.0) + t *vec3(0.5,0.7,1.0);  
        } 
    };
    
    std::vector<unsigned char> imgData;
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
        for(int i = 0 ; i < g_Width; ++i){
            vec3 color(0,0,0);
            for(int s = 0; s< g_RayNums ; ++s){
                float u = float(i + random_double())/ float(g_Width);
                float v = float(j + random_double())/ float(g_Height);
                ray r = cam.get_ray(u,v);
                color += getColor(r,world.get());
            }
            color /= float(g_RayNums);
            int ir = int(255.99 * color.r());  imgData.push_back(ir);
            int ig = int(255.99 * color.g());  imgData.push_back(ig);
            int ib = int(255.99 * color.b());  imgData.push_back(ib);
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
        rtwProgress.Refresh(g_Height - j);
    }

    imageFile.close();
    imgFilePath.replace(imgFilePath.find(".ppm"), 4, ".bmp");
    stbi_write_bmp(imgFilePath.c_str(), g_Width, g_Height, 3, imgData.data());
    return 0;
}

//Ch7: 模拟杂乱无章随机的漫反射
int Ch7_DiffuseMaterial(std::string imgFilePath){
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
        if(world -> hit(r,0.001,g_MAX_TmFloat,reco)){//Ch7: 0.001f,表示去除靠近0的浮点值，避免浮点精度带来的毛刺
            vec3 target = reco.p_ + reco.normal_ + random_in_unit_sphere();//Ch7:p_+normal_得到hit-point的球心，再随机选个方向作为反射关系
            return 0.5* getColor(ray(reco.p_,target-reco.p_),world);//Ch7:递归调用，即多次反射，直到hit-miss
        }
        else{
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5 * (unit_direction.y() + 1.0);
            return (1.0 - t) * vec3(1.0,1.0,1.0) + t *vec3(0.5,0.7,1.0);  
        } 
    };
    
    std::vector<unsigned char> imgData;
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
        for(int i = 0 ; i < g_Width; ++i){
            vec3 color(0,0,0);
            for(int s = 0; s< g_RayNums ; ++s){
                float u = float(i + random_double())/ float(g_Width);
                float v = float(j + random_double())/ float(g_Height);
                ray r = cam.get_ray(u,v);
                color += getColor(r,world.get());
            }
            color /= float(g_RayNums);
            int ir = int(255.99 * color.r());   imgData.push_back(ir);
            int ig = int(255.99 * color.g());   imgData.push_back(ig);
            int ib = int(255.99 * color.b());   imgData.push_back(ib);
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    imgFilePath.replace(imgFilePath.find(".ppm"), 4, ".bmp");
    stbi_write_bmp(imgFilePath.c_str(), g_Width, g_Height, 3, imgData.data());
    return 0;
}

// 物体的材质讨论, 这里加入2种: 金属材质的高光,粗糙表面Lambert表面的漫反射
int Ch8_MaterialMetal(std::string imgFilePath){
    RtwProgress rtwProgress(imgFilePath, g_Height);
   //Ch8: modify getColor()
   using getColorFuncType = std::function<vec3(const ray&r,hitable *world,int depth)>;
    getColorFuncType getColor = [&](const ray&r,hitable *world,int depth) -> vec3{
        hit_record reco;
        //Ch6:根据光线击中的最近点，进行渲染着色
        if(world -> hit(r,0.001,g_MAX_TmFloat,reco)){//Ch7: 0.001f,表示去除靠近0的浮点值，避免浮点精度带来的毛刺
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
    
    std::vector<unsigned char> imgData;
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
    int nSphereNum = 5;//球的个数
    using hitableRef = std::shared_ptr<hitable>;
    std::vector<hitableRef> list;
    list.resize(nSphereNum);
    list[0] = std::make_shared<sphere>(vec3(0,0,-1),        0.5, std::make_shared<lambertian>(vec3(0.8,0.3,0.3)));
    list[1] = std::make_shared<sphere>(vec3(0,-100.5,-1), 100.0, std::make_shared<lambertian>(vec3(0.8, 0.8, 0.3)));
    list[2] = std::make_shared<sphere>(vec3(1,0,-1),        0.4, std::make_shared<metal>(vec3(0.8, 0.6, 0.2)));
    list[3] = std::make_shared<sphere>(vec3(-1, 0, -1),     0.5, std::make_shared<metal>(vec3(0.8, 0.8, 0.8)));
    list[4] = std::make_shared<sphere>(vec3(-0.5,-0.4,-0.5),0.1, std::make_shared<lambertian>(vec3(0.2,1.0,1.0)));

    std::shared_ptr<hitable> world = std::make_shared<hitable_list>(list.data(), nSphereNum);
    camera cam;//多条光线打向同一个pixel，模拟MSAA进行抗混叠
    for(int j = g_Height -1 ; j >= 0 ; --j){
        for(int i = 0 ; i < g_Width; ++i){
            vec3 color(0,0,0);
            for(int s = 0; s< g_RayNums ; ++s){
                float u = float(i + random_double())/ float(g_Width);
                float v = float(j + random_double())/ float(g_Height);
                ray r = cam.get_ray(u,v);
                color += getColor(r,world.get(),0);
            }
            color /= float(g_RayNums);
            int ir = int(255.99 * color.r()); imgData.push_back(ir);
            int ig = int(255.99 * color.g()); imgData.push_back(ig);
            int ib = int(255.99 * color.b()); imgData.push_back(ib);
            imageFile << ir << " " << ig << " " << ib << "\n"; 
        }
        rtwProgress.Refresh(g_Height - j);
    }
    imageFile.close();
    imgFilePath.replace(imgFilePath.find(".ppm"), 4, ".bmp");
    stbi_write_bmp(imgFilePath.c_str(), g_Width, g_Height, 3, imgData.data());
    return 0;
}




/* Ch9 透明介质(dielectrics),涉及光的折射refract
 * 复习初中物理知识: 空气折射率1.0 ,玻璃1.3~1.7，钻石2.4
 * 折射定律(也叫斯涅尔定律): 入射光与折射光线位于分界线两侧,分界线垂直于法线
 *                         2种材质的折射率与角度n1*sin(x1) = n2*sin(x2)
                           相对折射率:n21 = sin(x1)/sin(x2) (1介质到2介质)
*/
int Ch9_Dielectrics(std::string imgFilePath) {
    RtwProgress rtwProgress(imgFilePath, g_Height);
    //Ch8: modify getColor()
    using getColorFuncType = std::function<vec3(const ray&r, hitable *world, int depth)>;
    getColorFuncType getColor = [&](const ray&r, hitable *world, int depth) -> vec3 {
        hit_record reco;
        if (depth > g_DepthThreshold) {
            return color(0, 0, 0);
        }
        //Ch6:根据光线击中的最近点，进行渲染着色
        if (world->hit(r, 0.001, g_MAX_TmFloat, reco)) {//Ch7: 0.001f,表示去除靠近0的浮点值，避免浮点精度带来的毛刺
            ray scattered;
            vec3 attenuation;//Ch8 : 材料属性,反射率,吸光率
            if (reco.mate_ptr->scatter(r, reco, attenuation, scattered)) {
                return attenuation * getColor(scattered, world, depth + 1);//递归,继续反射
            }
            return color(0, 0, 0);
        }

        vec3 unit_direction = unit_vector(r.direction());
        float t = 0.5 * (unit_direction.y() + 1.0);
        return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
    };

    std::vector<unsigned char> imgData;
    if (access(imgFilePath.c_str(), 0) == 0) {
        std::remove(imgFilePath.c_str());
    }
    std::ofstream imageFile(imgFilePath);
    imageFile << "P3\n" << g_Width << " " << g_Height << "\n255\n";
    vec3 lower_left_corner_P(-2.0, -1.0, -1.0);
    vec3 horizontalDir(4.0, 0.0, 0.0);
    vec3 verticalDir(0.0, 2.0, 0.0);
    vec3 originP(0.0, 0.0, 0.0);

    // multi-object
    using hitableRef = std::shared_ptr<hitable>;
    std::vector<hitableRef> list;
    list.push_back(std::make_shared<sphere>(vec3(0, 0, -1),       0.5, std::make_shared<lambertian>(vec3(0.1, 0.2, 0.5))));
    list.push_back(std::make_shared<sphere>(vec3(0, -100.5,-1),   100, std::make_shared<lambertian>(vec3(0.8, 0.8, 0.0))));
    list.push_back(std::make_shared<sphere>(vec3(0.5,-0.4,-0.5),  0.1, std::make_shared<lambertian>(vec3(0.2, 1.0, 1.0))));
    list.push_back(std::make_shared<sphere>(vec3(1 , 0, -1),      0.5, std::make_shared<metal>(vec3(0.8, 0.6, 0.2))));
    list.push_back(std::make_shared<sphere>(vec3(-1, 0, -1),      0.5, std::make_shared<dielectric>(1.5)));
    list.push_back(std::make_shared<sphere>(vec3(-1 , 0, -1),   -0.45, std::make_shared<dielectric>(1.5)));

    std::shared_ptr<hitable> world = std::make_shared<hitable_list>(list.data(), list.size());
    camera cam;//多条光线打向同一个pixel，模拟MSAA进行抗混叠
    for (int j = g_Height - 1; j >= 0; --j) {
        for (int i = 0; i < g_Width; ++i) {
            vec3 color(0, 0, 0);
            for (int s = 0; s < g_RayNums; ++s) {
                float u = float(i + random_double()) / float(g_Width);
                float v = float(j + random_double()) / float(g_Height);
                ray r = cam.get_ray(u, v);
                color += getColor(r, world.get(), 0);
            }
            color /= float(g_RayNums);
            int ir = int(255.99 * color.r()); imgData.push_back(ir);
            int ig = int(255.99 * color.g()); imgData.push_back(ig);
            int ib = int(255.99 * color.b()); imgData.push_back(ib);
            imageFile << ir << " " << ig << " " << ib << "\n";
        }
        rtwProgress.Refresh(g_Height - j);
    }

    imageFile.close();
    imgFilePath.replace(imgFilePath.find(".ppm"), 4, ".bmp");
    stbi_write_bmp(imgFilePath.c_str(), g_Width, g_Height, 3, imgData.data());
    return 0;
}