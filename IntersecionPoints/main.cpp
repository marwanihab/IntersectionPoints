//
//  main.cpp
//  Assignment(1)
//
//  Created by Marwan Ihab on 9/23/18.
//  Copyright © 2018 Marwan Ihab. All rights reserved.
//
#include <cmath>
#include <iostream>
#define PI 3.14159265

class vector{
public: double x;
public: double y;
public: double z;
    

};

double dotProduct(vector vectorA, vector vectorB){
    return vectorA.x*vectorB.x+vectorA.y*vectorB.y+vectorA.z*vectorB.z;
}

void gettingIntersectionpointOfaRayCone(double p0x, double p0y,double p0z , double p1x,double p1y, double p1z, double radius , double centreX , double centreY,  double centreZ, double height){
    
    double a,b,c,numerator,cosine,Vmagnitude,Dmagnitude;
    
//    vertexY = centreY + height;
    
    //ray vector
    vector D;
    D.x=p1x-p0x;
    D.y=p1y-p0y;
    D.z=p1z-p0z;
    Dmagnitude = sqrt(pow(D.x, 2) + pow(D.y, 2) + pow(D.z, 2));
    D.x=D.x/Dmagnitude;
    D.y=D.y/Dmagnitude;
    D.z=D.z/Dmagnitude;
    
    
    //cone vector
    vector V;
    V.x=0;
    V.y=-height;
    V.z=0;
    Vmagnitude = sqrt(pow(V.x, 2) + pow(V.y, 2) + pow(V.z, 2));
    V.x=0;
    V.y=V.y/Vmagnitude;
    V.z=0;
    
    //vector from center to p0
    vector CO;
    CO.x = p0x-centreX;
    CO.y = p0y-(centreY+height);
    CO.z = p0z-centreZ;
    
    numerator = sqrt(pow(height, 2)+pow(radius, 2));
    
    cosine = height/numerator;
    
    double ceta;
    //ceta1 = acos(cosine);
    ceta = acos(cosine) * 180.0 / PI;;
    //double test = cos(ceta1);
    
    a= pow((dotProduct(D, V)), 2)-pow(cosine, 2);
    b= 2*((dotProduct(D, V)*dotProduct(CO, V))- (dotProduct(D, CO)*pow(cosine, 2)));
    c= pow((dotProduct(CO, V)),2) - (dotProduct(CO, CO)*pow(cosine, 2));
    
    double decider,tFactor1, tFactor2, realPart, imaginaryPart,interPointx,interPointy,interPointz;
    a=a;
    b=b;
    c=c;
    decider= pow(b,2)-4.0 * a * c ;
    
    
    
    
    
    vector check;
    
    
    
    if (decider > 0) {
        tFactor1 = (-b + sqrt(decider)) / (2*a);
        tFactor2 = (-b - sqrt(decider)) / (2*a);
        std::cout << "Roots are real and different.\n" << std::endl;
        
       
        interPointx= p0x+(D.x*tFactor1);
        interPointy= p0y+(D.y*tFactor1);
        interPointz= p0z+(D.z*tFactor1);
        
        check.x = interPointx-centreX;
        check.y = interPointy-centreY;
        check.z = interPointz-centreZ;
        
        //if (ceta < 90) {
            //double result = dotProduct(check, V);
//            if (result > 0) {
            
                std::cout << "t1 = \n" << tFactor1 << std::endl;
                std::cout << "1st possible intersection point:\n" << std::endl;
                std::cout <<"X ="<< interPointx <<std::endl;
                std::cout <<"Y ="<< interPointy <<std::endl;
                std::cout <<"Z ="<< interPointz <<std::endl;
                
        //}
        //}
       
        
        
        
        
        
       
        interPointx= p0x+(D.x*tFactor2);
        interPointy= p0y+(D.y*tFactor2);
        interPointz= p0z+(D.z*tFactor2);
        
        check.x = interPointx-centreX;
        check.y = interPointy-centreY;
        check.z = interPointz-centreZ;
        
        //if (ceta < 90) {
            //double result = dotProduct(check, V);
            //if (result > 0) {
                
                std::cout << "t2 = \n" << tFactor2 << std::endl;
                std::cout << "2nd possible intersection point:\n" << std::endl;
                std::cout <<"X ="<< interPointx <<std::endl;
                std::cout <<"Y ="<< interPointy <<std::endl;
                std::cout <<"Z ="<< interPointz <<std::endl;
                
            //}
        //}
        
        
        
        
    }
    
    else if (decider == 0) {
        std::cout << "Roots are real and same.\n" << std::endl;
        tFactor1 = (-b + sqrt(decider)) / (2*a);
       
        interPointx= p0x+(D.x*tFactor1);
        interPointy= p0y+(D.y*tFactor1);
        interPointz= p0z+(D.z*tFactor1);
        
        check.x = interPointx-centreX;
        check.y = interPointy-centreY;
        check.z = interPointz-centreZ;
        
        //if (ceta < 90) {
            double result = dotProduct(check, V);
            //if (result > 0) {
                
                std::cout << "t1 = t2 =\n" << tFactor1 << std::endl;
                std::cout << "possible intersection point:\n" << std::endl;
                std::cout <<"X ="<< interPointx <<std::endl;
                std::cout <<"Y ="<< interPointy <<std::endl;
                std::cout <<"Z ="<< interPointz <<std::endl;
                
            //}
       // }
        
    }
    
    else {
        realPart = -b/(2*a);
        imaginaryPart =sqrt(-decider)/(2*a);
        std::cout << "Roots are complex and different."  << std::endl;
        std::cout << "t1 = " << realPart << "+" << imaginaryPart << "i" << std::endl;
        std::cout << "t2 = " << realPart << "-" << imaginaryPart << "i" << std::endl;
        std::cout <<"no intersection"<<std::endl;
    }

    
    
    
    
    
}

void gettingIntersectionPointGeometriv(double p0x, double p0y,double p0z , double p1x,double p1y, double p1z, double radius , double centreX , double centreY,  double centreZ){
    double interPointX,interPointY,interPointZ,directionVectorMagnitude;
    
    vector A;
    A.x = centreX - p0x;
    A.y = centreY - p0y;
    A.z = centreZ - p0z;
    
    vector V;
    V.x=p1x-p0x;
    V.y=p1y-p0y;
    V.z=p1z-p0z;
    
    directionVectorMagnitude= sqrt(pow(V.x, 2) + pow(V.y, 2) + pow(V.z, 2));
    
    double AV = dotProduct(A, V);
    double AA = dotProduct(A, A);
    double decider= (radius*radius)-(AA-pow((AV/directionVectorMagnitude), 2));
    
    if (decider > 0) {
        std::cout <<"2 intersection points \n"<<std::endl;
        std::cout <<"1st point: \n"<<std::endl;
        
        interPointX = p0x + ((AV/directionVectorMagnitude) - sqrt((radius*radius)-(AA-pow((AV/directionVectorMagnitude), 2))))*(V.x/directionVectorMagnitude);
        interPointY = p0y + ((AV/directionVectorMagnitude) - sqrt((radius*radius)-(AA-pow((AV/directionVectorMagnitude), 2))))*(V.y/directionVectorMagnitude);
        interPointZ = p0z + ((AV/directionVectorMagnitude) - sqrt((radius*radius)-(AA-pow((AV/directionVectorMagnitude), 2))))*(V.z/directionVectorMagnitude);
        
        std::cout <<"X ="<< interPointX <<std::endl;
        std::cout <<"Y ="<< interPointY <<std::endl;
        std::cout <<"Z ="<< interPointZ <<std::endl;
        
        std::cout <<"2nd point: \n"<<std::endl;
        
        interPointX = p0x + ((AV/directionVectorMagnitude) + sqrt((radius*radius)-(AA-pow((AV/directionVectorMagnitude), 2))))*(V.x/directionVectorMagnitude);
        interPointY = p0y + ((AV/directionVectorMagnitude) + sqrt((radius*radius)-(AA-pow((AV/directionVectorMagnitude), 2))))*(V.y/directionVectorMagnitude);
        interPointZ = p0z + ((AV/directionVectorMagnitude) + sqrt((radius*radius)-(AA-pow((AV/directionVectorMagnitude), 2))))*(V.z/directionVectorMagnitude);
        
        std::cout <<"X ="<< interPointX <<std::endl;
        std::cout <<"Y ="<< interPointY <<std::endl;
        std::cout <<"Z ="<< interPointZ <<std::endl;
        
    }
    else{
        if (decider ==0) {
            std::cout <<"1 intersection point \n"<<std::endl;
            std::cout <<"The point: \n"<<std::endl;
            
            interPointX = p0x + ((AV/directionVectorMagnitude) - sqrt((radius*radius)-(AA-pow((AV/directionVectorMagnitude), 2))))*(V.x/directionVectorMagnitude);
            interPointY = p0y + ((AV/directionVectorMagnitude) - sqrt((radius*radius)-(AA-pow((AV/directionVectorMagnitude), 2))))*(V.y/directionVectorMagnitude);
            interPointZ = p0z + ((AV/directionVectorMagnitude) - sqrt((radius*radius)-(AA-pow((AV/directionVectorMagnitude), 2))))*(V.z/directionVectorMagnitude);
            
            std::cout <<"X ="<< interPointX <<std::endl;
            std::cout <<"Y ="<< interPointY <<std::endl;
            std::cout <<"Z ="<< interPointZ <<std::endl;
            
            
        }else{
            std::cout <<"no solution \n"<<std::endl;
        }
    }
    
    
    
    
   
    
    
    
    
    
}

void gettingIntersectionPointAlgebric(double p0x, double p0y,double p0z , double p1x,double p1y, double p1z, double radius , double centreX , double centreY,  double centreZ){
    //Vx,Vy,Vz
    double a,b,c,directionVectorMagnitude,vectorFromOriginToP0MagnitudeSquared;
    
    //vector v coordinates
    vector V;
    V.x=p1x-p0x;
    V.y=p1y-p0y;
    V.z=p1z-p0z;
    
    
    //[P0-Origin]
    vector vectorFromOriginToP0;
    vectorFromOriginToP0.x=p0x-centreX;
    vectorFromOriginToP0.y=p0y-centreY;
    vectorFromOriginToP0.z=p0z-centreZ;
    
    ////b coefficient
    b= 2*dotProduct(V, vectorFromOriginToP0);
    
    
    directionVectorMagnitude= sqrt(pow(V.x, 2) + pow(V.y, 2) + pow(V.z, 2));
    
    ////a coefficient
    a=pow(directionVectorMagnitude,2);
    
    vectorFromOriginToP0MagnitudeSquared= pow(vectorFromOriginToP0.x, 2) + pow(vectorFromOriginToP0.y, 2) + pow(vectorFromOriginToP0.z, 2);
    
    ////c coefficient
    c= vectorFromOriginToP0MagnitudeSquared-(radius*radius);
    
//    double tFactor = solvingQuadraticEqaution(a,b,c);
    
    double decider,tFactor1, tFactor2, realPart, imaginaryPart,interPointx,interPointy,interPointz;
    
    decider= (b*b)-4*a*c;
    
    if (decider > 0) {
        tFactor1 = (-b + sqrt(decider)) / (2*a);
        tFactor2 = (-b - sqrt(decider)) / (2*a);
        std::cout << "Roots are real and different.\n" << std::endl;
        std::cout << "t1 = \n" << tFactor1 << std::endl;
        std::cout << "1st possible intersection point:\n" << std::endl;
        interPointx= p0x+(V.x*tFactor1);
        interPointy= p0y+(V.y*tFactor1);
        interPointz= p0z+(V.z*tFactor1);
        std::cout <<"X ="<< interPointx <<std::endl;
        std::cout <<"Y ="<< interPointy <<std::endl;
        std::cout <<"Z ="<< interPointz <<std::endl;
        
        std::cout << "t2 = \n" << tFactor2 << std::endl;
        std::cout << "2nd possible intersection point:\n" << std::endl;
        interPointx= p0x+(V.x*tFactor2);
        interPointy= p0y+(V.y*tFactor2);
        interPointz= p0z+(V.z*tFactor2);
        std::cout <<"X ="<< interPointx <<std::endl;
        std::cout <<"Y ="<< interPointy <<std::endl;
        std::cout <<"Z ="<< interPointz <<std::endl;
        
        
        
    }
    
    else if (decider == 0) {
        std::cout << "Roots are real and same.\n" << std::endl;
        tFactor1 = (-b + sqrt(decider)) / (2*a);
        std::cout << "t1 = t2 =\n" << tFactor1 << std::endl;
        std::cout << "possible intersection point:\n" << std::endl;
        interPointx= p0x+(V.x*tFactor1);
        interPointy= p0y+(V.y*tFactor1);
        interPointz= p0z+(V.z*tFactor1);
        std::cout <<"X ="<< interPointx <<std::endl;
        std::cout <<"Y ="<< interPointy <<std::endl;
        std::cout <<"Z ="<< interPointz <<std::endl;
    }
    
    else {
        realPart = -b/(2*a);
        imaginaryPart =sqrt(-decider)/(2*a);
        std::cout << "Roots are complex and different."  << std::endl;
        std::cout << "t1 = " << realPart << "+" << imaginaryPart << "i" << std::endl;
        std::cout << "t2 = " << realPart << "-" << imaginaryPart << "i" << std::endl;
        std::cout <<"no intersection"<<std::endl;
    }

    
    
    
}






int main(int argc, const char * argv[]) {
    
    std::cout << "Please enter '1' for ray-sphere intersection and '2' for ray-cone intersection !!! \n";
    int interDeterminant;
    std::cin  >> interDeterminant;
    
    switch (interDeterminant) {
        case 1:
          
    std::cout << "Please enter the center of the sphere coordinates( x then y then z ) !!! \n";
    double centreX,centreY,centreZ;
    std::cin  >> centreX;
    std::cin  >> centreY;
    std::cin  >> centreZ;
    
    std::cout << "Please enter the radius !!! \n";
    double radius;
    std::cin  >> radius;
    
    std::cout << "Please enter P0 coordinates( x then y then z )!!! \n";
    double p0x,p0y,p0z;
    std::cin  >> p0x;
    std::cin  >> p0y;
    std::cin  >> p0z;
    
    std::cout << "Please enter P1 coordinates( x then y then z )!!! !!! \n";
    double p1x,p1y,p1z;
    std::cin  >> p1x;
    std::cin  >> p1y;
    std::cin  >> p1z;
    
    std::cout << "Please enter '1' for geometric method and '2' for algebric method !!! \n";
    int methodInput;
    std::cin  >> methodInput;
    
    switch (methodInput) {
        case 2:
            gettingIntersectionPointAlgebric(p0x, p0y, p0z, p1x, p1y, p1z, radius, centreX, centreY, centreZ);
            break;
        case 1:
            gettingIntersectionPointGeometriv(p0x, p0y, p0z, p1x, p1y, p1z, radius, centreX, centreY, centreZ);
            break;
        default:
            break;
    }   break;
    
       case 2:
            
            std::cout << "Please enter the center of the cone coordinates( x then y then z ) !!! \n";
            double centreXX,centreYY,centreZZ;
            std::cin  >> centreXX;
            std::cin  >> centreYY;
            std::cin  >> centreZZ;
            
            std::cout << "Please enter the radius !!! \n";
            double radiusCone;
            std::cin  >> radiusCone;
            
            std::cout << "Please enter P0 coordinates( x then y then z )!!! \n";
            double p0xx,p0yy,p0zz;
            std::cin  >> p0xx;
            std::cin  >> p0yy;
            std::cin  >> p0zz;
            
            std::cout << "Please enter P1 coordinates( x then y then z )!!! !!! \n";
            double p1xx,p1yy,p1zz;
            std::cin  >> p1xx;
            std::cin  >> p1yy;
            std::cin  >> p1zz;
            
            std::cout << "Please enter the height of the cone !!! \n";
            double height;
            std::cin  >> height;
            
            gettingIntersectionpointOfaRayCone( p0xx,  p0yy,  p0zz,  p1xx,  p1yy,  p1zz,  radiusCone,  centreXX,  centreYY,  centreZZ, height);
            
        break;
            
    
   default:
    break;
    }
    
    while (true) {
        
    }
    
    
    return 0;
}
