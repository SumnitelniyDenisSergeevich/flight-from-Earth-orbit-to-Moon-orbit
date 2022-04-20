// Earth_Moon.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <iostream>
#include <fstream>
#include "vector.h"
#include "integrate.h"
#include "matrix.h"

#include "change_coord_system.h"

#include "integrate.h"//delete
#include "constants.h" //delete



using namespace std;

int main()
{
    dph::EphemerisRelease de405("JPLEPH");
    if (!de405.isReady())
    {
        throw runtime_error("open file error");
    }
   
    //J200toGCS({ -6960.29650278241, 35.40570314626, 15.38416914911 }, { .02180701116, 1.00761760916, 7.50001035085 }, 2459935.375);  
    //J200toGCS({ 1000, 9000, 3000 }, { 3,1,7 }, 2451544.5);
    //J200toGCS({ 8000, 0, 0 }, { 0,4,6 }, 2451544.5); 

    Vector radius = { -6960.29650278241, 35.40570314626, 15.38416914911 };
    Vector speed = { .02180701116 ,  1.00761760916  , 7.50001035085 };
    double time = 2459935.375000000; // JED
    RungeCutte8(radius, speed, time + 14/*10 min*/, 10, de405);

    //Vector radius = { -5066.2156143621, -3568.3146432419, -3148.7370464944 };
    //double time = 2458823.375; // JED
    //Vector r_grinvich = J200toGCS(radius, time);

    //Vector harmonic_accel_earth_grinvich = CalcHarmonicAcceleratin(mu, Rz, r_grinvich, 33, Cnn, Snn);


}
