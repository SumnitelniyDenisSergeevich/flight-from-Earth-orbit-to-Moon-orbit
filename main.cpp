#include "vector.h"//delete
#include "change_coord_system.h" // delete
#include "integrate.h"
#include <iostream>

#include "log_duration.h"


using namespace std;

int main()
{
    dph::EphemerisRelease de405("JPLEPH");
    if (!de405.isReady())
    {
        throw runtime_error("open file error");
    }
    cout.precision(15);
   
    //J200toGCS({ -6960.29650278241, 35.40570314626, 15.38416914911 }, { .02180701116, 1.00761760916, 7.50001035085 }, 2459935.375);  
    //J200toGCS({ 1000, 9000, 3000 }, { 3,1,7 }, 2451544.5);
    //J200toGCS({ 8000, 0, 0 }, { 0,4,6 }, 2451544.5); 
    /*cout << ConvertDateToUliane(23,8,2025,19,42,32,3) - 2460905.63401414 << endl;
    system("pause");*/
    // Земная орбита
    /*Vector radius = { -6960.29650278241, 35.40570314626 , 15.38416914911 };
    Vector speed = { .02180701116 ,  1.00761760916  , 7.50001035085 };
    double time = 2459935.3750000;*/
    //На Луну 2460905.63401414
    Vector radius = { 6274.036317, 3774.014834, 6577.356959 };
    Vector speed = { -2.211116, 5.5672, 6.566456 };
    double time = 2460905.63401414;

    {
        LOG_DURATION("runge");
        RungeCutte8(radius, speed, time, time + 6, 1, de405);
    }

    //Vector radius = { 2,2,5 };
    //Vector planet = { 6,4,6 };
    //double planet_radius = 1.5;
    //CalcVisiblePlanetAngle(radius, planet, planet_radius);
    
    
    /*Vector r_grinvich = J200toGCS(radius, time);
    Vector harmonic_accel_earth_grinvich = CalcHarmonicAcceleratin(mu, Rz, r_grinvich, 9, Cnn, Snn);
    Vector harmonic_accel_earth = GCStoJ2000(harmonic_accel_earth_grinvich, time);
    harmonic_accel_earth.Print();*/

}
