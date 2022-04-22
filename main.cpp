#include "vector.h"//delete
#include "integrate.h"



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

    Vector radius = { -6960.29650278241, 35.40570314626 , 15.38416914911 };
    Vector speed = { .02180701116 ,  1.00761760916  , 7.50001035085 };
    double time = 2459935.3750000;
    RungeCutte8(radius, speed, time + 2400/86400./*10 min*/, 1, de405);

    
    
    /*Vector r_grinvich = J200toGCS(radius, time);
    Vector harmonic_accel_earth_grinvich = CalcHarmonicAcceleratin(mu, Rz, r_grinvich, 9, Cnn, Snn);
    Vector harmonic_accel_earth = GCStoJ2000(harmonic_accel_earth_grinvich, time);
    harmonic_accel_earth.Print();*/

}
