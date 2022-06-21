//---------------------------------------------------------------------------
#pragma hdrstop
#include <math.h>
#include "coordinate_system.h"
#include "mathematics.h"
#include "constants.h"
#include "time.h"
#include "matrix.h"
#include "jpl_eph.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)

double DegToRad(const double deg) {
	return deg / 180 * pi;
}

void moonSCS(double flr[3], double r[3]){

double f=1./800.;
double e=sqrt(f*(2-f));
double x=1/sqrt(1-e*e*sin(flr[0])*sin(flr[0]));
double R=1738;
r[0] = (R*x+flr[2])*cos(flr[0])*cos(flr[1]);
r[1] = (R*x+flr[2])*cos(flr[0])*sin(flr[1]);
r[2] = (R*x*(1+e*e)+flr[2])*sin(flr[0]);



}
//----------------------------------------------------------------------------
//Преобразование сферических координат в декартовы
void philamTOxyz (double flr[3], double r[3])
				 //flr[0] - широта в радианах
				 //flr[1] - долгота в радианах
				 //flr[2] - радиус в километрах
				 //r[0] - проекция вектора на ось х
				 //r[1] - проекция вектора на ось y
				 //r[3] - проекция вектора на ось z
{
	r[0] = flr[2]*cos(flr[0])*cos(flr[1]);
	r[1] = flr[2]*cos(flr[0])*sin(flr[1]);
	r[2] = flr[2]*sin(flr[0]);
}
//----------------------------------------------------------------------------
//Преобразование сферических координат в декартовы
void philamTOxyz (	const double &phi,	//широта в радианах
					const double &lam,  //долгота в радианах
					const double &r,    //радиус в километрах
					double &x,          //проекция вектора на ось х
					double &y,          //проекция вектора на ось y
					double &z)          //проекция вектора на ось z
{
	x = r*cos(phi)*cos(lam);
	y = r*cos(phi)*sin(lam);
	z = r*sin(phi);
}
//-------------------------------------------------------------------------
// Преобразование декартовых координат в сферических
void xyzTOphilam (double r[3], double flr[3])
				 //flr[0] - широта в радианах
				 //flr[1] - долгота в радианах
				 //flr[2] - радиус в километрах
				 //r[0] - проекция вектора на ось х
				 //r[1] - проекция вектора на ось y
				 //r[3] - проекция вектора на ось z
{
	flr[0] = atan(r[2]/sqrt(r[0]*r[0] + r[1]*r[1]));
	flr[1] = atan2(r[1], r[0]);
	flr[2] = norm(r);
}
//-------------------------------------------------------------------------
// Преобразование декартовых координат в сферических
void xyzTOphilam (double r[3], double flr[3], double DS_DX[9])
				 //flr[0] - широта в радианах
				 //flr[1] - долгота в радианах
				 //flr[2] - радиус в километрах
				 //r[0] - проекция вектора на ось х
				 //r[1] - проекция вектора на ось y
				 //r[3] - проекция вектора на ось z
				 //DS_DX - матрица частных производных сферических координат flr по компонентам вектора r
{
	flr[0] = atan(r[2]/sqrt(r[0]*r[0] + r[1]*r[1]));
	flr[1] = atan2(r[1], r[0]);
	flr[2] = norm(r);


double R=norm(r);
double xy=sqrt(r[0]*r[0]+r[1]*r[1]);

   //радиус
   DS_DX[0]=r[0]/R;
   DS_DX[1]=r[1]/R;
   DS_DX[2]=r[2]/R;
   //долгота
   DS_DX[3]=-r[1]/(xy*xy);
   DS_DX[4]=r[0]/(xy*xy);
   DS_DX[5]=0.0;
   //широта
   DS_DX[6]=-(r[2]*r[0])/(R*R*xy);
   DS_DX[7]=-(r[2]*r[1])/(R*R*xy);
   DS_DX[8]=xy/(R*R);








}
// Преобразование декартовых координат в сферических
void xyzTOphilam (	const double &x,   //проекция вектора на ось х
					const double &y,   //проекция вектора на ось y
					const double &z,   //проекция вектора на ось z
					double &phi,       //широта в радианах
					double &lam,       //долгота в радианах
					double &r)         //радиус в километрах
{
	phi = atan(z/sqrt(x*x + y*y));
	lam = atan2(y, x);
	r = norm(x,y,z);
}
//---------------------------------------------------------------------------
/*Перевод декартовых координат в эллипсоидальные
Входные параметры:
		  Dec - массив декартовых координат (x, y, z) (километры)
Выходные параметры:
			Ell - массив эллипсоидальных координат (fi, lambda, h)
				  fi     - широта  (радианы)
				  lambda - долгота (радианы)
				  h      - высота  (километры)
Подключаемые модули:
		<math.h>
*/
void DecToEll(const double Dec[3], double Ell[3])
{
	double x  = Dec[0]*1000.;
	double y  = Dec[1]*1000.;
	double z  = Dec[2]*1000.;

	double a = 6378137.;
	double f = 1./298.257223563;
	double e = sqrt(2*f-f*f);

	double lambda;
	if( x == 0.)
		{
			if( y > 0.) lambda = M_PI/2.;
			if( y < 0.) lambda = M_PI*3./2.;
		}
	else
		{
			lambda = atan2( y, x );
		}

	double fi = atan(z/sqrt(x*x+y*y)); //стартовое значение для цикла
	double fi1, N, h;
	fi1=fi+1;
   //	for( ;  fabs(fi-fi1) > 10e-10  ; )
   while(fabs(fi-fi1)>10e-8)
		{
			fi1=fi;
			N = a/sqrt(1. - e*e*sin(fi1)*sin(fi1));
			h = ( sqrt(x*x+y*y) / cos(fi1) ) - N;
			fi = atan(z/sqrt(x*x+y*y) / (1.-e*e* N/(N+h)));
		}

	Ell[0] = fi;
	Ell[1] = lambda;
	Ell[2] = h/1000.;
}
//---------------------------------------------------------------------------
/*Перевод эллипсоидальных координат в декартовы
Входные параметры:
			Ell - массив эллипсоидальных координат (fi, lambda, h)
				  fi     - широта  (радианы)
				  lambda - долгота (радианы)
				  h      - высота  (километры)
Выходные параметры:
		  Dec - массив декартовых координат (x, y, z) (километры)
Подключаемые модули:
		<math.h>
*/
void EllToDec(const double Ell[3], double Dec[3])
{
	double fi      = Ell[0];
	double lambda  = Ell[1];
	double h       = Ell[2]*1000.;

	double a = 6378137.;
	double f = 1./298.257223563;
	double e = sqrt(2*f-f*f);

	double N       = a/sqrt(1.-e*e*sin(fi)*sin(fi));

	double x       = (N + h) * cos(fi) * cos(lambda);
	double y       = (N + h) * cos(fi) * sin(lambda);
	double z       = ( N*( 1.-e*e ) + h ) * sin(fi) ;

	Dec[0] = x/1000.;
	Dec[1] = y/1000.;
	Dec[2] = z/1000.;
}
//---------------------------------------------------------------------------
/*Преобразование Гринвической во Вторую Экваториальную
(перевод только вектора)
Входные параметры:
		Ar1[3] - трехмерный вектор (x, y, z) в километрах
		t - время в юлиянских сутках
Выходные параметры:
		Ar2[3] - трехмерный вектор (x, y, z) в километрах
Используемые процедуры:
		precessionVSnutation - поик звездного времени
Подключаемые модули:
		<math.h>
*/
void GEOtoGEI(const double Ar1[],const double t, double Ar2[]){   //  там ли минус????
	 double cs, ss, s,
			S[3][3], NP[3][3],
			PNS[3][3], SNP[3][3];
precessionVSnutation(t, s, NP);
cs = cos(s);
ss = sin(s);
S[0][0]=cs;		S[0][1]=-ss; 	S[0][2]=0;
S[1][0]=ss;		S[1][1]=cs; 	S[1][2]=0;//  там ли минус????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
S[2][0]=0;		S[2][1]=0; 	 	S[2][2]=1;
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++)   {
		PNS[i][c]=0;
		for(int j=0; j<3; j++)
			PNS[i][c]+=NP[j][i]*S[j][c];
	}
for(int i=0; i<3; i++){
	Ar2[i]=0;
	for(int j=0; j<3; j++)
		Ar2[i]+=PNS[i][j]*Ar1[j];
}
}
//-------------------------------------------------------------------
/*Преобразование Второй Экваториальной в Гринвическую
(перевод только вектора)
Входные параметры:
		Ar1[3] - трехмерный вектор (x, y, z) в километрах
		t - время в юлиянских сутках
Выходные параметры:
		Ar2[3] - трехмерный вектор (x, y, z) в километрах
Используемые процедуры:
		precessionVSnutation - поик звездного времени
Подключаемые модули:
		<math.h>
*/
void GEItoGEO(const double Ar1[],const double t, double Ar2[]){
	double cs, ss, s,
		   S[3][3], NP[3][3],
		   SNP[3][3];
precessionVSnutation(t, s, NP);
cs = cos(s);
ss = sin(s);
S[0][0]=cs;		S[0][1]=ss; 	S[0][2]=0;
S[1][0]=-ss;	S[1][1]=cs; 	S[1][2]=0;
S[2][0]=0;		S[2][1]=0; 	 	S[2][2]=1;
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++)  {
		SNP[i][c]=0;
		for(int j=0; j<3; j++)
			SNP[i][c]+=S[i][j]*NP[j][c];
	}
for(int i=0; i<3; i++){
	Ar2[i]=0;
	for(int j=0; j<3; j++)
		Ar2[i]+=SNP[i][j]*Ar1[j];
}
}
//----------------------------------------------------------------------------
/*Процедура вычисления истинного звездного времени и матрицы прицесси и нутации
Входные параметры:
		t - время в юлиянских сутках
Выходные параметры:
		s - звездное время в радианах
		NP - матрица 3х3 результат произведения матрицы нутации на матрицу прицессии
Используемые процедуры:
		UTCtoTDB - перевод UTC в динамическое время
		GMST - процедура вычисления среднего звездного времени
Подключаемые модули:
		<math.h>
		#include "time.h"
*/
void precessionVSnutation(const double t,double &s, double NP[3][3]){
double ksi, zeta, teta, psi, eps, eps0;
double 	c_ksi, s_ksi,
		c_zeta, s_zeta,
		c_teta, s_teta,
		c_psi, s_psi,
		c_eps, s_eps,
		c_eps0, s_eps0;
double T, T2, T3;
double Ml,Ms,ul,Ds,Uzl;
double P[3][3], N[3][3];
T   = (UTCtoTDB(t) - 2451545.0)/36525.0;
T2=T*T;
T3=T2*T;
ksi  = 0.011180860*T + 1.464E-6*T2 + 8.70E-8*T3;
zeta = 0.011180860*T + 5.308E-6*T2 + 8.90E-8*T3;
teta = 0.009717173*T - 2.068E-6*T2 - 2.02E-7*T3;
eps0 = 0.4090928042 - 0.2269655E-3*T - 2.86E-9*T2 + 8.80E-9*T3;

Ml   = 2.355548394 + (1325*pi_2 + 3.470890873)*T + 1.517952E-4*T2 + 3.103E-7*T3;
Ms   = 6.240035940 +   (99*pi_2 + 6.266610600)*T - 2.797400E-6*T2 - 5.820E-8*T3;
ul   = 1.627901934 + (1342*pi_2 + 1.431476084)*T - 6.427170E-5*T2 + 5.340E-8*T3;
Ds   = 5.198469514 + (1236*pi_2 + 5.360106500)*T - 3.340860E-5*T2 + 9.220E-8*T3;
Uzl  = 2.182438624 -    (5*pi_2 + 2.341119397)*T + 3.614290E-5*T2 + 3.880E-8*T3;

eps = 0; psi = 0;
for(int i=0;i<=105;i++){
	psi +=   (coef_A[i] + coef_B[i]*T)*sin(coef_a1[i]*Ml + coef_a2[i]*Ms +coef_a3[i]*ul +coef_a4[i]*Ds +coef_a5[i]*Uzl);
	eps +=   (coef_C[i] + coef_D[i]*T)*cos(coef_a1[i]*Ml + coef_a2[i]*Ms +coef_a3[i]*ul +coef_a4[i]*Ds +coef_a5[i]*Uzl);
}
psi  = psi*1E-4/3600./180.*pi;
eps  = eps*1E-4/3600./180.*pi;
eps += eps0;

c_ksi=cos(ksi);
s_ksi=sin(ksi);
c_zeta=cos(zeta);
s_zeta=sin(zeta);
c_teta=cos(teta);
s_teta=sin(teta);
c_psi=cos(psi);
s_psi=sin(psi);
c_eps=cos(eps);
s_eps=sin(eps);
c_eps0=cos(eps0);
s_eps0=sin(eps0);

s=GMST(t)+psi*c_eps;

P[0][0] =  c_ksi*c_zeta*c_teta - s_ksi*s_zeta;
P[0][1] = -s_ksi*c_zeta*c_teta - c_ksi*s_zeta;
P[0][2] = -c_zeta*s_teta;

P[1][0] =  c_ksi*s_zeta*c_teta + s_ksi*c_zeta;
P[1][1] = -s_ksi*s_zeta*c_teta + c_ksi*c_zeta;
P[1][2] = -s_zeta*s_teta;

P[2][0] =  c_ksi*s_teta;
P[2][1] = -s_ksi*s_teta;
P[2][2] =  c_teta;

N[0][0] =  c_psi;
N[0][1] = -s_psi*c_eps0;
N[0][2] = -s_psi*s_eps0;

N[1][0] =  s_psi*c_eps;
N[1][1] =  c_psi*c_eps*c_eps0+s_eps*s_eps0;
N[1][2] =  c_psi*c_eps*s_eps0-s_eps*c_eps0;

N[2][0] =  s_psi*s_eps;
N[2][1] =  c_psi*s_eps*c_eps0-c_eps*s_eps0;
N[2][2] =  c_psi*s_eps*s_eps0+c_eps*c_eps0;

for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		NP[i][c]=0;
		for(int j=0; j<3; j++)
			NP[i][c]+=N[i][j]*P[j][c];
	}

}
//---------------------------------------------------------------------------
/*Преобразование Второй Экваториальной в истинную Гринвическую
Входные параметры:
		r[3] - радиус-вектор (x, y, z) [км]
		v[3] - вектор скорости (vx, vy, vz) [км/с]
		t - время в юлиянских сутках
Выходные параметры:
		r1[3] - радиус-вектор (x, y, z) [км]
		v1[3] - вектор скорости (vx, vy, vz) [км/с]
Используемые процедуры:
		UTCtoTDB - перевод UTC в динамическое время
		precessionVSnutation - поик звездного времени
Подключаемые модули:
		<math.h>
		"constants.h"
		#include "time.h"
*/
void GEItoGEO(double r[3], double v[3], double t, double r1[3], double v1[3]){
double ksi, zeta, teta, psi, eps, eps0;
double 	c_ksi, s_ksi,
		c_zeta, s_zeta,
		c_teta, s_teta,
		c_psi, s_psi,
		c_eps, s_eps,
		c_eps0, s_eps0,
		c_l, s_l;
double T, T2, T3;
double xp, yp;
double Ml,Ms,ul,Ds,Uzl,s;
double w=7.292115e-5;   // что это?
double A[3][3], B[3][3], B1[3][3], P[3][3], N[3][3], a[3],
		AB[3][3]={0}, AB1[3][3]={0},
		ABN[3][3]={0}, AB1N[3][3]={0},
		ABNP[3][3]={0}, AB1NP[3][3]={0};
T   = (UTCtoTDB(t) - 2451545.0)/36525.0;
T2=T*T;
T3=T2*T;
ksi  = 0.011180860*T + 1.464E-6*T2 + 8.70E-8*T3;
zeta = 0.011180860*T + 5.308E-6*T2 + 8.90E-8*T3;
teta = 0.009717173*T - 2.068E-6*T2 - 2.02E-7*T3;
eps0 = 0.4090928042 - 0.2269655E-3*T - 2.86E-9*T2 + 8.80E-9*T3;

Ml   = 2.355548394 + (1325*pi_2 + 3.470890873)*T + 1.517952E-4*T2 + 3.103E-7*T3;
Ms   = 6.240035940 +   (99*pi_2 + 6.266610600)*T - 2.797400E-6*T2 - 5.820E-8*T3;
ul   = 1.627901934 + (1342*pi_2 + 1.431476084)*T - 6.427170E-5*T2 + 5.340E-8*T3;
Ds   = 5.198469514 + (1236*pi_2 + 5.360106500)*T - 3.340860E-5*T2 + 9.220E-8*T3;
Uzl  = 2.182438624 -    (5*pi_2 + 2.341119397)*T + 3.614290E-5*T2 + 3.880E-8*T3;

eps = 0; psi = 0;
for(int i=0;i<=105;i++){
	psi +=   (coef_A[i] + coef_B[i]*T)*sin(coef_a1[i]*Ml + coef_a2[i]*Ms +coef_a3[i]*ul +coef_a4[i]*Ds +coef_a5[i]*Uzl);
	eps +=   (coef_C[i] + coef_D[i]*T)*cos(coef_a1[i]*Ml + coef_a2[i]*Ms +coef_a3[i]*ul +coef_a4[i]*Ds +coef_a5[i]*Uzl);
}
psi  = psi*1E-4/3600./180.*pi;
eps  = eps*1E-4/3600./180.*pi;
eps += eps0;

c_ksi=cos(ksi);
s_ksi=sin(ksi);
c_zeta=cos(zeta);
s_zeta=sin(zeta);
c_teta=cos(teta);
s_teta=sin(teta);
c_psi=cos(psi);
s_psi=sin(psi);
c_eps=cos(eps);
s_eps=sin(eps);
c_eps0=cos(eps0);
s_eps0=sin(eps0);



P[0][0] =  c_ksi*c_zeta*c_teta - s_ksi*s_zeta;
P[0][1] = -s_ksi*c_zeta*c_teta - c_ksi*s_zeta;
P[0][2] = -c_zeta*s_teta;

P[1][0] =  c_ksi*s_zeta*c_teta + s_ksi*c_zeta;
P[1][1] = -s_ksi*s_zeta*c_teta + c_ksi*c_zeta;
P[1][2] = -s_zeta*s_teta;

P[2][0] =  c_ksi*s_teta;
P[2][1] = -s_ksi*s_teta;
P[2][2] =  c_teta;

N[0][0] =  c_psi;
N[0][1] = -s_psi*c_eps0;
N[0][2] = -s_psi*s_eps0;

N[1][0] =  s_psi*c_eps;
N[1][1] =  c_psi*c_eps*c_eps0+s_eps*s_eps0;
N[1][2] =  c_psi*c_eps*s_eps0-s_eps*c_eps0;

N[2][0] =  s_psi*s_eps;
N[2][1] =  c_psi*s_eps*c_eps0-c_eps*s_eps0;
N[2][2] =  c_psi*s_eps*s_eps0+c_eps*c_eps0;

s=GMST(t)+psi*c_eps;
c_l=cos(s);
s_l=sin(s);
B[0][0]=c_l;  B[0][1]=s_l; B[0][2]=0;
B[1][0]=-s_l; B[1][1]=c_l; B[1][2]=0;
B[2][0]=0;    B[2][1]=0;   B[2][2]=1;


B1[0][0]=-w*s_l; B1[0][1]=w*c_l;  B1[0][2]=0;
B1[1][0]=-w*c_l; B1[1][1]=-w*s_l; B1[1][2]=0;
B1[2][0]=0;      B1[2][1]=0;      B1[2][2]=0;

xp=0;//0.10157/3600*M_PI/180;
yp=0;//0.49738/3600*M_PI/180;

A[0][0]=1;   A[0][1]=0;  A[0][2]=xp;
A[1][0]=0;   A[1][1]=1;  A[1][2]=-yp;
A[2][0]=-xp; A[2][1]=yp; A[2][2]=1;



for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB[i][c]=0;
		for(int j=0; j<3; j++)
			AB[i][c]+=A[i][j]*B[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		ABN[i][c]=0;
		for(int j=0; j<3; j++)
			ABN[i][c]+=AB[i][j]*N[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		ABNP[i][c]=0;
		for(int j=0; j<3; j++)
			ABNP[i][c]+=ABN[i][j]*P[j][c];
	}
/////////////////////////////////
	for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB1[i][c]=0;
		for(int j=0; j<3; j++)
			AB1[i][c]+=A[i][j]*B1[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB1N[i][c]=0;
		for(int j=0; j<3; j++)
			AB1N[i][c]+=AB1[i][j]*N[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB1NP[i][c]=0;
		for(int j=0; j<3; j++)
			AB1NP[i][c]+=AB1N[i][j]*P[j][c];
	}
//////////////////////////////////////////////
for(int i=0; i<3; i++){
	r1[i]=0;
	for(int j=0; j<3; j++)
		r1[i]+=ABNP[i][j]*r[j];
}

for(int i=0; i<3; i++){
	v1[i]=0;
	for(int j=0; j<3; j++)
		v1[i]+=ABNP[i][j]*v[j]+AB1NP[i][j]*r[j];
}
}
//--------------------------------------------------------------------------
/*Преобразование истинной Гринвической во Вторую Экваториальную
Входные параметры:
		r[3] - радиус-вектор (x, y, z) [км]
		v[3] - вектор скорости (vx, vy, vz) [км/с]
		t - время в юлиянских сутках
Выходные параметры:
		r1[3] - радиус-вектор (x, y, z) [км]
		v1[3] - вектор скорости (vx, vy, vz) [км/с]
Используемые процедуры:
		UTCtoTDB - перевод UTC в динамическое время
		precessionVSnutation - поик звездного времени
Подключаемые модули:
		<math.h>
		"constants.h"
		#include "time.h"
*/
void GEOtoGEI(double r[3], double v[3], double t, double r1[3], double v1[3]){
double ksi, zeta, teta, psi, eps, eps0;
double 	c_ksi, s_ksi,
		c_zeta, s_zeta,
		c_teta, s_teta,
		c_psi, s_psi,
		c_eps, s_eps,
		c_eps0, s_eps0,
		c_l, s_l;
double T, T2, T3;
double xp, yp;
double Ml,Ms,ul,Ds,Uzl,s;
double w=7.292115e-5;
double A[3][3], B[3][3], B1[3][3], P[3][3], N[3][3], a[3],
		AB[3][3]={0}, AB1[3][3]={0},
		ABN[3][3]={0}, AB1N[3][3]={0},
		ABNP[3][3]={0}, AB1NP[3][3]={0};
T   = (UTCtoTDB(t) - 2451545.0)/36525.0;  //Ttdb
T2=T*T;
T3=T2*T;
ksi  = 0.011180860*T + 1.464E-6*T2 + 8.70E-8*T3;
zeta = 0.011180860*T + 5.308E-6*T2 + 8.90E-8*T3;
teta = 0.009717173*T - 2.068E-6*T2 - 2.02E-7*T3;
eps0 = 0.4090928042 - 0.2269655E-3*T - 2.86E-9*T2 + 8.80E-9*T3;

Ml   = 2.355548394 + (1325*pi_2 + 3.470890873)*T + 1.517952E-4*T2 + 3.103E-7*T3;
Ms   = 6.240035940 +   (99*pi_2 + 6.266610600)*T - 2.797400E-6*T2 - 5.820E-8*T3;
ul   = 1.627901934 + (1342*pi_2 + 1.431476084)*T - 6.427170E-5*T2 + 5.340E-8*T3;
Ds   = 5.198469514 + (1236*pi_2 + 5.360106500)*T - 3.340860E-5*T2 + 9.220E-8*T3;
Uzl  = 2.182438624 -    (5*pi_2 + 2.341119397)*T + 3.614290E-5*T2 + 3.880E-8*T3;

eps = 0; psi = 0;
for(int i=0;i<=105;i++){
	psi +=   (coef_A[i] + coef_B[i]*T)*sin(coef_a1[i]*Ml + coef_a2[i]*Ms +coef_a3[i]*ul +coef_a4[i]*Ds +coef_a5[i]*Uzl);
	eps +=   (coef_C[i] + coef_D[i]*T)*cos(coef_a1[i]*Ml + coef_a2[i]*Ms +coef_a3[i]*ul +coef_a4[i]*Ds +coef_a5[i]*Uzl);
}
psi  = psi*1E-4/3600./180.*pi;
eps  = eps*1E-4/3600./180.*pi;
eps += eps0;

c_ksi=cos(ksi);
s_ksi=sin(ksi);
c_zeta=cos(zeta);
s_zeta=sin(zeta);
c_teta=cos(teta);
s_teta=sin(teta);
c_psi=cos(psi);
s_psi=sin(psi);
c_eps=cos(eps);
s_eps=sin(eps);
c_eps0=cos(eps0);
s_eps0=sin(eps0);



P[0][0] =  c_ksi*c_zeta*c_teta - s_ksi*s_zeta;
P[0][1] = -s_ksi*c_zeta*c_teta - c_ksi*s_zeta;
P[0][2] = -c_zeta*s_teta;

P[1][0] =  c_ksi*s_zeta*c_teta + s_ksi*c_zeta;
P[1][1] = -s_ksi*s_zeta*c_teta + c_ksi*c_zeta;
P[1][2] = -s_zeta*s_teta;

P[2][0] =  c_ksi*s_teta;
P[2][1] = -s_ksi*s_teta;
P[2][2] =  c_teta;

N[0][0] =  c_psi;
N[0][1] = -s_psi*c_eps0;
N[0][2] = -s_psi*s_eps0;

N[1][0] =  s_psi*c_eps;
N[1][1] =  c_psi*c_eps*c_eps0+s_eps*s_eps0;
N[1][2] =  c_psi*c_eps*s_eps0-s_eps*c_eps0;

N[2][0] =  s_psi*s_eps;
N[2][1] =  c_psi*s_eps*c_eps0-c_eps*s_eps0;
N[2][2] =  c_psi*s_eps*s_eps0+c_eps*c_eps0;

s=GMST(t)+psi*c_eps;//GMST считает Smean
c_l=cos(s);
s_l=sin(s);
B[0][0]=c_l;  B[0][1]=s_l; B[0][2]=0;
B[1][0]=-s_l; B[1][1]=c_l; B[1][2]=0;
B[2][0]=0;    B[2][1]=0;   B[2][2]=1;


B1[0][0]=-w*s_l; B1[0][1]=w*c_l;  B1[0][2]=0;
B1[1][0]=-w*c_l; B1[1][1]=-w*s_l; B1[1][2]=0;
B1[2][0]=0;      B1[2][1]=0;      B1[2][2]=0;

xp=0;
yp=0;

A[0][0]=1;   A[0][1]=0;  A[0][2]=xp;
A[1][0]=0;   A[1][1]=1;  A[1][2]=-yp;
A[2][0]=-xp; A[2][1]=yp; A[2][2]=1;

for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB[i][c]=0;
		for(int j=0; j<3; j++)
			AB[i][c]+=A[i][j]*B[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		ABN[i][c]=0;
		for(int j=0; j<3; j++)
			ABN[i][c]+=AB[i][j]*N[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		ABNP[i][c]=0;
		for(int j=0; j<3; j++)
			ABNP[i][c]+=ABN[i][j]*P[j][c];
	}
/////////////////////////////////
	for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB1[i][c]=0;
		for(int j=0; j<3; j++)
			AB1[i][c]+=A[i][j]*B1[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB1N[i][c]=0;
		for(int j=0; j<3; j++)
			AB1N[i][c]+=AB1[i][j]*N[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB1NP[i][c]=0;
		for(int j=0; j<3; j++)
			AB1NP[i][c]+=AB1N[i][j]*P[j][c];
	}
//////////////////////////////////////////////
for(int i=0; i<3; i++){
	r1[i]=0;
	for(int j=0; j<3; j++)
		r1[i]+=ABNP[j][i]*r[j];
}
for(int i=0; i<3; i++){
	v1[i]=0;
	for(int j=0; j<3; j++)
		v1[i]+=ABNP[j][i]*v[j]+AB1NP[j][i]*r[j];
}

}


//---------------------------------------------------------------------
//процедура перевода из J2000 в селенографическую СК
void J2000toSG(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]){
double r[3], v[3];
J2000toSC(r_in,v_in,JD, r, v);
SCtoSG(r, v, JD, r_out, v_out);
}
//---------------------------------------------------------------------
//процедура перевода из селенографической СК в J2000
void SGtoJ2000(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]){
double r[3], v[3];

SGtoSC(r_in,v_in,JD, r, v);
SCtoJ2000(r, v, JD, r_out, v_out);
}
//---------------------------------------------------------------------
//процедура перевода из J2000 в селеноцентрическую СК
void J2000toSC(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]){
int	NT=10, 	// moon
	NC=3; 	 //earth
double r_moon[3], v_moon[3];
	PLEPH_RV(JD, NT, NC, r_moon,v_moon);
for(int i=0; i<3; i++){
	r_out[i]=r_in[i]-r_moon[i];
	v_out[i]=v_in[i]-v_moon[i];
}

}

//---------------------------------------------------------------------
//процедура перевода из селеноцентрической СК в J2000
void SCtoJ2000(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]){
int	NT=10, 	// moon
	NC=3; 	//earth
double r_moon[3], v_moon[3];
	PLEPH_RV(JD, NT, NC, r_moon,v_moon);
for(int i=0; i<3; i++){
	r_out[i]=r_in[i]+r_moon[i];
	v_out[i]=v_in[i]+v_moon[i];
}

}
//---------------------------------------------------------------------
//процедура перевода из селеноцентрической СК в селенографическую СК
void SCtoSG(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]){
double R_in[6];
double R_out[6];
for(int i=0; i<3; i++){
	R_in[i]=r_in[i];
	R_in[i+3]=v_in[i];
}
SCtoSG(R_in,JD, R_out);
for(int i=0; i<3; i++){
	 r_out[i]=R_out[i];
	 v_out[i]=R_out[i+3];
}
}
//---------------------------------------------------------------------
//процедура перевода из селеноцентрической СК в селенографическую СК
 void SCtoSG(double R_in[6],double JD, double R_out[6]){

	double dt = JD + (dUTC+dUT1)/86400. - 2451545.;
	double T = dt/36525.;
	double E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13;
	double sin_E1,sin_E2,sin_E3,sin_E4,sin_E6,sin_E10,sin_E13;

	E1= DegToRad(125.045- 0.0529921*dt);
	E2= DegToRad(250.089- 0.1059842*dt);
	E3= DegToRad(260.008+13.0120009*dt);
	E4= DegToRad(176.625+13.3407154*dt);
	E5= DegToRad(357.529+ 0.9856003*dt);
	E6= DegToRad(311.589+26.4057084*dt);
	E7= DegToRad(134.963+13.0649930*dt);
	E8= DegToRad(276.617+ 0.3287146*dt);
	E9= DegToRad( 34.226+ 1.7484877*dt);
	E10=DegToRad( 15.134- 0.1589763*dt);
	E11=DegToRad(119.743+ 0.0036096*dt);
	E12=DegToRad(239.961+ 0.1643573*dt);
	E13=DegToRad( 25.053+12.9590088*dt);

	sin_E1  = sin(E1);
	sin_E2  = sin(E2);
	sin_E3  = sin(E3);
    sin_E4  = sin(E4);
	sin_E6  = sin(E6);
    sin_E10 = sin(E10);
    sin_E13 = sin(E13);

	double alpha = 269.9949 + 0.0031*T      - 3.8787*sin_E1	    - 0.1204*sin_E2
				            + 0.07*sin_E3	- 0.0172*sin_E4	    + 0.0072*sin_E6
											- 0.0052*sin_E10	+ 0.0043*sin_E13;

    double beta  = 66.5392  + 0.013*T	    + 1.5419*cos(E1)	+ 0.0239*cos(E2)
							- 0.0278*cos(E3)+ 0.0068*cos(E4)    - 0.0029*cos(E6)
				            + 0.0009*cos(E7)+ 0.0008*cos(E10)   - 0.0009*cos(E13);

    double om = 38.3213	+ 13.17635815*dt	- 1.4e-12*dt*dt		+ 3.5610*sin_E1
                        + 0.1208*sin_E2     - 0.0642*sin_E3 	+ 0.0158*sin_E4
						+ 0.0252*sin(E5)	- 0.0066*sin_E6  	- 0.0047*sin(E7)
                        - 0.0046*sin(E8)	+ 0.0028*sin(E9)	+ 0.0052*sin_E10
                        + 0.0040*sin(E11)	+ 0.0019*sin(E12)	- 0.0044*sin_E13;

	alpha=DegToRad(alpha);
	beta=DegToRad(beta);
	om=DegToRad(om);

	double cos_om    = cos(om);
	double cos_alpha = cos(alpha);
	double cos_beta  = cos(beta);
	double sin_om    = sin(om);
    double sin_alpha = sin(alpha);
    double sin_beta  = sin(beta);


	double M[3][3]={0};

		M[0][0] = -sin_alpha*cos_om - cos_alpha*sin_beta*sin_om;
		M[1][0] = cos_alpha*cos_om - sin_alpha*sin_beta*sin_om;
        M[2][0] = sin_om*cos_beta;
        M[0][1] = sin_alpha*sin_om - cos_alpha*sin_beta*cos_om;
		M[1][1] = -cos_alpha*sin_om - sin_alpha*sin_beta*cos_om;
		M[2][1] = cos_om*cos_beta;
		M[0][2] = cos_alpha*cos_beta;
		M[1][2] = sin_alpha*cos_beta;
		M[2][2] = sin_beta;


	for(int i=0;i<6;i++) R_out[i]=0;
	for(int i=0;i<3;++i){
		for(int j=0;j<3;++j){
			R_out[i]   += R_in[j]   *M[j][i];
			R_out[i+3] += R_in[j+3] *M[j][i];
		}
	}
}
//---------------------------------------------------------------------
//процедура перевода из селенографической СК в селеноцентрическую СК
void SGtoSC(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]){
double R_in[6];
double R_out[6];
for(int i=0; i<3; i++){
	R_in[i]=r_in[i];
	R_in[i+3]=v_in[i];
}
SGtoSC(R_in,JD, R_out);
for(int i=0; i<3; i++){
	 r_out[i]=R_out[i];
	 v_out[i]=R_out[i+3];
}
}
//---------------------------------------------------------------------
//процедура перевода из селенографической СК в селеноцентрическую СК
 void SGtoSC(double R_in[6],double JD, double R_out[6]){

	double dt = JD + (dUTC+dUT1)/86400. - 2451545.;
	double T = dt/36525.;
	double E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13;
	double sin_E1,sin_E2,sin_E3,sin_E4,sin_E6,sin_E10,sin_E13;

    E1= DegToRad(125.045- 0.0529921*dt);
	E2= DegToRad(250.089- 0.1059842*dt);
	E3= DegToRad(260.008+13.0120009*dt);
	E4= DegToRad(176.625+13.3407154*dt);
	E5= DegToRad(357.529+ 0.9856003*dt);
	E6= DegToRad(311.589+26.4057084*dt);
	E7= DegToRad(134.963+13.0649930*dt);
	E8= DegToRad(276.617+ 0.3287146*dt);
	E9= DegToRad( 34.226+ 1.7484877*dt);
	E10=DegToRad( 15.134- 0.1589763*dt);
	E11=DegToRad(119.743+ 0.0036096*dt);
	E12=DegToRad(239.961+ 0.1643573*dt);
	E13=DegToRad( 25.053+12.9590088*dt);

	sin_E1  = sin(E1);
    sin_E2  = sin(E2);
    sin_E3  = sin(E3);
    sin_E4  = sin(E4);
	sin_E6  = sin(E6);
    sin_E10 = sin(E10);
    sin_E13 = sin(E13);

	double alpha = 269.9949 + 0.0031*T      - 3.8787*sin_E1	    - 0.1204*sin_E2
				            + 0.07*sin_E3	- 0.0172*sin_E4	    + 0.0072*sin_E6
                                            - 0.0052*sin_E10	+ 0.0043*sin_E13;

    double beta  = 66.5392  + 0.013*T	    + 1.5419*cos(E1)	+ 0.0239*cos(E2)
							- 0.0278*cos(E3)+ 0.0068*cos(E4)    - 0.0029*cos(E6)
				            + 0.0009*cos(E7)+ 0.0008*cos(E10)   - 0.0009*cos(E13);

    double om = 38.3213	+ 13.17635815*dt	- 1.4e-12*dt*dt		+ 3.5610*sin_E1
                        + 0.1208*sin_E2     - 0.0642*sin_E3 	+ 0.0158*sin_E4
						+ 0.0252*sin(E5)	- 0.0066*sin_E6  	- 0.0047*sin(E7)
                        - 0.0046*sin(E8)	+ 0.0028*sin(E9)	+ 0.0052*sin_E10
                        + 0.0040*sin(E11)	+ 0.0019*sin(E12)	- 0.0044*sin_E13;

    alpha=DegToRad(alpha);
	beta=DegToRad(beta);
    om=DegToRad(om);

    double cos_om    = cos(om);
    double cos_alpha = cos(alpha);
	double cos_beta  = cos(beta);
    double sin_om    = sin(om);
    double sin_alpha = sin(alpha);
    double sin_beta  = sin(beta);


    double M[3][3]={0};

        M[0][0] = -sin_alpha*cos_om - cos_alpha*sin_beta*sin_om;
        M[0][1] = cos_alpha*cos_om - sin_alpha*sin_beta*sin_om;
		M[0][2] = sin_om*cos_beta;
        M[1][0] = sin_alpha*sin_om - cos_alpha*sin_beta*cos_om;
        M[1][1] = -cos_alpha*sin_om - sin_alpha*sin_beta*cos_om;
        M[1][2] = cos_om*cos_beta;
        M[2][0] = cos_alpha*cos_beta;
		M[2][1] = sin_alpha*cos_beta;
		M[2][2] = sin_beta;


	for(int i=0;i<6;i++) R_out[i]=0;
	for(int i=0;i<3;++i){
		for(int j=0;j<3;++j){
			R_out[i]   += R_in[j]   *M[j][i];
			R_out[i+3] += R_in[j+3] *M[j][i];
		}
	}
}












