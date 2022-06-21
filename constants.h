#pragma once
//---------------------------------------------------------------------------
#include <array>
#include <map>
#include <string>
#include "vector.h"
extern const double pi;     // Pi
extern const double pi_2;   // 2 Pi
extern const double AE;     // Астроном. единица (км)
extern const double Rz;		// Экв. радиус Земли (км)
extern const double R_Sun;  // радиус Солнца (км)
extern const double R_Moon; // радиус Луны (км)
extern const double mu;    	// гравитационная постоянная Земли
extern const double mum;    // гравитационная постоянная Луны
extern const double mus;    // гравитационная постоянная Солнца
extern const double rtd;
extern const double dtr;
/*Массив гравитационных постоянных планет.
ниже приведены номера планет для планетария,
сответственно номера в массиве на единицу меньше
 1-Меркурий, 2-Венера, 3-Земля, 4-Марс, 5-Юпитер,
 6-Сатурн, 7-Уран, 8-Нептун, 9-Плутон, 10-Луна, 11-Солнце*/
extern const double mu_planet[11];



extern const double Alpha_e;    //Сжатие Земли
extern const double REb; 		//Малая полуось земного эллипсоида
extern const double w_earth;    //Угловая скорость вращения Земли, 1/с
extern const double MP_lat;     //Коширота магнитного полюса
extern const double MP_lon;    	//Долгота магнитного полюса
extern const double CONST_P0;	//давление солнечного света на орбите Земли
extern const double CONST_A_EARTH; //Большая полуось орбиты Земли в км.


extern const int    N_garm;  	// ранг матрицы разложения гравитационного потенциала по сферическим гармоникам
extern const int    Ks;  		// признак учета грав. возмущения Солнца при интегрировании диф. ур-й движения
extern const int    Km;  		// признак учета грав. возмущения Луны при интегрировании диф. ур-й движения

extern const double dUTC; 		// временная поправка между атомным и всемирным координированным временем
extern const double dUT1;  		//временная поправка между земным временем (TT) и всемирным временем (UT).


extern const double rkf8_alpha[13];
extern const double rkf8_ce[13];
extern const double rkf8_ce_cap[13];
extern const double rkf8_beta[13][12];

//Коэффициенты для численного интегрирования методом Дормана-Принса
extern const double dp5_alpha[7];
extern const double dp5_beta[7][6];
extern const double dp5_ce[7];
extern const double dp5_ce_cap[7];

/*Коэффициенты ряда разложения нецентральности гравитационного поля Земли
до 32 гормоники */
//extern const double Alpha[13];
//extern const double Ck[13];
//extern const double Cks[13];
//extern const double Bk[13][12];

/*Коэффициенты ряда разложения нецентральности гравитационного поля Земли
до 32 гормоники */
extern const std::array<std::array<double, 33>, 33> Cnn;
extern const std::array<std::array<double, 33>, 33> Snn;
extern const std::array<std::array<double, 76>, 76> Snn_moon;
extern const std::array<std::array<double, 76>, 76> Cnn_moon;

//коэффициенты для учета прицессии и нутации ->
extern const double coef_a1[106];
extern const double coef_a2[106];
extern const double coef_a3[106];
extern const double coef_a4[106];
extern const double coef_a5[106];
extern const double coef_A[106];
extern const double coef_B[106];
extern const double coef_C[106];
extern const double coef_D[106];
//  <- коэффициенты для учета прицессии и нутации

//коэффициенты для определения плотности атмосферы ->
extern const double atm_k120[4][4];
extern const double atm_F0[7];
extern const double atm_a[2][8][7];
extern const double atm_b[2][6][7];
extern const double atm_c[2][6][7];
extern const double atm_n[2][4][7];
extern const double atm_fi[2][2][7];
extern const double atm_d[2][6][7];
extern const double atm_e[2][14][7];
extern const double atm_l[2][6][7];
extern const double atm_Ad[9];
//  <- коэффициенты для определения плотности атмосферы
extern const std::map<std::string, Vector> CONST_NIP_COOR;
