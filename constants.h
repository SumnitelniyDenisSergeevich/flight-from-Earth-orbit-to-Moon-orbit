#pragma once
//---------------------------------------------------------------------------
#include <array>
#include <map>
#include <string>
#include "vector.h"
extern const double pi;     // Pi
extern const double pi_2;   // 2 Pi
extern const double AE;     // ��������. ������� (��)
extern const double Rz;		// ���. ������ ����� (��)
extern const double R_Sun;  // ������ ������ (��)
extern const double R_Moon; // ������ ���� (��)
extern const double mu;    	// �������������� ���������� �����
extern const double mum;    // �������������� ���������� ����
extern const double mus;    // �������������� ���������� ������
extern const double rtd;
extern const double dtr;
/*������ �������������� ���������� ������.
���� ��������� ������ ������ ��� ����������,
������������� ������ � ������� �� ������� ������
 1-��������, 2-������, 3-�����, 4-����, 5-������,
 6-������, 7-����, 8-������, 9-������, 10-����, 11-������*/
extern const double mu_planet[11];



extern const double Alpha_e;    //������ �����
extern const double REb; 		//����� ������� ������� ����������
extern const double w_earth;    //������� �������� �������� �����, 1/�
extern const double MP_lat;     //�������� ���������� ������
extern const double MP_lon;    	//������� ���������� ������
extern const double CONST_P0;	//�������� ���������� ����� �� ������ �����
extern const double CONST_A_EARTH; //������� ������� ������ ����� � ��.


extern const int    N_garm;  	// ���� ������� ���������� ��������������� ���������� �� ����������� ����������
extern const int    Ks;  		// ������� ����� ����. ���������� ������ ��� �������������� ���. ��-� ��������
extern const int    Km;  		// ������� ����� ����. ���������� ���� ��� �������������� ���. ��-� ��������

extern const double dUTC; 		// ��������� �������� ����� ������� � ��������� ���������������� ��������
extern const double dUT1;  		//��������� �������� ����� ������ �������� (TT) � ��������� �������� (UT).


extern const double rkf8_alpha[13];
extern const double rkf8_ce[13];
extern const double rkf8_ce_cap[13];
extern const double rkf8_beta[13][12];

//������������ ��� ���������� �������������� ������� �������-������
extern const double dp5_alpha[7];
extern const double dp5_beta[7][6];
extern const double dp5_ce[7];
extern const double dp5_ce_cap[7];

/*������������ ���� ���������� ��������������� ��������������� ���� �����
�� 32 ��������� */
//extern const double Alpha[13];
//extern const double Ck[13];
//extern const double Cks[13];
//extern const double Bk[13][12];

/*������������ ���� ���������� ��������������� ��������������� ���� �����
�� 32 ��������� */
extern const std::array<std::array<double, 33>, 33> Cnn;
extern const std::array<std::array<double, 33>, 33> Snn;
extern const std::array<std::array<double, 76>, 76> Snn_moon;
extern const std::array<std::array<double, 76>, 76> Cnn_moon;

//������������ ��� ����� ��������� � ������� ->
extern const double coef_a1[106];
extern const double coef_a2[106];
extern const double coef_a3[106];
extern const double coef_a4[106];
extern const double coef_a5[106];
extern const double coef_A[106];
extern const double coef_B[106];
extern const double coef_C[106];
extern const double coef_D[106];
//  <- ������������ ��� ����� ��������� � �������

//������������ ��� ����������� ��������� ��������� ->
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
//  <- ������������ ��� ����������� ��������� ���������
extern const std::map<std::string, Vector> CONST_NIP_COOR;
