#ifndef DEPHEM_HELP_HPP
#define DEPHEM_HELP_HPP

namespace dph
{

	// ************************************************************************** //
	//                                   Body                                     //
	//                                                                            //
	//                                ������� ���                                 //
	// -------------------------------------------------------------------------- //
	//                                 ��������                                   //
	// -------------------------------------------------------------------------- //
	// ��������������� �����, �������� �������� ���������� ��� ������             //
	// dph::EphemerisReelase::calculateBody(...).                                 //
	//                                                                            //
	// �������� �������� �������� ������������� �������� ���, ��� �������         //
	// ����� �������� ��������� ����������.                                       //
	//                                                                            //
	// ************************************************************************** //
	class Body
	{
	public:

		static const unsigned MERCURY = 1;
		static const unsigned VENUS = 2;
		static const unsigned EARTH = 3;
		static const unsigned MARS = 4;
		static const unsigned JUPITER = 5;
		static const unsigned SATURN = 6;
		static const unsigned URANUS = 7;
		static const unsigned NEPTUNE = 8;
		static const unsigned PLUTO = 9;
		static const unsigned MOON = 10;
		static const unsigned SUN = 11;
		static const unsigned SSBARY = 12;	// ��������� ��������� �������.
		static const unsigned EMBARY = 13;	// ��������� ������� �����-����.

	private:
		Body(); // ������ �� �������� ������� ���� Body.
	};

	// ************************************************************************** //
	//                                   Other                                    //
	//                                                                            //
	//                         ������� ������ ���������                           //
	// -------------------------------------------------------------------------- //
	//                                 ��������                                   //
	// -------------------------------------------------------------------------- //
	// ��������������� �����, �������� �������� ���������� ��� ������             //
	// dph::EphemerisReelase::calculateOther(...).                                //
	//                                                                            //
	// �������� �������� �������� ������������� �������� ������ ���������, ���    //
	// ������� ����� �������� ��������� ����������.                               //
	//                                                                            //
	// ************************************************************************** //
	class Other
	{
	public:

		static const unsigned EARTH_NUTATIONS = 14;
		static const unsigned LUNAR_MANTLE_LIBRATION = 15;
		static const unsigned LUNAR_MANTLE_ANGULAR_VELOCITY = 16;
		static const unsigned TTmTDB = 17;
	private:
		Other(); // ������ �� �������� ������� ���� Other.
	};

	// ************************************************************************** //
	//                                 Calculate                                  //
	//                                                                            //
	//                      ������� ����������� ����������                        //
	// -------------------------------------------------------------------------- //
	//                                 ��������                                   //
	// -------------------------------------------------------------------------- //
	// ��������������� �����, �������� �������� ���������� ��� �������:           //
	//    dph::EphemerisReelase::calculateBody(...),                              //
	//    dph::EphemerisReelase::calculateOther(...).                             //
	//                                                                            //
	// �������� �������� �������� ������������� �������� ����������� ����������,  // 
	// ������� ����� ��������.                                                    //
	//                                                                            //
	// ************************************************************************** //
	class Calculate
	{
	public:

		static const unsigned POSITION = 0;
		static const unsigned STATE = 1;

	private:
		Calculate(); // ������ �� �������� ������� ���� Calculate.
	};

} // namespace dph

#endif // DEPHEM_HELP_HPP