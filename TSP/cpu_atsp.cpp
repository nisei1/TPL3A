/************************************************
 *    ��Ώ̏���Z�[���X�}�����v���O����      *
 *�@�@City-block exchange + block insertion     *
 *                 CPU Version                  *
 *  �@            Version 0.10                  *
 ************************************************/

/* SFMT (SSE2�Ή���)
InitMt(S)				�����̎�r�ɂ�鏉����
InitMtEx(K,L)			�����k�̔z��j�ɂ�鏉����
NextMt()				32�r�b�g�����Ȃ������̗���
NextUnif()				�O�ȏ�P�����̗���(53bit���x)
NextInt(N)				�O�ȏ�m�����̐�������
NextIntEx(N)			�ۂߌ덷�̂Ȃ��O�ȏ�m�����̐�������

NextChisq(N)			���R�x�˂̃J�C�Q�敪�z
NextGamma(A)			�p�����[�^�`�̃K���}���z
NextGeometric(P)		�m���o�̊􉽕��z
NextTriangle()			�O�p���z
NextExp()				���ςP�̎w�����z
NextNormal()			�W�����K���z(�ő�8.57��)
NextUnitVect(V,N)		�m�����̃����_���P�ʃx�N�g��
NextBinomial(N,P)		�p�����[�^�m,�o�̂Q�����z
NextBinormal(R,&X,&Y)	���֌W���q�̂Q�ϗʐ��K���z
NextBeta(A,B)			�p�����[�^�`,�a�̃x�[�^���z
NextPower(N)			�p�����[�^�m�̗ݏ敪�z
NextLogistic()			���W�X�e�B�b�N���z
NextCauchy()			�R�[�V�[���z
NextFDist(A,B)			���R�x�`,�a�̂e���z
NextPoisson(L)			���σɂ̃|�A�\�����z
NextTDist(N)			���R�x�m�̂����z
NextWeibull(A)			�p�����[�^���̃��C�u�����z
*/

#pragma warning(disable : 4995)

/****************************************
 *                                      *
 *       �w�b�_�[�t�@�C���̐錾         *
 *                                      *
 ****************************************/

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "element_operate.h"

#include "zsfmt.c"

/****************************************
 *                                      *
 *      �@�@�@ �萔�̐錾�@�@�@         *
 *                                      *
 ****************************************/

#define InName "rbg443a.txt"
#define OutName "rbg443a_z1_tour_cpu03.csv"
#define CITYNO 443
#define SEARCH_COUNT 100
#define REPEAT 10	 //���s��
#define LEARN_MODE 0 //�w�K���ԁF0=�����_���A9=����
#define TOUR_MODE 5	 // 0=�����_���A5=�ŒZ��p

#define PATCHING 0

#define BCESIZE 2
#define BIMSIZE 2

#define MSGMODE 1	  //���b�Z�[�W�̉�ʏo�͂��邩�H
#define TOURSEARCH 0  //�K�⏇���Ō������邩�H�j���[�����ԍ��Ō������邩�H \
					  //1��2-opt�K�⏇���C0�ŗ����j���[�����ԍ�
#define TWOOPTRESET 1 //1�ŏ��œK���ɒu������
#define BlockMinTourReset 1

#define TIMEEXP 0 //1�ɂ���ƌv�Z���Ԃ��v��

#define ALPHA 1.0
#define BETA 75.0
#define THETA 0.5
#define KR 0.99
#define R THETA *(1.0 - KR)
#define W 0.15
#define EPSILON 0.002

#define CNG_ALPHA 3.0
#define CNG_A 2.5

using namespace std;

/****************************************
 *                                      *
 *        �@�@ �}�N���̐錾�@�@�@�@     *
 *                                      *
 ****************************************/

/* SSE2�e�N�m���W�[���g�p���� */
#ifdef HAVE_SSE2
DefaultMt->sse2;
#endif

#ifndef MEXP
#define MEXP 19937L
#endif

/* 16bit �R���p�C���̏ꍇ�A�������ɐ��������邽�ߌn��̐������炷 */
#if INT_MAX == 32767 && MEXP == 11213L
#define MSEQ2 22
#elif INT_MAX == 32767 && MEXP == 19937L
#define MSEQ2 14
#elif INT_MAX == 32767 && MEXP == 44497L
#define MSEQ 6
#elif INT_MAX == 32767 && MEXP == 86243L
#define MSEQ 3
#elif INT_MAX == 32767 && MEXP == 132049L
#define MSEQ 2
#elif INT_MAX == 32767 && MEXP == 216091L
#define MSEQ 1
#elif MEXP <= 19937L
#define MSEQ2 32
#else
#define MSEQ 32
#endif

unsigned long Loop(void)
{
	long i;
	unsigned long ret = 0;

	for (i = 0; i < 100000000L; i++) /*ret+=*/
		NextMt();
	return ret;
}

#ifndef MSEQ
#define MSEQ 1
#endif

#ifdef LSI_C
#undef CLOCKS_PER_SEC
#endif

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000
#include <dos.h>
static long clock(void)
{
	union REGS t;
	t.h.ah = 0x2c;
	intdos(&t, &t);
	return t.h.dl * 10L + t.h.dh * 1000L + t.h.cl * 60000L + t.h.ch * 3600000L;
}
#endif

/****************************************
 *                                      *
 *        �N���X�E�\���̂̐錾�@�@�@    *
 *                                      *
 ****************************************/

class CityData
{
public:
	int CityX[CITYNO];
	int CityY[CITYNO];
	int CityDistance[CITYNO * CITYNO];
	int Tour[CITYNO + 1];
	int MinTour[CITYNO + 1];
};

class ChaosNN
{
public:
	double xi[CITYNO * SEARCH_COUNT];
	double zeta[CITYNO * SEARCH_COUNT];
	double eta[CITYNO * SEARCH_COUNT];
	double x[CITYNO * SEARCH_COUNT];
};

class ChangeChaosNN
{
public:
	double xi[SEARCH_COUNT];
	double zeta[SEARCH_COUNT];
	double eta[SEARCH_COUNT];
	double x[SEARCH_COUNT];
	int BeforeTour;
	int AfterTour;
};

/****************************************
 *                                      *
 *         �O���[�o���ϐ��̐錾�@�@     *
 *                                      *
 ****************************************/
CityData CD;
ChaosNN CNN;
ChangeChaosNN CNG;

bool typeexchange = true;

int twoopttimes = 0;
int blocktimes = 0;

int temp_j;

int h_Tour[CITYNO * CITYNO];

ifstream in(InName);
ofstream out(OutName);

/****************************************
 *                                      *
 *        �@�@ �@�֐��̐錾�@�@�@�@     *
 *                                      *
 ****************************************/
int calc_distance(int *, int); //����H�̑���p�v�Z
void make_tour(int);		   //��������H�쐬
void initialize(void);		   //�ϐ��̏������ƒn�_�f�[�^�Ǎ�
void patching(void);
void bce_calc_distance(int, int, int *); //2-opt�����Ƃ��̏���H�̑���p�v�Z
void bim_calc_distance(int, int, int *); //2-opt�����Ƃ��̏���H�̑���p�v�Z
double bce_max(int, int);				 //xi�v�Z��max_j�̌v�Z
double bim_max(int, int);				 //xi�v�Z��max_j�̌v�Z
double sum_eta(int, int);				 //eta�v�Z��sum_k�̌v�Z
double sigmoid(double);					 //�V�O���C�h�֐�
void bce_chaosNN(int, int);				 //�u���b�N�V�t�g�����J�I�X�j���[�����l�b�g���[�N
void bim_chaosNN(int, int);				 //�u���b�N�V�t�g�����J�I�X�j���[�����l�b�g���[�N
void bce_changeTour(int, int);			 //���ۂɑI�����ꂽ�n�_���������{
void bim_changeTour(int, int);			 //���ۂɑI�����ꂽ�n�_���������{
void ChangeChaosNeuralNet(int);

int main(int argc, char **argv)
{
	int i;
	int length = 0, rlength = 0;
	int Temp[CITYNO];
	clock_t start, end;
	/*
	try{
		CD.CityX=new int[CITYNO];
		CD.CityY=new int[CITYNO];
		CD.CityDistance=new int[CITYNO*CITYNO];

		CNN.xi=new double[CITYNO*SEARCH_COUNT];
		CNN.zeta=new double[CITYNO*SEARCH_COUNT];
		CNN.eta=new double[CITYNO*SEARCH_COUNT];
		CNN.x=new double[CITYNO*SEARCH_COUNT];

		CNG.xi=new double[SEARCH_COUNT];
		CNG.zeta=new double[SEARCH_COUNT];
		CNG.eta=new double[SEARCH_COUNT];
		CNG.x=new double[SEARCH_COUNT];
	}
	catch(bad_alloc){
	    cout << "error\n"; 
        abort(); // �I�������� 
    }
*/
	cout << InName << endl;
	out << InName << endl;
	out << "ALPHA=" << ALPHA << " BETA=" << BETA << " THETA=" << THETA << " KR=" << KR << " R=" << R << " W=" << W << " EPSILON=" << EPSILON << endl;
	out << "CNG_ALPHA=" << CNG_ALPHA << " A=" << CNG_A << endl;

	while (!in.eof())
	{
		for (i = 0; i < CITYNO; i++)
		{
			for (int j = 0; j < CITYNO; j++)
			{
				CD.CityX[j] = 0; //�@������
				CD.CityY[j] = 0; //�@������
				in >> CD.CityX[j] >> CD.CityY[j] >> CD.CityDistance[i * CITYNO + j];
			}
		}
	}

	for (i = 0; i < REPEAT; i++)
	{
		if (TIMEEXP == 1)
		{
			start = clock();
		}

		initialize();

		cout << "Initial Tour Cost:" << calc_distance(CD.Tour, 0) << endl;
		CD.MinTour[CITYNO] = calc_distance(CD.Tour, 0);

		for (int t = 0; t < SEARCH_COUNT; t++)
		{

			CNG.BeforeTour = CD.MinTour[CITYNO];
			CNG.AfterTour = 100000;

			if (CNG.x[t] > THETA)
			{
				typeexchange = !typeexchange;
			}

			if (typeexchange)
			{ // �u���b�N�V�t�g�@�����{����
				for (int p = 0; p < CITYNO; p++)
				{

					bce_chaosNN(CD.Tour[p], t);

					if (CNN.x[CITYNO * t + CD.Tour[p]] > THETA)
					{

						/*				���ۂ̃c�A�[�̑g�ݑւ� */
						bce_changeTour(CD.Tour[p], temp_j);

						length = 0;
						rlength = 0;

						/*						������ */
						length = calc_distance(CD.Tour, 0);

						/*						�t���� */
						for (int x = CITYNO - 1; x > 0; x--)
						{
							rlength += (int)CD.CityDistance[CITYNO * CD.Tour[x % CITYNO] + CD.Tour[(x - 1) % CITYNO]];
						}

						/*						�t�����̔�p���������Ƃ��̃c�A�[�̔��] */
						if (length > rlength)
						{
							for (int x = 0; x < CITYNO; x++)
							{
								Temp[x] = CD.Tour[x];
							}
							for (int x = 0; x < CITYNO; x++)
							{
								CD.Tour[CITYNO - x - 1] = Temp[x];
							}
							if (CD.MinTour[CITYNO] > rlength)
							{
								CD.MinTour[CITYNO] = (int)rlength;
							}
						}
					}
					else
					{
						CNG.AfterTour = CNG.BeforeTour;
					}
				}
			}
			else
			{
				for (int p = 0; p < CITYNO; p++)
				{

					bim_chaosNN(CD.Tour[p], t);

					if (CNN.x[CITYNO * t + CD.Tour[p]] > THETA)
					{

						/*				���ۂ̃c�A�[�̑g�ݑւ� */
						bim_changeTour(CD.Tour[p], temp_j);

						length = 0;
						rlength = 0;

						/*						������ */
						length = calc_distance(CD.Tour, 0);

						/*						�t���� */
						for (int x = CITYNO - 1; x > 0; x--)
						{
							rlength += (int)CD.CityDistance[CITYNO * CD.Tour[x % CITYNO] + CD.Tour[(x - 1) % CITYNO]];
						}

						/*						�t�����̔�p���������Ƃ��̃c�A�[�̔��] */
						if (length > rlength)
						{
							for (int x = 0; x < CITYNO; x++)
							{
								Temp[x] = CD.Tour[x];
							}
							for (int x = 0; x < CITYNO; x++)
							{
								CD.Tour[CITYNO - x - 1] = Temp[x];
							}
						}
					}
					else
					{
						CNG.AfterTour = CNG.BeforeTour;
					}
				}
			}
			ChangeChaosNeuralNet(t);
		}

		if (TIMEEXP == 1)
		{
			end = clock();
		}

		if (TIMEEXP == 0)
		{
			cout << endl
				 << endl
				 << "Final Tour: " << CD.MinTour[CITYNO] << endl
				 << endl
				 << endl;
			//			out << endl << endl << "Final Tour:" << CD.MinTour[CITYNO] << endl << endl << endl;
			out << CD.MinTour[CITYNO] << endl;
		}
		else
		{
			cout << endl
				 << endl
				 << "Processing Time: " << end - start << "(ms)" << endl
				 << endl
				 << endl;
			out << end - start << endl;
		}
	}
	/*
	delete [] CD.CityX;
	delete [] CD.CityY;
	delete [] CD.CityDistance;

	delete [] CNN.eta;
	delete [] CNN.zeta;
	delete [] CNN.x;
	delete [] CNN.xi;
	
	delete [] CNG.eta;
	delete [] CNG.zeta;
	delete [] CNG.x;
	delete [] CNG.xi;
*/
	return 0;
}

int calc_distance(int *Tour, int pitch)
{
	int length;

	length = 0;

	for (int n = 0; n < CITYNO; n++)
	{
		length += (int)CD.CityDistance[Tour[pitch + n % CITYNO] * CITYNO + Tour[pitch + (n + 1) % CITYNO]];
	}
	//   cout << "Tour;" << length << endl;
	return length;
}

void make_tour(int i)
{
	int n;
	int distance = 1000000;
	int *lrn_neuron;
	int j;
	vector<int> vec;

	lrn_neuron = new int[CITYNO];

	if (vec.empty() == false)
	{
		vec.clear();
	}
	distance = 1000000;

	if (i != 0)
	{
		for (int p = 0; p < CITYNO; p++)
		{
			lrn_neuron[p] = 1;
		}

		n = NextMt() % CITYNO;
		CD.Tour[0] = n;
		lrn_neuron[n] = 0;

		for (int p = 0; p < CITYNO; p++)
		{
			for (int i = 0; i < CITYNO; i++)
			{
				if (distance > CD.CityDistance[CD.Tour[p] * CITYNO + i] && lrn_neuron[i] == 1)
				{
					distance = CD.CityDistance[CD.Tour[p] * CITYNO + i];
					j = i;
				}
			}
			lrn_neuron[j] = 0;
			CD.Tour[p + 1] = j;
			CD.MinTour[p + 1] = j;
			distance = 1000000;
		}
	}
	else
	{
		if (vec.empty() == false)
		{
			vec.clear();
		}

		for (i = 0; i < CITYNO; i++)
		{
			vec.push_back(i);
		}

		random_shuffle(vec.begin(), vec.end());

		for (i = 0; i < CITYNO; i++)
		{
			CD.Tour[i] = vec[i];
			CD.MinTour[i] = vec[i];
		}
	}
	delete lrn_neuron;
}

void initialize()
{
	int j;

	InitMt((unsigned long)time(NULL)); //�����̏�����
	srand((unsigned long)time(NULL));

	if (PATCHING == 0)
	{
		make_tour(TOUR_MODE);
	}
	else
	{
		patching();
	}

	for (j = 0; j < CITYNO; j++)
	{
		CNN.xi[j] = NextUnif();
		CNN.zeta[j] = NextUnif();
		CNN.eta[j] = NextUnif();
		CNN.x[j] = 0.0;
	}

	CNG.xi[0] = NextUnif();
	CNG.zeta[0] = NextUnif();
	CNG.eta[0] = NextUnif();
	CNG.x[0] = NextUnif();
}

void patching(void)
{
	vector<int> c1, c2, length;
	int current_tour, min_tour;

	for (int n = 0; n < CITYNO; n++)
	{
		c2.push_back(n);
	}

	for (unsigned int m = 0; m < c2.size(); m++)
	{
		for (unsigned int n = 0; n < c2.size(); n++)
		{
			if (c2.at(m) == c2.at(n))
			{
				length.push_back(999999999);
			}
			else
			{
				length.push_back(CD.CityDistance[c2.at(m) * CITYNO + c2.at(n)] + CD.CityDistance[c2.at(m) * CITYNO + c2.at(n)]);
			}
		}
	}
	min_tour = 0;
	min_tour = find(length.begin(), length.end(), *min_element(length.begin(), length.end())) - length.begin();
	c1.push_back(c2.at(min_tour / c2.size()));
	c1.push_back(c2.at(min_tour % c2.size()));

	for (unsigned int n = 0; n < c1.size(); n++)
	{
		c2.erase(find(c2.begin(), c2.end(), c1.at(n)));
	}
	random_shuffle(c2.begin(), c2.end());

	current_tour = 0;
	current_tour = CD.CityDistance[c1.front() + CITYNO + c1.back()] + CD.CityDistance[c1.back() * CITYNO + c1.front()];

	while (!c2.empty())
	{
		length.clear();

		for (unsigned int m = 0; m < c1.size(); m++)
		{
			for (unsigned int n = 0; n < c2.size(); n++)
			{
				length.push_back(current_tour - CD.CityDistance[c1.at(m % c1.size()) * CITYNO + c1.at((m + 1) % c1.size())] + CD.CityDistance[c1.at(m % c1.size()) * CITYNO + c2.at(n)] + CD.CityDistance[c2.at(n) * CITYNO + c1.at((m + 1) % c1.size())]);
			}
		}
		min_tour = 0;
		min_tour = find(length.begin(), length.end(), *min_element(length.begin(), length.end())) - length.begin();
		current_tour = length.at(min_tour);

		c1.insert(c1.begin() + (min_tour / c2.size()) + 1, c2.at(min_tour % c2.size()));
		c2.erase(find(c2.begin(), c2.end(), c2.at(min_tour % c2.size())));
	}

	for (unsigned int i = 0; i < c1.size(); i++)
	{
		CD.Tour[i] = c1.at(i);
		CD.MinTour[i] = c1.at(i);
	}
}

void bce_calc_distance(int i, int j, int *d_Tour)
{
	int pitch_x;

	pitch_x = j * CITYNO;

	if (i < j)
	{
		element_swap(d_Tour, pitch_x + i, pitch_x + j);
		element_move(d_Tour, CITYNO, pitch_x + (i + 1) % CITYNO, pitch_x + j);
	}
	else if (i > j && i != CITYNO - 1)
	{
		element_swap(d_Tour, pitch_x + i, pitch_x + j);
		element_move(d_Tour, CITYNO, pitch_x + (i + 1) % CITYNO, pitch_x + (j + 1) % CITYNO);
	}
	else if (i > j && i == CITYNO - 1)
	{
		element_swap(d_Tour, pitch_x + i, pitch_x + j);
		element_move(d_Tour, CITYNO, pitch_x + i, pitch_x + j);
	}
}

void bim_calc_distance(int i, int j, int *d_Tour)
{
	int pitch_x;

	pitch_x = j * CITYNO;

	for (int n = 0; n < BIMSIZE; n++)
	{
		if (i < j)
		{
			for (int x = i; x < j; x++)
			{
				element_move(d_Tour, CITYNO, pitch_x + x % CITYNO, pitch_x + (x + 1) % CITYNO);
			}
		}
		else if (i > j)
		{
			for (int x = i; x < j + (CITYNO - 1); x++)
			{
				element_move(d_Tour, CITYNO, pitch_x + x % CITYNO, pitch_x + (x + 1) % CITYNO);
			}
		}
	}
}

double bce_max(int i, int t)
{
	double d0;
	double temp_max;
	int tour_length;

	temp_max = -999999.0;
	d0 = (double)calc_distance(CD.Tour, 0);

	for (int m = 0; m < CITYNO; m++)
	{
		for (int n = 0; n < CITYNO; n++)
		{
			h_Tour[m * CITYNO + n] = CD.Tour[n];
		}
	}

	for (int n = 0; n < CITYNO; n++)
	{
		bce_calc_distance(i, n, h_Tour);
	}

	/*
	for(int iz=0;iz<CITYNO;iz++){
	//	cout << min_length[iz] << endl;
		for(int jz=0;jz<CITYNO;jz++){
			cout << h_Tour[iz*CITYNO+jz] << " ";
		}
		cout << endl;
	}
	cout << endl;
*/
	//		out << d0 << " " << i << " " << j << " " << Calc_BlockShift_Distance(i,j) << endl;
	for (int j = 0; j < CITYNO; j++)
	{
		tour_length = calc_distance(h_Tour, j * CITYNO);
		//			cout << tour_length << "\t";
		if (temp_max < CNN.zeta[(t + 1) * CITYNO + j] + BETA * (d0 - tour_length))
		{ //&& (i!=j || j!=(i+1)%CITYNO)){
			temp_max = CNN.zeta[(t + 1) * CITYNO + j] + BETA * (d0 - tour_length);
			temp_j = j;
		}
	}
	//	cout << endl;
	return temp_max;
}

double bim_max(int i, int t)
{
	double d0;
	double temp_max;
	int tour_length;

	temp_max = -999999.0;
	d0 = (double)calc_distance(CD.Tour, 0);

	for (int m = 0; m < CITYNO; m++)
	{
		for (int n = 0; n < CITYNO; n++)
		{
			h_Tour[m * CITYNO + n] = CD.Tour[n];
		}
	}

	for (int n = 0; n < CITYNO; n++)
	{
		bim_calc_distance(i, n, h_Tour);
	}

	/*
	for(int iz=0;iz<CITYNO;iz++){
	//	cout << min_length[iz] << endl;
		for(int jz=0;jz<CITYNO;jz++){
			cout << h_Tour[iz*CITYNO+jz] << " ";
		}
		cout << endl;
	}
	cout << endl;
*/
	//		out << d0 << " " << i << " " << j << " " << Calc_BlockShift_Distance(i,j) << endl;
	for (int j = 0; j < CITYNO; j++)
	{
		tour_length = calc_distance(h_Tour, j * CITYNO);
		//			cout << tour_length << "\t";
		if (temp_max < CNN.zeta[(t + 1) * CITYNO + j] + BETA * (d0 - tour_length))
		{ //&& (i!=j || j!=(i+1)%CITYNO)){
			temp_max = CNN.zeta[(t + 1) * CITYNO + j] + BETA * (d0 - tour_length);
			temp_j = j;
		}
	}
	//	cout << endl;
	return temp_max;
}

double sum_eta(int i, int t)
{
	double sumX;

	sumX = 0.0;

	for (int k = 0; k < CITYNO; k++)
	{
		if (k != i)
		{
			sumX += CNN.x[t * CITYNO + k];
		}
	}
	return sumX;
}

double sigmoid(double x)
{
	return 1.0 / (1.0 + exp(-x / EPSILON));
}

void bce_chaosNN(int i, int t)
{
	//	t=1.0/(1.0+exp(-i/EPSILON));
	CNN.zeta[(t + 1) * CITYNO + i] = KR * CNN.zeta[t * CITYNO + i] - ALPHA * CNN.x[t * CITYNO + i] + R;

	CNN.xi[(t + 1) * CITYNO + i] = bce_max(i, t);
	CNN.eta[(t + 1) * CITYNO + i] = -W * sum_eta(i, t) + W;

	CNN.x[(t + 1) * CITYNO + i] = sigmoid(CNN.xi[(t + 1) * CITYNO + i] + CNN.eta[(t + 1) * CITYNO + i] + CNN.zeta[(t + 1) * CITYNO + i]);
}

void bce_changeTour(int i, int j)
{
	int length = 0;

	for (int n = 0; n < CITYNO; n++)
	{
		CD.Tour[n] = h_Tour[j * CITYNO + n];
		//		cout << h_Tour[j*CITYNO+n] << "\t";
	}
	//	cout << endl;

	/*
	//i���u���b�N�̐擪�Cj�͌�������
	if(i<j){
		element_swap(CD.Tour,i,j);
		element_move(CD.Tour,CITYNO,(i+1),j);
	}
	else if(i>j && i!=CITYNO-1){
		element_swap(CD.Tour,i,j);
		element_move(CD.Tour,CITYNO,(i+1),(j+1));
	}
	else if(i>j && i==CITYNO-1){
		element_swap(CD.Tour,i,j);
		element_move(CD.Tour,CITYNO,i,j);
	}
*/
	length = calc_distance(CD.Tour, 0);
	CNG.AfterTour = length;

	if (CD.MinTour[CITYNO] >= length)
	{
		if (CD.MinTour[CITYNO] > length && MSGMODE == 1)
		{
			cout << BCESIZE << "BCE Tour:" << length << endl;
		}
		//		out << "Tour:" << length << endl;
		for (int n = 0; n < CITYNO; n++)
		{
			CD.MinTour[n] = CD.Tour[n];
		}
		CD.MinTour[CITYNO] = (int)length;
	}
	else
	{
		if (BlockMinTourReset == 1)
		{
			for (int n = 0; n < CITYNO; n++)
			{
				CD.Tour[n] = CD.MinTour[n];
			}
		}
	}
}

void bim_chaosNN(int i, int t)
{
	//	t=1.0/(1.0+exp(-i/EPSILON));
	CNN.zeta[(t + 1) * CITYNO + i] = KR * CNN.zeta[(t + 1) * CITYNO + i] - ALPHA * CNN.x[(t + 1) * CITYNO + i] + R;

	CNN.xi[(t + 1) * CITYNO + i] = bim_max(i, t);
	CNN.eta[(t + 1) * CITYNO + i] = -W * sum_eta(i, t) + W;

	CNN.x[(t + 1) * CITYNO + i] = sigmoid(CNN.xi[(t + 1) * CITYNO + i] + CNN.eta[(t + 1) * CITYNO + i] + CNN.zeta[(t + 1) * CITYNO + i]);
}

void bim_changeTour(int i, int j)
{
	int length = 0;

	for (int n = 0; n < CITYNO; n++)
	{
		CD.Tour[n] = h_Tour[j * CITYNO + n];
		//		cout << h_Tour[j*CITYNO+n] << "\t";
	}
	//	cout << endl;
	/*

	//i���u���b�N�̐擪�Cj�͌�������
	if(i<j){
		element_swap(CD.Tour,i,j);
		element_move(CD.Tour,CITYNO,(i+1),j);
	}
	else if(i>j && i!=CITYNO-1){
		element_swap(CD.Tour,i,j);
		element_move(CD.Tour,CITYNO,(i+1),(j+1));
	}
	else if(i>j && i==CITYNO-1){
		element_swap(CD.Tour,i,j);
		element_move(CD.Tour,CITYNO,i,j);
	}
*/
	length = calc_distance(CD.Tour, 0);
	CNG.AfterTour = length;

	if (CD.MinTour[CITYNO] >= length)
	{
		if (CD.MinTour[CITYNO] > length && MSGMODE == 1)
		{
			cout << BIMSIZE << "BIM Tour:" << length << endl;
		}
		//		out << "Tour:" << length << endl;
		for (int n = 0; n < CITYNO; n++)
		{
			CD.MinTour[n] = CD.Tour[n];
		}
		CD.MinTour[CITYNO] = (int)length;
	}
	else
	{
		if (BlockMinTourReset == 1)
		{
			for (int n = 0; n < CITYNO; n++)
			{
				CD.Tour[n] = CD.MinTour[n];
			}
		}
	}
}

void ChangeChaosNeuralNet(int t)
{
	CNG.zeta[t + 1] = KR * CNG.zeta[t] - CNG_ALPHA * CNG.x[t] + R;
	if (typeexchange)
	{
		CNG.xi[t + 1] = BETA * (CNG.AfterTour - CNG.BeforeTour) + CNG_A;
		twoopttimes++;
	}
	else
	{
		CNG.xi[t + 1] = BETA * (CNG.AfterTour - CNG.BeforeTour) + CNG_A;
		blocktimes++;
	}
	CNG.x[t + 1] = sigmoid(CNG.xi[t + 1] + CNG.zeta[t + 1]);
}