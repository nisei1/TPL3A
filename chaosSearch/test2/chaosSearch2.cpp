#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>
#include <random>
#include <algorithm>

using namespace std;

/***定数の宣言***/
//ファイル入出力用の定数宣言
#define IN_NAME "../../TSP/rbg443a.txt" //読み込みたいTSP
#define OUT_NAME "out.txt"				//書き出したいファイル名
//TSP用定数宣言
#define CITY_NUM 443 //TSPの都市数
//ランダムに2-optする時用の定数宣言
#define ENABLE_TWO_OPT_RANDOM true //ランダム2-opt 有効= true ,無効 = false
#define TWO_OPT_TIMES 10		   //2optで何回最小値を出すか,最小値を出すまでループで減らない
//カオスサーチで使う定数の宣言
#define T_TIMES 100				  //TODO:時刻tがどこまでふえるのかわからないので宣言してみた定数
#define ENABLE_CHAOS_SEARCH false //カオスサーチするか 有効= true ,無効 = false
#define ALPHA 1.0
#define BETA 75.0
#define THETA 0.5
#define KR 0.99
#define R THETA *(1.0 - KR)
#define WEIGHT 0.15
#define EPSILON 0.002

/***クラス・構造体の宣言***/
//カオスサーチ用の構造体
struct ChaosNN
{
	ChaosNN() : zai(CITY_NUM, 0.0) {}
	vector<double> zai;
	ChaosNN() : zeta(CITY_NUM, 0.0) {}
	vector<double> zeta;
	ChaosNN() : eta(CITY_NUM, 0.0) {}
	vector<double> eta;
	ChaosNN() : x(CITY_NUM, 0.0) {}
	vector<double> x;
};
//Δijが最大にコストを下げた時のijを格納するための構造体
struct MaxDelta_I_J
{
	int i;
	int j;
};

/***グローバル変数の宣言***/
//2-optの時にも利用
vector<vector<int>> edge(CITY_NUM, vector<int>(CITY_NUM)); //TSPを表す2次元vector、例) 1 -> 2 のコストを要素へ記録 edge[0][1] = 1から2へ行くためのコスト格納
vector<int> city;										   //巡回路用vector、edgeの要素番号へ入れるのに使用
ofstream out(OUT_NAME);									   //ファイル出力用変数
ifstream in(IN_NAME);									   //ファイル入力用変数
//ランダム用儀式
time_t t = time(NULL);
random_device seed_gen;
mt19937 engine(seed_gen() * t);
//カオスサーチに利用
vector<ChaosNN> cnn(T_TIMES); //カオスサーチ用のvector<ChaosNN>
MaxDelta_I_J maxIJ;			  //Δijが最大にコストを下げた時のijを格納

/***テンプレート関数の宣言***/
void inputTSP(void);						  //TSPをedgeへ格納する関数
void makeFirstTour(void);					  //初回巡回路をランダムに作成する関数
inline int calcDistance(void);				  //巡回路の総コスト計算関数(戻り値:(int)巡回路の総コスト)
inline void twoOptRandom(void);				  //ランダムな2点を選んで2-opt交換する関数
inline bool twoOptPermission(int p1, int p2); //2-opt可能な2点かどうか判定(引数:都市1,都市2)(戻り値:true or false)
inline void twoOptSwap(int p1, int p2);		  //2-opt交換実行関数(引数:都市1,都市2)
void initialize(void);						  //TODO:時刻tの時の初期値
inline double sigmoid(double x);			  //シグモイド関数
inline double calcZai(int t, int i);		  //(3)式関数
inline double calcEta(int t, int i);		  //(4)式関数
inline int calcDelta(int i, int j);			  //(3)式のΔij関数
inline double calcZeta(int t, int i);		  //(5)式関数
inline double calcX(int t, int i);			  //(6)式関数

/***main関数***/
int main(int argc, char const *argv[])
{
	inputTSP();

	if (ENABLE_TWO_OPT_RANDOM)
	{
		out << "EnableTwoOptRandom" << endl;
		makeFirstTour();
		//最適化前に巡回路出力
		out << "<before>\t";
		for (auto i : city)
		{
			out << city.at(i) + 1 << " ";
		}
		out << endl;

		twoOptRandom();

		//最適化後に巡回路出力
		out << "<After>\t";
		for (auto i : city)
		{
			out << city.at(i) + 1 << " ";
		}
		out << endl;
	}

	if (ENABLE_CHAOS_SEARCH)
	{
		out << "EnableChaosSearch" << endl;
		makeFirstTour();

		//最適化前に巡回路出力
		out << "<before>\t";
		for (auto i : city)
		{
			out << city.at(i) + 1 << " ";
		}
		out << endl;

		for (int t = 0; t < CITY_NUM; t++)
		{
			for (int i = 0; i < T_TIMES; i++)
			{
				cnn[t].x[i] = calcX(t, i);
			}
		}
	}

	return 0;
}

//TSPをedgeへ格納する関数
void inputTSP(void)
{
	int d1, d2, x;
	if (in.fail())
	{
		cerr << "Cannot open file\n";
		exit(0);
	}
	string str;
	while (getline(in, str)) //txtの中身一行ずつループ
	{
		stringstream ss(str);
		ss >> d1 >> d2 >> x;	  //1行に空白が入る度、入れる変数が変わる
		edge[d1 - 1][d2 - 1] = x; //TSPを表す2次元vectorへ格納、例) 1 -> 2 のコストを要素へ記録 edge[0][1] = 1から2へ行くためのコスト格納
	}

	//入力されたデータを一旦出力
	for (int i = 0; i < CITY_NUM; i++)
	{
		for (int j = 0; j < CITY_NUM; j++)
		{
			out << "data[" << i << "]"
				<< "[" << j << "]"
				<< "="
				<< " "
				<< edge.at(i).at(j) << endl;
		}
	}
	out << endl;
}

//初回巡回路をランダムに作成する関数
void makeFirstTour(void)
{
	//巡回路vectorの中身を空にする
	if (city.empty() == false)
	{
		city.clear();
	}

	//巡回路vectorに0~citynum-1までの数を順番に入れる
	for (int i = 0; i < CITY_NUM; i++)
	{
		city.push_back(i);
	}

	//巡回路シャッフル
	shuffle(city.begin(), city.end(), engine);
}

//巡回路の総コスト計算関数
inline int calcDistance(void)
{
	int sum = 0;
	// out << "<Total Distance>\t";
	for (int i = 0; i < CITY_NUM - 1; i++)
	{
		sum += edge[city[i]][city[i + 1]];
	}
	sum += edge[city[CITY_NUM - 1]][city[0]];
	// out << sum << endl
	// 	<< endl;

	return sum;
}

//ランダムな2点を選んで2-opt交換する関数
inline void twoOptRandom(void)
{
	int times = TWO_OPT_TIMES; //2optで何回最小値を出すか,最小値を出すまでループで減らない
	//2opt用2点ランダム生成範囲
	uniform_int_distribution<> dist1(0, CITY_NUM - 1);
	uniform_int_distribution<> dist2(0, CITY_NUM - 1);

	do
	{
		int r1, r2; //2点保存用変数
		//2opt point set
		do //2optができる2点をランダムに生成
		{
			r1 = dist1(engine); //ランダム生成点1
			r2 = dist2(engine); //ランダム生成点2
		} while (!(twoOptPermission(r1, r2)));

		vector<int> oldCity = city; //最短ルート保存用vector、2optの前後で合計のコストと比較し2opt後でコストが増えればこの変数を利用し、ロールバックする

		int distance1 = calcDistance(); //合計距離比較用変数1

		//繋ぎ変えてみる
		twoOptSwap(r1, r2);

		int distance2 = calcDistance(); //合計距離比較用変数2

		//2opt前後で合計距離比較。大きければロールバック
		if (distance1 < distance2)
		{
			city = oldCity; //ロールバック
		}
		else //小さいまたは同じコストの場合、現在のコストを出力し、残り2opt回数を減らす
		{
			out << "debug:After 2opt Total Distance:\t" << distance2 << endl
				<< "debug:Remaining 2opt times:\t" << TWO_OPT_TIMES - times + 1 << endl;

			times--;
		}

	} while (times > 0);
}

inline bool twoOptPermission(int p1, int p2)
{
	int range1, //繋ぎ変える前の2辺のコスト
		range2; //繋ぎ変えた後の2辺のコスト

	//TODO:whileの条件がわかりにくいからif文で書きたい
	if (p1 > p2)
	{
		return false;
	}
	else if (abs(p2 - p1) <= 1)
	{
		return false;
	}
	else if (p1 == 0 && p2 == CITY_NUM - 1)
	{
		return false;
	}
	else
	{
		if (p2 == CITY_NUM - 1)
		{
			range1 = edge[p1][p1 + 1] + edge[p2][0]; //繋ぎ変える前
			range2 = edge[p1][p2] + edge[p1 + 1][0]; //繋ぎ変えた後
			if (range1 >= range2)					 // TODO:">"だけだと無限ループになるのか要検証
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			range1 = edge[p1][p1 + 1] + edge[p2][p2 + 1]; //繋ぎ変える前
			range2 = edge[p1][p2] + edge[p1 + 1][p2 + 1]; //繋ぎ変えた後
			if (range1 >= range2)						  // TODO:">"だけだと無限ループになるのか要検証
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}
}

inline void twoOptSwap(int p1, int p2)
{
	if (twoOptPermission(p1, p2))
	{
		int swapP1 = p1 + 1,
			swapP2 = p2;

		while (swapP1 < swapP2)
		{
			swap(city[swapP1++], city[swapP2--]);
		}
	}
	else
	{
		cout << "ERROR:2-optSwap, Not twoOptPermission Point" << endl;
		exit(0);
	}
}

inline double sigmoid(double x)
{
	return 1.0 / (1.0 + exp(-x / EPSILON));
}

inline double calcZai(int t, int i)
{
	double max = 0.0;
	bool isFirst = true;		//初期max代入時に利用
	vector<int> oldCity = city; //最短ルート保存用vector、2optの前後で合計のコストと比較し2opt後でコストが増えればこの変数を利用し、ロールバックする
	for (int j = 0; j < CITY_NUM; j++)
	{
		if (!(twoOptPermission(i, j)))
		{
		}
		else
		{
			double sumZetaBetaDelta = calcZeta(t, j) + BETA * calcDelta(i, j);
			if (isFirst)
			{
				max = sumZetaBetaDelta;
				isFirst = false;
			}
			else if (max <= sumZetaBetaDelta)
			{
				max = sumZetaBetaDelta;
				maxIJ.i = i;
				maxIJ.j = j;
			}
			else
			{
			}
		}
		city = oldCity; //Δijの計算のたびにロールバック
	}
	return max;
}

inline double calcEta(int t, int i)
{
	t = t - 1;
	double sum = 0.0;
	for (int k = 0; k < CITY_NUM; k++) //kはiが0から始まるため、0でよい
	{
		if (k != i)
		{
			sum += cnn[t].x[k];
		}
	}
	return -WEIGHT * sum + WEIGHT;
}

inline double calcZeta(int t, int i)
{
	t = t - 1;
	double sum = 0.0;

	for (int d = 0; d <= t; d++)
	{
		sum += pow(KR, d) * cnn[t - d].x[i];
	}

	return -ALPHA * sum + THETA;
}

inline int calcDelta(int i, int j)
{
	if (twoOptPermission(i, j))
	{
		//TODO:
		int oldDistance = 0,			  //i-j間の2-opt前巡回路コスト総計
			newDistance = 0;			  //i-j間の2-opt後巡回路コスト総計
		oldDistance = calcDistance();	  //swap前にコストを入れる
		twoOptSwap(i, j);				  //swap
		newDistance = calcDistance();	  //swap後のコストを入れる
		return oldDistance - newDistance; //TODO:戻り値は？
	}
	else
	{
		cout << "ERROR:calcDelta twoOptPermission is false" << endl;
		exit(0);
	}
}

inline double calcX(int t, int i)
{
	cnn[t].eta[i] = calcEta(t, i);
	cnn[t].zeta[i] = calcZeta(t, i);
	cnn[t].zai[i] = calcZai(t, i);
	return sigmoid(cnn[t].zai[i] + cnn[t].eta[i] + cnn[t].zeta[i]);
}

void initialize(void)
{
	//TODO:時刻tの時の初期値
}