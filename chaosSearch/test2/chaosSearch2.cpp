#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>
#include <random>
#include <algorithm>
#include <cmath>

/***定数の宣言***/
//ファイル入出力用の定数宣言
#define IN_NAME "../../TSP/rbg443a.txt" //読み込みたいTSP
#define OUT_NAME "out.txt"				//書き出したいファイル名
//TSP用定数宣言
#define CITY_NUM 443 //TSPの都市数
//ランダムに2-optする時用の定数宣言
#define ENABLE_TWO_OPT_RANDOM false //ランダム2-opt 有効= true ,無効 = false
#define TWO_OPT_TIMES 10			//2optで何回最小値を出すか,最小値を出すまでループで減らない
//カオスサーチで使う定数の宣言
#define T_TIMES 10				 //時刻tがどこまで増やすのか適当に決めてよい
#define ENABLE_CHAOS_SEARCH true //カオスサーチするか 有効= true ,無効 = false
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
	ChaosNN() : zai(CITY_NUM, 0.0),
				eta(CITY_NUM, 0.0),
				zeta(CITY_NUM, 0.0),
				x(CITY_NUM, 0.0),
				delta_i(CITY_NUM, 0),
				delta_j(CITY_NUM, 0) {}
	std::vector<double> zai; //ξ = zai とする。xiと書くのが一般的だがニューロンxの都市iの出力と混同しないようにするため。
	std::vector<double> eta;
	std::vector<double> zeta;
	std::vector<double> x;
	std::vector<int> delta_i; //時刻tニューロンiのときのiを格納する(TODO:必要ない気がするけど良い書き方が思い浮かばない)
	std::vector<int> delta_j; //時刻tニューロンiときのjを格納する
							  // ChaosNN() : isMaxX(CITY_NUM, false) {}
							  // vector<bool> isMaxX; //時刻tのニューロンx[i]が最大値かどうか
};

/***グローバル変数の宣言***/
//2-optの時にも利用
std::vector<std::vector<int>> edge(CITY_NUM, std::vector<int>(CITY_NUM)); //TSPを表す2次元vector、例) 1 -> 2 のコストを要素へ記録 edge[0][1] = 1から2へ行くためのコスト格納
std::vector<int> city;													  //巡回路用vector、edgeの要素番号へ入れるのに使用
std::ofstream out(OUT_NAME);											  //ファイル出力用変数
std::ifstream in(IN_NAME);												  //ファイル入力用変数
//ランダム用儀式
time_t t = time(NULL);
std::random_device seed_gen;
std::mt19937 engine(seed_gen() * t);
//カオスサーチに利用
std::vector<ChaosNN> cnn(T_TIMES); //カオスサーチ用のvector<ChaosNN>
// MaxDelta_I_J maxIJ;			  //Δijが最大にコストを下げた時のijを格納

/***テンプレート関数の宣言***/
void inputTSP(void);						  //TSPをedgeへ格納する関数
void makeFirstTour(void);					  //初回巡回路をランダムに作成する関数
inline int calcDistance(void);				  //巡回路の総コスト計算関数(戻り値:(int)巡回路の総コスト)
inline void twoOptRandom(void);				  //ランダムな2点を選んで2-opt交換する関数
inline bool twoOptPermission(int p1, int p2); //ランダム2-opt可能な2点かどうか判定(引数:都市1,都市2)(戻り値:true or false)
inline void twoOptSwap(int p1, int p2);		  //ランダム2-opt交換実行関数(引数:都市1,都市2)
inline bool swapPermission(int p1, int p2);	  //2-opt可能な2点かどうか判定(引数:都市1,都市2)(戻り値:true or false)
void initializeChaosNN(void);				  //時刻tの時の初期値
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
		out << "EnableTwoOptRandom" << std::endl;
		makeFirstTour();
		//最適化前に巡回路出力
		out << "<before>\t";
		for (auto i : city)
		{
			out << city.at(i) + 1 << " ";
		}
		out << std::endl;

		twoOptRandom();

		//最適化後に巡回路出力
		out << "<After>\t";
		for (auto i : city)
		{
			out << city.at(i) + 1 << " ";
		}
		out << std::endl;
	}

	if (ENABLE_CHAOS_SEARCH)
	{
		out << "EnableChaosSearch" << std::endl;
		makeFirstTour();

		//最適化前に巡回路出力
		out << "<before>\t";
		for (auto i : city)
		{
			out << city.at(i) + 1 << " ";
		}
		out << std::endl;

		initializeChaosNN();

		for (int t = 1; t < T_TIMES; t++) //tが0回目のときはinitializeChaosNN()で初期化した値とする。tが1回目から0回目の情報を使ってカオスニューラルネットワークの状態を更新していく
		{
			for (int i = 0; i < CITY_NUM; i++)
			{
				cnn[t].x[i] = calcX(t, i);

				if (cnn[t].x[i] >= 0.5)
				{
					if (swapPermission(cnn[t].delta_i[i], cnn[t].delta_j[i])) //TODO:エラーをスルーするようなif文。本来いらないはず-> 最大値のij
					{
						twoOptSwap(cnn[t].delta_i[i], cnn[t].delta_j[i]);
						int distance = calcDistance();
						out << "debug:After Chaos Search Total Distance:\t" << distance << std::endl;
					}
					else
					{
						out << "debug:ERROR main swapPermission false" << std::endl
							<< "t = " << t
							<< ", i = " << i << std::endl
							<< "cnn[t].delta_i[i] = " << cnn[t].delta_i[i] << ", cnn[t].delta_j[i] = " << cnn[t].delta_j[i] << std::endl
							<< "city[cnn[t].delta_i[i]] = " << city[cnn[t].delta_i[i]] << ", city[cnn[t].delta_j[i]] = " << city[cnn[t].delta_j[i]]
							<< std::endl;
					}
				}
				else
				{
					continue;
				}
			}
			out << "debug:Remaining t times:\t" << T_TIMES - t - 1 << std::endl;
		}

		//最適化後に巡回路出力
		out << "<After>\t";
		for (auto i : city)
		{
			out << city.at(i) + 1 << " ";
		}
		out << std::endl;
		int distance = calcDistance();
		out << "debug:After Chaos Search Total Distance:\t" << distance << std::endl;
	}
	return 0;
}

//TSPをedgeへ格納する関数
void inputTSP(void)
{
	int d1, d2, x;
	if (in.fail())
	{
		std::cerr << "inputTSP:Cannot open file\n";
		exit(0);
	}
	std::string str;
	while (getline(in, str)) //txtの中身一行ずつループ
	{
		std::stringstream ss(str);
		ss >> d1 >> d2 >> x;	  //1行に空白が入る度、入れる変数が変わる
		edge[d1 - 1][d2 - 1] = x; //TSPを表す2次元vectorへ格納、例) 1 -> 2 のコストを要素へ記録 edge[0][1] = 1から2へ行くためのコスト格納
	}

	//入力されたデータを一旦出力
	// for (int i = 0; i < CITY_NUM; i++)
	// {
	// 	for (int j = 0; j < CITY_NUM; j++)
	// 	{
	// 		out << "data[" << i << "]"
	// 			<< "[" << j << "]"
	// 			<< "="
	// 			<< " "
	// 			<< edge.at(i).at(j) << std::endl;
	// 	}
	// }
	// out << std::endl;
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
	for (int i = 0; i < CITY_NUM - 1; i++)
	{
		sum += edge[city[i]][city[i + 1]];
	}
	sum += edge[city[CITY_NUM - 1]][city[0]];

	return sum;
}

//ランダムな2点を選んで2-opt交換する関数
inline void twoOptRandom(void)
{
	int times = TWO_OPT_TIMES; //2optで何回最小値を出すか,最小値を出すまでループで減らない
	//2opt用2点ランダム生成範囲
	std::uniform_int_distribution<> dist1(0, CITY_NUM - 1);
	std::uniform_int_distribution<> dist2(0, CITY_NUM - 1);

	do
	{
		int r1, r2; //2点保存用変数
		//2opt point set
		do //2optができる2点をランダムに生成
		{
			r1 = dist1(engine); //ランダム生成点1
			r2 = dist2(engine); //ランダム生成点2
		} while (!(twoOptPermission(r1, r2)));

		std::vector<int> oldCity = city; //最短ルート保存用vector、2optの前後で合計のコストと比較し2opt後でコストが増えればこの変数を利用し、ロールバックする

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
			out << "debug:After 2opt Total Distance:\t" << distance2 << std::endl
				<< "debug:Remaining 2opt times:\t" << TWO_OPT_TIMES - times + 1 << std::endl;

			times--;
		}

	} while (times > 0);
}

inline bool twoOptPermission(int p1, int p2)
{
	int range1, //繋ぎ変える前の2辺のコスト
		range2; //繋ぎ変えた後の2辺のコスト

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
			if (range1 >= range2)
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
			if (range1 >= range2)
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

inline bool swapPermission(int p1, int p2)
{
	if (abs(p2 - p1) <= 1)
	{
		return false;
	}
	else if (p1 == 0 && p2 == CITY_NUM - 1)
	{
		return false;
	}
	else if (p1 == CITY_NUM - 1 && p2 == 0)
	{
		return false;
	}
	else
	{
		return true;
	}
}

inline void twoOptSwap(int p1, int p2)
{
	if (swapPermission(p1, p2))
	{
		if (p1 > p2)
		{
			int temp = p1;
			p1 = p2;
			p2 = temp;
		}

		int swapP1 = p1 + 1,
			swapP2 = p2;

		while (swapP1 < swapP2)
		{
			std::swap(city[swapP1++], city[swapP2--]);
		}
	}
	else
	{
		std::cout << "ERROR:2-optSwap, Not swapPermission Point" << std::endl;
		exit(0);
	}
}

inline double sigmoid(double x)
{
	return 1.0 / (1.0 + std::exp(-x / EPSILON));
}

inline double calcZai(int t, int i)
{
	double max = 0.0;
	bool isFirst = true; //初期max代入時に利用
	int cityIndex_i = 0; //都市番号とニューロン番号を一致させるため、for文でcity[]の中身とiが一致->都市番号iとニューロン番号iが一致した時のcity[]の要素番号を格納(もっと効率良い方法ありそう)
	int cityIndex_j = 0;
	for (int k = 0; k < CITY_NUM; k++)
	{
		if (city[k] == i)
		{
			cityIndex_i = k;
			break;
		}
	}

	cnn[t].delta_i[i] = cityIndex_i;

	for (int j = 0; j < CITY_NUM; j++)
	{
		for (int k = 0; k < CITY_NUM; k++)
		{
			if (city[k] == j)
			{
				cityIndex_j = k;
				break;
			}
		}

		if (!(swapPermission(cityIndex_i, cityIndex_j)))
		{
			continue;
		}
		else
		{
			double sumZetaBetaDelta = calcZeta(t, j) + BETA * calcDelta(cityIndex_i, cityIndex_j);
			if (isFirst)
			{
				max = sumZetaBetaDelta;
				cnn[t].delta_j[i] = cityIndex_j;
				isFirst = false;
				// out << "debug: calcZai First i = " << i
				// 	<< " j = " << j << std::endl
				// 	<< "city[cityIndex_i] = " << city[cityIndex_i]
				// 	<< " city[cityIndex_j] = " << city[cityIndex_j] << std::endl
				// 	<< "cityIndex_i = " << cityIndex_i
				// 	<< " cityIndex_j = " << cityIndex_j
				// 	<< std::endl;
			}
			else if (max <= sumZetaBetaDelta)
			{
				max = sumZetaBetaDelta;
				cnn[t].delta_j[i] = cityIndex_j;
				// out << "debug: calcZai i = " << i
				// 	<< " j = " << j << std::endl
				// 	<< "city[cityIndex_i] = " << city[cityIndex_i]
				// 	<< " city[cityIndex_j] = " << city[cityIndex_j] << std::endl
				// 	<< "cityIndex_i = " << cityIndex_i
				// 	<< " cityIndex_j = " << cityIndex_j
				// 	<< std::endl;
			}
		}
	}

	return max;
}

inline double calcEta(int t, int i)
{
	t = t - 1; //tはt+1から入力されるため、前の時刻tを利用するためにtから1を引いている
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

inline double calcZeta(int t, int i) //Zetaはdを含まない形で実装(7)
{
	t = t - 1; //tはt+1から入力されるため、前の時刻tを利用するためにtから1を引いている
	return (KR * cnn[t].zeta[i]) - (ALPHA * cnn[t].x[i]) + R;
}

inline int calcDelta(int i, int j)
{
	if (swapPermission(i, j))
	{
		std::vector<int> oldCity = city;  //最短ルート保存用vector、2optの前後で合計のコストと比較し2opt後でコストが増えればこの変数を利用し、ロールバックする
		int oldDistance = 0,			  //i-j間の2-opt前巡回路コスト総計
			newDistance = 0;			  //i-j間の2-opt後巡回路コスト総計
		oldDistance = calcDistance();	  //swap前にコストを入れる
		twoOptSwap(i, j);				  //swap
		newDistance = calcDistance();	  //swap後のコストを入れる
		city = oldCity;					  //Δijの計算のたびにロールバック
		return oldDistance - newDistance; //戻り値は削減されればプラスの値を取る
	}
	else
	{
		std::cout << "ERROR:calcDelta swapPermission is false" << std::endl;
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

void initializeChaosNN(void)
{
	//xの初期値もその後の計算後も0か1にしかならない -> それでいい、そうなってしまうもの
	std::uniform_real_distribution<> distr(0.0, 1.0);
	for (int t = 0; t < T_TIMES; t++)
	{
		for (int i = 0; i < CITY_NUM; i++)
		{
			cnn[t].zai[i] = distr(engine);
			cnn[t].eta[i] = distr(engine);
			cnn[t].zeta[i] = distr(engine);
			cnn[t].x[i] = sigmoid(cnn[t].zai[i] + cnn[t].eta[i] + cnn[t].zeta[i]);
		}
	}
}
