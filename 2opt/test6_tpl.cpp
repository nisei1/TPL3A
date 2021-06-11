#include <iostream>
#include <vector>
#include <random>
// #include <ctime>
#include <algorithm>
#include <fstream>
#include <sstream>

using namespace std;

#define inName "br17a.txt"
#define cityNum 17
#define twoOptTimes 50 //2optで何回最小値を出すか,最小値を出すまでループで減らない

int main(void)
{
	// time_t t;
	int times = twoOptTimes; //2optで何回最小値を出すか,最小値を出すまでループで減らない

	vector<vector<int>> edge(cityNum, vector<int>(cityNum)); //inputdata用2次元配列、例) 1 -> 2 のコストを要素へ記録 edge[0][1] = cost;
	int d1, d2, x;											 //edgeへinputdataを入れるための一時的な変数
	vector<int> city;										 //巡回路用vector、edgeの要素番号へ入れるのに使用
	// do
	// {
	// 	cin >> d1 >> d2 >> x;
	// 	edge[d1 - 1][d2 - 1] = x;
	// } while (!(d1 == 17 && d2 == 17));

	//ファイル入力
	ifstream ifs(inName);

	if (ifs.fail())
	{
		cerr << "Cannot open file\n";
		exit(0);
	}
	string str;
	while (getline(ifs, str)) //txtの中身一行ずつループ
	{
		stringstream ss(str);
		ss >> d1 >> d2 >> x;	  //1行に空白が入る度、入れる変数を変える
		edge[d1 - 1][d2 - 1] = x; //inputdata用2次元配列へ、例) 1 -> 2 のコストを要素へ記録 edge[0][1] = cost;
	}
	//入力されたデータを一旦出力
	for (int i = 0; i < cityNum; i++)
	{
		for (int j = 0; j < cityNum; j++)
		{
			cout << "data[" << i << "]"
				 << "[" << j << "]"
				 << "="
				 << " "
				 << edge.at(i).at(j) << endl;
		}
	}
	cout << endl;

	//巡回路vectorに0~citynum-1までの数を順番に入れる
	for (int i = 0; i < cityNum; i++)
	{
		city.push_back(i);
	}

	//ランダム用儀式
	std::random_device seed_gen;
	std::mt19937 engine(seed_gen());
	std::shuffle(city.begin(), city.end(), engine); //巡回路シャッフル

	//最適化前に出力
	cout << "<before>\t";
	for (int i = 0; i < cityNum; i++)
	{
		cout << city.at(i) << " ";
	}
	cout << endl;
	int beforeSum = 0;
	cout << "<Total Distance>\t";
	for (int i = 0; i < cityNum - 1; i++)
	{
		beforeSum += edge[city[i]][city[i + 1]];
	}
	cout << beforeSum << endl;
	cout << "\n";

	//2opt用2点ランダム生成範囲
	std::uniform_int_distribution<> dist1(0, cityNum - 1);
	std::uniform_int_distribution<> dist2(0, cityNum - 1);

	time_t i = 0; //何回r1,r2を振ったかを表示,以下のdoWhile用のループ変数
	do
	{
		int r1, r2; //2点保存用変数
		//2opt point set
		do //2optができる2点をランダムに生成
		{
			r1 = dist1(engine);						   //ランダム生成点1
			r2 = dist2(engine);						   //ランダム生成点2
			cout << "Set point TIME:" << i + 1 << "\t" //何回ランダム生成したかを出力
				 << "point:"
				 << "\t" << r1 << "\t" << r2 << endl; //出力された点を表示
			i++;
		} while (r1 == r2 || (r1 == 0 && r2 == cityNum - 1) || (r2 == 0 && r1 == cityNum - 1) || std::abs(r1 - r2) == 1); //同じ点でなく、始点と終点でなく、2点の差の絶対値が1より大きくなるまでランダム生成

		//常にr2が大きいようにする
		if (r1 > r2)
		{
			int tem = r1;
			r1 = r2;
			r2 = tem;
		}

		//切られた2組の辺のコスト(2次元配列edgeの要素)を組み替える前後と比較して2optすべきか判定
		int range1, range2;
		range1 = edge[r1][r1 + 1] + edge[r2][r2 + 1]; //繋ぎ変える前
		range2 = edge[r1][r2] + edge[r1 + 1][r2 + 1]; //繋ぎ変えた後
		if (range1 > range2)						  //繋ぎ変えた方がよかったら(繋ぎ変える前が大きいなら)2optする
		{
			vector<int> oldCity = city; //最短ルート保存用vector、2optの前後で合計のコストと比較し2opt後でコストが増えればこの変数を利用し、ロールバックする
			int sum = 0;				//合計一時保存用変数
			// cout << "<before>\t";
			// for (int i = 0; i < cityNum; i++)
			// {
			// 	cout << city.at(i) << " ";
			// }
			// cout << endl;

			// cout << "<Before Total Distance>\t";
			//2opt前合計距離算出
			for (int i = 0; i < cityNum - 1; i++)
			{
				sum += edge[city[i]][city[i + 1]];
			}
			int distance1 = sum; //合計距離比較用変数1
			cout << "debug:before swap distance:" << distance1 << endl;

			//繋ぎ変えてみる
			while (r1 < r2)
			{
				cout << "debug:while 1" << endl;
				std::swap(city[r1++], city[r2--]);
			}
			// while (r1 > r2)
			// {
			// 	cout << "debug:while 2" << endl;
			// 	std::swap(city[r1--], city[r2++]);
			// }

			// cout << "\n"
			// 	 << "<after> ";
			// for (int i = 0; i < cityNum; i++)
			// {
			// 	cout << city.at(i) << " ";
			// }
			// cout << endl;

			sum = 0;
			// cout << "<Total Distance> ";
			//2opt後合計値算出
			for (int i = 0; i < cityNum - 1; i++)
			{
				sum += edge[city[i]][city[i + 1]];
			}
			int distance2 = sum; //合計距離比較用変数2
			cout << "debug:After swap distance:" << distance2 << endl;

			//2opt前後で合計距離比較。大きければロールバック
			if (distance1 < distance2)
			{
				city = oldCity; //ロールバック
				cout << "debug:after is larger than before" << endl;
			}
			else //小さければ現在のコストを出力し、残り2opt回数を減らす
			{
				cout << "debug:After 2opt Total Distance:\t" << distance2 << endl
					 << "debug:Remaining 2opt times:\t" << times - 1 << endl;

				times--;
			}
		}
		else //繋ぎ変えなくていいならループ
		{
			cout << "debug:range1 <= range2" << endl;
		}

	} while (times > 0);

	//最後に巡回路vector出力
	cout << "<after>\t";
	for (int i = 0; i < cityNum; i++)
	{
		cout << city.at(i) << " ";
	}
	cout << endl;
	return 0;
}
