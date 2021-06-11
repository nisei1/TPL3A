#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <sstream>

#define inName "br17a.txt"
#define cityNum 17
#define twoOptTimes 10000

int main(void)
{
	time_t t;//乱数生成のために現在時刻を格納
	int r1, r2;//ランダムに選んだ地点を格納
    int d1, d2, x;//txtデータを一時的に格納
    int bestDist = 0;//最適な経路の合計コストを格納
    int candDist = 0;//新たな経路のコストを格納
	std::vector<int> bestCity;//最適な経路を格納
	std::vector<int> candCity;//新たな経路を格納
	std::vector<std::vector<int>> edge(cityNum, std::vector<int>(cityNum));//txtデータのコストを格納

	std::ifstream ifs(inName);

	if (ifs.fail())
	{
		std::cerr << "Cannot open file\n";
		exit(0);
	}
	std::string str;
	
	while (getline(ifs, str))
	{
		std::stringstream ss(str);
		ss >> d1 >> d2 >> x;
		edge[d1 - 1][d2 - 1] = x;
	}

	/* for (int i = 0; i < cityNum; i++)
	{
		for (int j = 0; j < cityNum; j++)
		{
			std::cout << "data[" << i << "]"
				      << "[" << j << "]"
				      << "="
				      << " "
				      << edge.at(i).at(j) << std::endl;
		}
	}
	std::cout << std::endl; */

	for (int i = 0; i < cityNum; i++)
	{
		bestCity.push_back(i);
	}

	t = time(NULL);
	std::random_device seed_gen;
	std::mt19937 engine(seed_gen() * t);
	std::shuffle(bestCity.begin(), bestCity.end(), engine);

	/* for(int i = 0; i < cityNum; i++)
	{
		candCity.push_back(bestCity[i]);
	} */

	std::cout << "<before> ";
	for (int i = 0; i < cityNum; i++)
	{
		std::cout << bestCity.at(i) << " ";
	}
	std::cout << std::endl;

	std::cout << "<Total Distance> ";
	for (int i = 0; i < cityNum - 1; i++)
	{
		bestDist += edge[bestCity[i]][bestCity[i + 1]];
	}
    bestDist += edge[bestCity[cityNum - 1]][bestCity[0]];
	std::cout << bestDist << std::endl;
	std::cout << "\n";

	t = time(NULL);
	std::random_device rnd;
	std::mt19937 mt(rnd() * t);

	std::uniform_int_distribution<> dist1(0, cityNum - 1);
	std::uniform_int_distribution<> dist2(0, cityNum - 1);

	for (int i = 0; i < twoOptTimes; i++)
	{
		r1 = dist1(mt);
		r2 = dist2(mt);

		/* std::cout << "TIME:" << i + 1 << "\t"
			      << "point:" << r1 << " " << r2 << std::endl; */

		while (r1 == r2 || (r1 == 0 && r2 == cityNum - 1) || (r2 == 0 && r1 == cityNum - 1) || std::abs(r1 - r2) == 1)
		{
			r1 = dist1(mt);
			r2 = dist2(mt);
			/* std::cout << "TIME:" << i + 1 << "\t"
					  << "point:" << r1 << " " << r2 << std::endl; */
		}

		if (r1 > r2)
		{
			int tem = r1;
			r1 = r2;
			r2 = tem;
		}

		int range1, range2;
		range1 = edge[r1][r1 + 1] + edge[r2][r2 + 1];
		range2 = edge[r1][r2] + edge[r1 + 1][r2 + 1];

		if (range1 > range2)
		{
            candCity.clear();
            //最適な経路を新たな経路に代入
            for(int i = 0; i < cityNum; i++)
	        {
		        candCity.push_back(bestCity[i]);
	        }
            //新たな経路で入れ替えを行う
			while (r1 < r2)
			{
				std::swap(candCity[r1++], candCity[r2--]);
			}
            //コピーしたベクターの移動距離の和を出す
            candDist = 0;
			for (int i = 0; i < cityNum - 1; i++)
			{
				candDist += edge[candCity[i]][candCity[i + 1]];
			}
            candDist += edge[candCity[cityNum - 1]][candCity[0]];
			/* コピーしたベクターの移動距離の和が、コピー元のベクターの和より小さければ
            ベクターの中身と総和を入れ替える */
            if(candDist < bestDist)
            {
                bestCity.clear();
                bestCity = candCity;
                bestDist = candDist;
            }

			/* std::cout << "\n" << "<after> ";
			for (int i = 0; i < cityNum; i++)
			{
				std::cout << bestCity.at(i) << " ";
			}
            std::cout << std::endl;
			std::cout << "<Total Distance> ";
            std::cout << bestDist << std::endl;
            std::cout << "\n"; */
        }
    }

    std::cout << "--------------------" << std::endl;
	std::cout << "<bestAnswer> ";
	for (int i = 0; i < cityNum; i++)
	{
		std::cout << bestCity.at(i) << " ";
	}
    std::cout << std::endl;

    std::cout << "<Total Distance> "
	          << bestDist << std::endl
	          << "\n";

	return 0;
}
