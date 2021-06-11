#include <iostream>
#include <vector>
#include <random>
using namespace std;

int main(int argc, char const *argv[])
{
    // int seq[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    int n = 5;     //都市数
    vector<int> v; //都市を格納する配列

    for (size_t i = 1; i <= n; i++) //都市の始点は1,終点はn
    {
        v.push_back(i);
    }

    // std::random_device rd;
    // std::default_random_engine eng(rd());
    // std::uniform_int_distribution<int> distr(0, n);
    // int p1 = distr(eng);
    // int p2 = distr(eng);

    // if (abs(p1 - p2) > 1)
    // {
    //     if ((p1 == n && p2 == 1) || (p2 == n && p1 == 1))
    //     {
    //         while (p1 < p2)
    //         {
    //             swap(v[p1++], v[p2--]);
    //         }
    //         // while (p2 < p1)
    //         // {
    //         //     swap(v[p2++], v[p1--]);
    //         // }
    //     }
    //     else
    //     {
    //         return 0;
    //     }
    // }
    // else
    // {
    //     return 0;
    // }

    // while (p1 < p2)
    // {
    //     int t = p1;
    //     p1 = p2;
    //     p2 = t;
    // }

    // int x = d, y = b;
    // while (x < y)
    // {
    //     swap(seq[x++], seq[y--]);
    // }
    // int x = seq[3], y = seq[7];
    // while (x < y)
    // {
    //     swap(seq[x++], seq[y--]);
    // }


    //出力
    for (size_t i = 0; i < n; i++)
    {
        cout << v[i] << "\t";
    }

    return 0;
}
