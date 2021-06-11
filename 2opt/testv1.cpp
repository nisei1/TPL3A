#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdio>
#include <random>

using namespace std;

struct Point
{
    double x;
    double y;
};

static uint32_t _x = 2463534242;

static void xorshift32_seed(uint32_t seed)
{
    if (seed != 0)
    {
        _x = seed;
    }
}

static __inline uint32_t xorshift32(void)
{
    _x ^= (_x << 13);
    _x ^= (_x >> 17);
    return _x ^= (_x << 5);
}

static __inline uint32_t _rand()
{
    return xorshift32();
}

// connect a -> b and c -> d
__inline void swap_edges(vector<int> &seq, int a, int b, int c, int d)
{
    int x = d, y = b, N = (int)seq.size();
    if (y - x > a + N - c)
    {
        x = c, y = a + N;
    }
    while (x < y)
    {
        swap(seq[x < N ? x : x - N], seq[y < N ? y : y - N]);
        ++x, --y;
    }
}

void two_opt(const vector<Point> &points, vector<int> &seq, int turns)
{
    int N = (int)seq.size();
    vector<vector<double>> dist(N, vector<double>(N));
    for (int i = 0; i < N; ++i)
    {
        for (int j = i + 1; j < N; ++j)
        {
            double dx = points[i].x - points[j].x, dy = points[i].y - points[j].y;
            dist[i][j] = dist[j][i] = sqrt(dx * dx + dy * dy);
        }
    }

    for (int t = 0; t < turns; ++t)
    {
        int a = _rand() % N, b = _rand() % N;
        if (a == b)
            continue;
        if (a > b)
            swap(a, b);
        int d = (a + 1) % N, c = (b + 1) % N;
        if (dist[seq[a]][seq[b]] + dist[seq[c]][seq[d]] + 1e-9 < dist[seq[a]][seq[d]] + dist[seq[b]][seq[c]])
        {
            swap_edges(seq, a, b, c, d);
        }
    }
}

int main()
{
    // vector<int> seq;
    // for (size_t i = 1; i <= 5; i++)
    // {
    //     seq.push_back(i);
    // }

    // int x = d, y = b;
    // while (x < y)
    // {
    //     swap(seq[x++], seq[y--]);
    // }
    return 0;
}
