#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdio>
#include <random>

using namespace std;

//構造体double x, double y
struct Point
{
    double x;
    double y;
};

static uint32_t _x = 2463534242; //グローバル変数シード値

//seedを引数に受け取ったseedが０でないとき_xにseedを代入
static void xorshift32_seed(uint32_t seed)
{
    if (seed != 0) //受け取ったseedが０でないとき
    {
        _x = seed; //_xにseedを代入
    }
}

//引数無し、戻り値は_xをビットシフトしたりXORしたり？
static __inline uint32_t xorshift32(void)
{
    _x ^= (_x << 13); //_xに13個左ビットシフト(xx = 0x12345678 << 8;  // xx には 0x34567800 が代入される)後、XOR(xx = 0xffff0000 ^ 0xff00ff00;  // xx には 0x00ffff00 が代入される)
    _x ^= (_x >> 17);
    return _x ^= (_x << 5);
}

//引数無し,戻り値はxorshift32関数を呼び出す
static __inline uint32_t _rand()
{
    return xorshift32();
}

// connect a -> b and c -> d
//配列をスワップ
__inline void swap_edges(vector<int> &seq, int a, int b, int c, int d)
{
    int x = d,
        y = b,
        N = (int)seq.size(); //vector<int> seqの要素の長さを格納
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

void two_opt(const vector<Point> &points, //都市配列引数
             vector<int> &seq,            //巡回路引数
             int turns)                   //繰り返し引数?
{
    int N = (int)seq.size();                           //巡回路の要素数を取得
    vector<vector<double>> dist(N, vector<double>(N)); //巡回路の要素数で2次元配列を宣言dist[N][N]2点間のコストを格納
    for (int i = 0; i < N; ++i)
    {
        for (int j = i + 1; j < N; ++j)
        {
            double dx = points[i].x - points[j].x,
                   dy = points[i].y - points[j].y;
            dist[i][j] = dist[j][i] = sqrt(dx * dx + dy * dy); //2点間のコストを格納 -> 対称TSPの為1->2,2->1のコストは同じ値の為、ここで代入
        }
    }

    for (int t = 0; t < turns; ++t)
    {
        int a = _rand() % N, //要素数の範囲の中でランダムな値を取得1
            b = _rand() % N; //要素数の範囲の中でランダムな値を取得2
        if (a == b)
            continue; //ループをスキップ
        if (a > b)
            swap(a, b); //a > b の時 aとbの中身を交換 -> a < b にさせるため
        int d = (a + 1) % N,    //
            c = (b + 1) % N;
        if (dist[seq[a]][seq[b]] + dist[seq[c]][seq[d]] + 1e-9 < dist[seq[a]][seq[d]] + dist[seq[b]][seq[c]])
        {
            swap_edges(seq, a, b, c, d);
        }
    }
}

int main(int argc, const char *argv[])
{
    //const int N = 10, T = 100;
    const int N = 100, T = 25000;
    //const int N = 1000, T = 2000000;

    vector<Point> points(N); //N個のpoint型vectorを宣言->都市場所(x,y)配列

    //pointsの各要素x,yにランダムな値を代入
    for (int i = 0; i < N; ++i)
    {
        //points[i].x = cos(2.0 * M_PI / N * i) * 100;
        //points[i].y = sin(2.0 * M_PI / N * i) * 100;
        points[i].x = (_rand() % 10000) * 0.01;
        points[i].y = (_rand() % 10000) * 0.01;
    }

    //ランダムにpointsの順番をシャッフル
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::shuffle(points.begin(), points.end(), engine); //シャッフル

    // random_shuffle(points.begin(), points.end());

    FILE *fp = fopen("in.txt", "wt"); //in.txtを読み書き用に作成展開
    fprintf(fp, "%d\n", N);           //fpに書き込み
    for (int i = 0; i < N; ++i)
    {
        fprintf(fp, "%.10f %.10f\n", points[i].x, points[i].y);
    }
    fclose(fp); //fpを閉じる

    vector<int> seq(N); //N個のvector<int>宣言->巡回路
    //0~N-1の並んだ要素を入れる
    for (int i = 0; i < N; ++i)
    {
        seq[i] = i;
    }

    two_opt(points, seq, T); //two_optへ。都市配列,巡回路配列,twoptをやる回数T？

    //out.txtへ巡回路を出力
    fp = fopen("out.txt", "wt");
    for (int i = 0; i < N; ++i)
    {
        fprintf(fp, " %d", seq[i]);
    }
    fclose(fp);

    return 0;
}