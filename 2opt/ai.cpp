#include <iostream>
#include <vector>
using namespace std;

vector<int> da;
int ans = 0, all = 0;

int city(int cun, vector<int> map)
{
    for (int i = 0; i < da.size(); i++)
    {
        if (map[i] == 0)
        {
            map[i] = cun;
            city(cun + 1, map);
            map[i] = 0;
        }
    }
    if (cun == da.size() + 1)
    {
        all++;
        cun = 0;
        for (int i = 1; i < da.size(); i++)
        {
            if (abs(map[i - 1] - map[i]) == 1)
                cun++;
        }
        if (abs(map[0] - map[da.size() - 1]) == 1)
            cun++;
        if (cun >= da.size() - 2)
        {
            for (int i = 0; i < da.size(); i++)
            {
                cout << map[i] << " ";
            }
            cout << endl;
            ans++;
        }
    }
    return 0;
}
int main(void)
{
    int a, tmp;
    vector<int> map;
    cout << "都市数：";
    cin >> a;
    for (int i = 0; i < a; i++)
    {
        da.push_back(i + 1);
        map.push_back(0);
    }
    city(1, map);
    cout << ans << ":" << all;
}