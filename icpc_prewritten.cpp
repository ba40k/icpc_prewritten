#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <cstring>
#include <map>

using namespace std;

namespace icpc
{
class testing_system
{
    class test
    {
    public:
        // то в каком виде хранится тест, указать руками
        test()
        {

        }
        void show_test()
        {
            // логика отображения теста
        }
    };
    class test_generator
    {
    private:
        test current;

        test generate_test() // аргументы если надо как-то ограничить генерацию
        {
            // логика рандомно создающая тест


            //
            // присвоить потом это в current
        }
    public:
        test_generator()
        {

        }
        test get_test() // аргументы если надо как-то ограничить генерацию
        {
            return generate_test(); // не забыть про аргументы если они есть
        }

    };
    class solution
    {
    public:
        solution()
        {

        }
        int get_result_from_stupid(test current) // тип данных зависящий от задачи
        {
            int result;
            // логика stupid-решения


            return result;
        }
        int get_result_from_smart(test current)
        {
            int result;
            // логика smart-решения

            return result;
        }
    };
public:
    test_generator genereator;
    solution cur_solution;
    testing_system()
    {

    }
    void test_solution(int iterations) // сколько раз надо протестировать
    {
        for (int i = 0; i < iterations; i++) {
            test random_test;
            if (cur_solution.get_result_from_smart(random_test) != cur_solution.get_result_from_stupid(random_test)) {
                cout << "MISTAKE DETECTED ! \n";
                cout << "======\n";
                random_test.show_test();
                cout << "======\n";
            }

        }
    }
};

void EuelerTourForTree(int vertex,
                       int parent,
                       int depth,
                       std::vector<std::vector<int>> &graph,
                       std::vector<std::pair<int, int>> &eueler_tour)
{
    eueler_tour.emplace_back(vertex, depth);
    int cur_vertex = vertex;
    for (int next_vertex: graph[vertex]) {
        if (next_vertex != parent) {
            EuelerTourForTree(next_vertex, vertex, depth + 1, graph, eueler_tour);
            eueler_tour.emplace_back(vertex, depth);
        }
    }
}
class SqrtDec // пустой шаблон
{
private :
    const int k = 8;
    int block_size = 1 << k; // default

    struct node
    {

    };
    std::vector<node> dec;
    std::vector<int> vec;
    void build(std::vector<int> &input)
    {
        vec = input;
        dec.resize((input.size() >> k) + 1);
        for (int i = 0; i < input.size(); ++i) {
            // creating
        }
    }
public:
    SqrtDec(std::vector<int> &input)
    {
        build(input);
    }

};
class Dsu // базовая дсу
{
private:
    std::vector<int> parent;
    std::vector<int> len;
    void build(int n)
    {
        parent.resize(n);
        len.resize(n, 1);
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }
public:
    Dsu(int n)
    {
        build(n);
    }
    int get_parent(int vertex)
    {
        if (parent[vertex] == vertex) {
            return vertex;
        }
        return parent[vertex] = get_parent(parent[vertex]);
    }
    void set_union(int a, int b)
    {
        a = get_parent(a);
        b = get_parent(b);
        if (len[a] < len[b]) {
            len[b] += len[a];
            parent[a] = b;
            return;
        }
        len[a] += len[b];
        parent[b] = a;
    }
};

namespace segtree // icpc::segtree
{
class JiDriverRemainder
{
    struct Node
    {
        int64_t maxNum;
        int64_t sum;
    };

    std::vector<Node> tree;
    int size;

    void init(int n)
    {
        size = 1;
        while (size < n) {
            size *= 2;
        }
        tree.resize(size * 2 - 1);
    }

    inline void updateFromChildren(int x)
    {
        tree[x].maxNum = std::max(tree[2 * x + 1].maxNum, tree[2 * x + 2].maxNum);
        tree[x].sum = tree[2 * x + 1].sum + tree[2 * x + 2].sum;
    }

    void build(int x, int lx, int rx, const std::vector<int64_t> &v)
    {
        if (rx - lx == 1) {
            if (lx < v.size()) {
                tree[x].maxNum = tree[x].sum = v[lx];
            }
            return;
        }
        int m = (lx + rx) / 2;
        build(2 * x + 1, lx, m, v);
        build(2 * x + 2, m, rx, v);
        updateFromChildren(x);
    }

    void updateModEqual(int64_t value, int l, int r, int x, int lx, int rx)
    {
        bool breakCondition = (lx >= r) || (rx <= l) || (tree[x].maxNum < value);
        if (breakCondition) {
            return;
        }

        if (rx - lx == 1) {
            tree[x].maxNum %= value;
            tree[x].sum = tree[x].maxNum;
            return;
        }

        int m = (lx + rx) / 2;
        updateModEqual(value, l, r, 2 * x + 1, lx, m);
        updateModEqual(value, l, r, 2 * x + 2, m, rx);
        updateFromChildren(x);
    }

    void updateEqual(int index, int64_t value, int x, int lx, int rx)
    {
        if (rx - lx == 1) {
            tree[x].sum = tree[x].maxNum = value;
            return;
        }
        int m = (lx + rx) / 2;
        if (index < m) {
            updateEqual(index, value, 2 * x + 1, lx, m);
        }
        else {
            updateEqual(index, value, 2 * x + 2, m, rx);
        }
        updateFromChildren(x);
    }

    int64_t findSum(int l, int r, int x, int lx, int rx)
    {
        if (lx >= r || rx <= l) {
            return 0LL;
        }
        if (lx >= l && rx <= r) {
            return tree[x].sum;
        }
        int m = (lx + rx) / 2;
        return findSum(l, r, 2 * x + 1, lx, m) + findSum(l, r, 2 * x + 2, m, rx);
    }

public:
    JiDriverRemainder(const std::vector<int64_t> &v)
    {
        init(v.size());
        build(0, 0, size, v);
    }

    void updateModEqual(int64_t value, int l, int r)
    {
        updateModEqual(value, l, r, 0, 0, size);
    }

    void updateEqual(int index, int64_t value)
    {
        updateEqual(index, value, 0, 0, size);
    }

    int64_t findSum(int l, int r)
    {
        return findSum(l, r, 0, 0, size);
    }
};

class RsqSegtree
{
    std::vector<int64_t> tree;
    int size;
    void init(int n)
    {
        size = 1;
        while (size < n) size *= 2;
        tree.assign(2 * size - 1, 0LL);
    }
    void set(int i, int64_t v, int x, int lx, int rx)
    {
        if (rx - lx == 1) {
            tree[x] = v;
            return;
        }
        int m = (lx + rx) / 2;
        if (i < m) {
            set(i, v, 2 * x + 1, lx, m);
        }
        else {
            set(i, v, 2 * x + 2, m, rx);
        }
        tree[x] = tree[2 * x + 1] + tree[2 * x + 2];
    }
    int64_t sum(int l, int r, int x, int lx, int rx)
    {
        if (l >= rx || lx >= r) return 0;
        if (lx >= l && rx <= r) return tree[x];
        int m = (lx + rx) / 2;
        int64_t sum1 = sum(l, r, 2 * x + 1, lx, m);
        int64_t sum2 = sum(l, r, 2 * x + 2, m, rx);
        return sum1 + sum2;
    }

    void build(int x, int lx, int rx, const std::vector<int64_t> &v)
    {
        if (rx - lx == 1) {
            if (lx < v.size())
                tree[x] = v[lx];
            return;
        }
        int m = (lx + rx) / 2;
        build(2 * x + 1, lx, m, v);
        build(2 * x + 2, m, rx, v);
        tree[x] = tree[2 * x + 1] + tree[2 * x + 2];
    }
public:
    void set(int i, int64_t v)
    {
        set(i, v, 0, 0, size);
    }

    int64_t sum(int l, int r)
    {
        return sum(l, r, 0, 0, size);
    }

    RsqSegtree(const std::vector<int64_t> &v)
    {
        init(v.size());
        build(0, 0, size, v);
    }
};

class RmqSegtree
{
    std::vector<int64_t> tree;
    const int64_t INF = std::numeric_limits<int64_t>::max() / 2;
    int size;
    void init(int n)
    {
        size = 1;
        while (size < n) size *= 2;
        tree.assign(2 * size - 1, INF);
    }
    void set(int i, int64_t v, int x, int lx, int rx)
    {
        if (rx - lx == 1) {
            tree[x] = v;
            return;
        }
        int m = (lx + rx) / 2;
        if (i < m) {
            set(i, v, 2 * x + 1, lx, m);
        }
        else {
            set(i, v, 2 * x + 2, m, rx);
        }
        tree[x] = std::min(tree[2 * x + 1], tree[2 * x + 2]);
    }
    int64_t min(int l, int r, int x, int lx, int rx)
    {
        if (l >= rx || lx >= r) return INF;
        if (lx >= l && rx <= r) return tree[x];
        int m = (lx + rx) / 2;
        int64_t min1 = min(l, r, 2 * x + 1, lx, m);
        int64_t min2 = min(l, r, 2 * x + 2, m, rx);
        return std::min(min1, min2);
    }

    void build(int x, int lx, int rx, const std::vector<int64_t> &v)
    {
        if (rx - lx == 1) {
            if (lx < v.size())
                tree[x] = v[lx];
            return;
        }
        int m = (lx + rx) / 2;
        build(2 * x + 1, lx, m, v);
        build(2 * x + 2, m, rx, v);
        tree[x] = std::min(tree[2 * x + 1], tree[2 * x + 2]);
    }
public:
    void set(int i, int64_t v)
    {
        set(i, v, 0, 0, size);
    }

    int64_t min(int l, int r)
    {
        return min(l, r, 0, 0, size);
    }

    RmqSegtree(const std::vector<int64_t> &v)
    {
        init(v.size());
        build(0, 0, size, v);
    }
};

class SegtreeEqual
{
    struct Node
    {
        int64_t value;
        bool flag = false;
    };
    std::vector<Node> tree;
    int size;

    void init(int n)
    {
        size = 1;
        while (size < n) size *= 2;
        tree.assign(size * 2 - 1, {0, false});
    }

    inline void propagate(int x)
    {
        if (!tree[x].flag) return;
        tree[2 * x + 2] = tree[2 * x + 1] = tree[x];
        tree[x].flag = false;
    }
    void set(int l, int r, int64_t v, int x, int lx, int rx)
    {
        if (tree[x].flag && rx - lx != 1) propagate(x);
        if (lx >= l && rx <= r) {
            tree[x].value = v;
            tree[x].flag = true;
            return;
        }
        if (lx >= r || rx <= l) return;

        int m = (lx + rx) / 2;
        set(l, r, v, 2 * x + 1, lx, m);
        set(l, r, v, 2 * x + 2, m, rx);

    }
    int64_t get(int i, int x, int lx, int rx)
    {
        if (rx - lx != 1) propagate(x);
        else return tree[x].value;

        int m = (lx + rx) / 2;
        if (i < m) return get(i, 2 * x + 1, lx, m);
        else return get(i, 2 * x + 2, m, rx);
    }
    void build(int x, int lx, int rx, const std::vector<int64_t> &v)
    {
        if (rx - lx == 1) {
            if (lx < v.size()) {
                tree[x].value = v[lx];
                tree[x].flag = false;
            }
            return;
        }
        int m = (lx + rx) / 2;
        build(2 * x + 1, lx, m, v);
        build(2 * x + 2, m, rx, v);
    }
public:
    SegtreeEqual(const std::vector<int64_t> &v)
    {
        init(v.size());
        build(0, 0, size, v);
    }

    void set(int l, int r, int64_t value)
    {
        set(l, r, value, 0, 0, size);
    }

    int64_t get(int index)
    {
        return get(index, 0, 0, size);
    }
};
}
vector<int> z_function(string s)
{
    int n = (int) s.length();
    vector<int> z(n);
    for (int i = 1, l = 0, r = 0; i < n; ++i) {
        if (i <= r)
            z[i] = min(r - i + 1, z[i - l]);
        while (i + z[i] < n && s[z[i]] == s[i + z[i]])
            ++z[i];
        if (i + z[i] - 1 > r)
            l = i, r = i + z[i] - 1;
    }
    return z;
}

vector<int> prefix_function(string s)
{
    int n = (int) s.length();
    vector<int> pi(n);
    for (int i = 1; i < n; ++i) {
        int j = pi[i - 1];
        while (j > 0 && s[i] != s[j])
            j = pi[j - 1];
        if (s[i] == s[j]) ++j;
        pi[i] = j;
    }
    return pi;
}

void sufmas(string s, int n)
{
    const int alphabet = 26;
    const int maxlen = 1e6;
    int p[maxlen], cnt[maxlen], c[maxlen];
    memset(cnt, 0, alphabet * sizeof(int));
    for (int i = 0; i < n; ++i)
        ++cnt[s[i]];
    for (int i = 1; i < alphabet; ++i)
        cnt[i] += cnt[i - 1];
    for (int i = 0; i < n; ++i)
        p[--cnt[s[i]]] = i;
    c[p[0]] = 0;
    int classes = 1;
    for (int i = 1; i < n; ++i) {
        if (s[p[i]] != s[p[i - 1]]) ++classes;
        c[p[i]] = classes - 1;
    }

    int pn[maxlen], cn[maxlen];
    for (int h = 0; (1 << h) < n; ++h) {
        for (int i = 0; i < n; ++i) {
            pn[i] = p[i] - (1 << h);
            if (pn[i] < 0) pn[i] += n;
        }
        memset(cnt, 0, classes * sizeof(int));
        for (int i = 0; i < n; ++i)
            ++cnt[c[pn[i]]];
        for (int i = 1; i < classes; ++i)
            cnt[i] += cnt[i - 1];
        for (int i = n - 1; i >= 0; --i)
            p[--cnt[c[pn[i]]]] = pn[i];
        cn[p[0]] = 0;
        classes = 1;
        for (int i = 1; i < n; ++i) {
            int mid1 = (p[i] + (1 << h)) % n, mid2 = (p[i - 1] + (1 << h)) % n;
            if (c[p[i]] != c[p[i - 1]] || c[mid1] != c[mid2])
                ++classes;
            cn[p[i]] = classes - 1;
        }
        memcpy(c, cn, n * sizeof(int));
    }
}
namespace SuffixAutomaton
{
struct state
{
    int len, link;
    std::map<char, int> next;
};

constexpr int MAXLEN = 100000;

state st[MAXLEN * 2];

int sz, last;

void sa_init()
{
    sz = last = 0;
    st[0].len = 0;
    st[0].link = -1;
    ++sz;
    /*

    for (int i=0; i<MAXLEN*2; ++i)
        st[i].next.clear();
    */
}

void sa_extend(char c)
{
    int cur = sz++;
    st[cur].len = st[last].len + 1;
    int p;
    for (p = last; p != -1 && !st[p].next.count(c); p = st[p].link)
        st[p].next[c] = cur;
    if (p == -1)
        st[cur].link = 0;
    else {
        int q = st[p].next[c];
        if (st[p].len + 1 == st[q].len)
            st[cur].link = q;
        else {
            int clone = sz++;
            st[clone].len = st[p].len + 1;
            st[clone].next = st[q].next;
            st[clone].link = st[q].link;
            for (; p != -1 && st[p].next[c] == q; p = st[p].link)
                st[p].next[c] = clone;
            st[q].link = st[cur].link = clone;
        }
    }
    last = cur;
}
}
}

int main()
{
    std::ios_base::sync_with_stdio(0);
    std::cin.tie(0);
    std::cout.tie(0);
    //обращение к алгосам  icpc::

}
