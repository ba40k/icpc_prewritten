#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <limits>

namespace icpc
{
    void EuelerTourForTree(int vertex, int parent, int depth, std::vector<std::vector<int>>& graph, std::vector<std::pair<int, int>>& eueler_tour)
    {
        eueler_tour.push_back({ vertex,depth });
        int cur_vertex = vertex;
        for (int next_vertex : graph[vertex])
        {
            if (next_vertex != parent)
            {
                EuelerTourForTree(next_vertex, vertex, depth + 1, graph, eueler_tour);
                eueler_tour.push_back({ vertex,depth });
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
        std::vector<node>  dec;
        std::vector<int> vec;
        void build(std::vector<int>& input)
        {
            vec = input;
            dec.resize((input.size()>>k) + 1);
            for (int i = 0; i < input.size(); ++i)
            {
                // creating
            }
        }
    public:
        SqrtDec(std::vector<int>& input)
        {
            build(input);
        }
        
    };
    class Dsu // базовая дсу 
    {
        private:
            std::vector<int> parent;
            std::vector<int> len;
            void build (int n)
            {
                parent.resize(n);
                len.resize(n,1);
                for (int i =0;i<n;++i) 
                {
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
                if (parent[vertex] == vertex)
                {
                    return vertex;
                }
                 return parent[vertex] = get_parent(parent[vertex]);
            }
            void set_union(int a, int b)
            {
                 a = get_parent(a);
                 b = get_parent(b);
                 if (len[a] < len[b])
                 {
                     len[b] += len[a];
                    parent[a] = b;
                    return ;
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
    }

 // алгосы тут 
}

int main()
{
    std::ios_base::sync_with_stdio(0);
    std::cin.tie(0);
    std::cout.tie(0);
    //обращение к алгосам  icpc::
    
}
