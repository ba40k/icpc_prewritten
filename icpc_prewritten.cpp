#include <iostream>
#include <vector>
#include <algorithm>
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
 // алгосы тут 
}

int main()
{
    std::ios_base::sync_with_stdio(0);
    std::cin.tie(0);
    std::cout.tie(0);
    //обращение к алгосам  icpc::
    
}
