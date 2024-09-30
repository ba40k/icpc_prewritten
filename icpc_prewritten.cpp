#include <iostream>
#include <vector>
#include <algorithm>
namespace icpc
{
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
                return parent[vertex] = parent[vertex]==vertex?vertex:get_parent(parent[vertex]);
            }
            void set_union(int a, int b)
            {
                a = get_parent(a);
                b = get_parent(b);
                len[a]<len[b]?std::swap(a,b):std::swap(a,a);
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
