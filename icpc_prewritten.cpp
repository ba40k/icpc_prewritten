#include <iostream>
#include <vector>

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
 // алгосы тут 
}

int main()
{
    //обращение к алгосам  icpc::
    
}
