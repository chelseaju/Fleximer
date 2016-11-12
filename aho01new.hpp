#pragma once
#include <map>
#include <utility>
#include <vector>
#include <string>
#include <set>

#define HASH_BASE 7
#define HASH_BASE2 5

class MatchPair{
    public:
        std::size_t pos;
        // std::set<int> *value;
        std::vector<int> *value;
        MatchPair(std::size_t p, std::vector<int> *v): pos(p), value(v){};   // std::set<int>(v)){};
};


class Node{
    public:
        Node *fail;
        Node* _next[2];
        std::map<long long int, int> *hash;
//        std::map<long long int, int> *hash2;
        int depth;

        ~Node();
        Node(): fail(NULL), hash(NULL), depth(0){
            _next[0] = _next[1]  = NULL;
        };

        Node* build(int idx);
        Node* next(int idx);
};


class AC{
    public:
        AC(): num_patterns(0), root(new Node){};

        void add_pattern(const std::string &p, std::vector<int> *v);
       //  void add_pattern(const std::string &p, const std::set<int> &v);
        void construct();
    
        std::vector<MatchPair> query(const std::string &text);
        int get_num_patterns();

    private:
        int num_patterns;
        Node *root;
        std::vector<long long int> base_level;
//        std::vector<long long int> base_level2;
//        std::vector<long long int> pattern_hash;
        std::vector< std::vector<int>* > pattern_value;
        std::vector<long long int> hasht;
//        std::vector<long long int> hasht2;
};
