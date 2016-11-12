#include <string>
#include <iostream>
#include <map>
#include <queue>
#include <vector>
#include <utility>
#include <set>
#include "aho01new.hpp"

int num_nodes = 0;

inline int char_to_idx(char ch){
    if(ch == 'A') return 1;
    else if(ch == 'C') return 2;
    else if(ch == 'G') return 3;
    else if(ch == 'T') return 4;
    else return 0;
}


inline int char_to_01(char ch){
    if(ch == 'A') return 0;
    else if(ch == 'C') return 0;
    else if(ch == 'G') return 1;
    else if(ch == 'T') return 1;
    else return 0;
}



inline char idx_to_char(int idx){
    if(idx == 0) return 'A';
    else if(idx == 1) return 'C';
    else if(idx == 2) return 'G';
    else if(idx == 3) return 'T';
    return 'Z';
}


int AC::get_num_patterns(){ return num_patterns; }

//void AC::add_pattern(const std::string &p, const std::set<int> &v){
void AC::add_pattern(const std::string &p, std::vector<int> *v){

    Node *cur = root;
    
    long long int h = 0;    
//    long long int h2 = 0;    
    for(std::size_t i = 0; i < p.length(); i++){
        cur = cur->build(char_to_01(p[i]));
        h = h * HASH_BASE + char_to_idx(p[i]);
//        h2 = h2 * HASH_BASE2 + char_to_idx(p[i]);
    }
    
    if(cur->hash == NULL){
        cur->hash = new std::map<long long int, int>;
//        cur->hash2 = new std::map<long long int, int>;
    }

    cur->hash->insert(std::pair<long long int, int>(h, num_patterns));
//    cur->hash2->insert(std::pair<long long int, int>(h2, num_patterns));
        

//    pattern_hash.push_back(h);
    pattern_value.push_back(v);
    ++num_patterns;
}


void AC::construct(){

    std::cerr << num_nodes << " nodes." << std::endl;

    long long int h = 1;
//    long long int h2 = 1;
    for(int i = 0; i < 1024; i++){
        base_level.push_back(h);
//        base_level2.push_back(h2);
        h *= HASH_BASE;
//        h2 *= HASH_BASE2;
    }


    std::queue<Node*> que;

    que.push(root);

    Node *cur, *nxt, *ptr;

    while(!que.empty()){
        cur = que.front();
        que.pop();
        
        for(int i = 0; i < 2; i++){
            if(cur->_next[i] == NULL){
                continue;
            }
            nxt = cur->_next[i];

            if(cur == root){
                nxt->fail = root;
            }else{
                ptr = cur->fail;
                while(ptr != NULL){
                    if(ptr->_next[i] != NULL){
                        nxt->fail = ptr->_next[i];
                        break;
                    }
                    ptr = ptr->fail;
                }
                if(ptr == NULL){
                    nxt->fail = root;
                }
            }
            que.push(nxt);

        }

    }

}


std::vector<MatchPair> AC::query(const std::string &text){
    std::vector<MatchPair> output;

    long long int h = 0, th;
//    long long int h2 = 0, th2;
    // std::vector<long long int> hasht;
    //hasht.clear();
    while(hasht.size() < text.length()){
        hasht.push_back(0);
    }


    //    hasht2.clear();

    Node *cur = root;
    Node *ptr;

    std::size_t TL = text.length();
    for(std::size_t i = 0; i < TL; i++){
        h = h * HASH_BASE + char_to_idx(text[i]);
//        h2 = h2 * HASH_BASE2 + char_to_idx(text[i]);
        // hasht.push_back(h);
        hasht[i] = h;
//        hasht2.push_back(h2);

        int idx = char_to_01(text[i]);
        while(cur->next(idx) == NULL && cur != root){
            cur = cur->fail;
        }

        if(cur->next(idx) == NULL){
            continue;
        }

        cur = cur->next(idx);

        if(cur->hash != NULL){
            ptr = cur;
            while(ptr != root){
                th = h;
    //            th2 = h2;
                if((int)i - (ptr->depth) >= 0){
                    th = (th - hasht[i - (ptr->depth)] * base_level[ptr->depth]);
    //                th2 = (th2 - hasht2[i - (ptr->depth)] * base_level2[ptr->depth]);
                }

                //if(ptr->hash != NULL && ptr->hash->count(th) > 0){// && ptr->hash2->count(th2) > 0){
                if(ptr->hash != NULL && ptr->hash->find(th) != ptr->hash->end()){// && ptr->hash2->count(th2) > 0){
                    output.push_back(MatchPair(i, pattern_value[ptr->hash->find(th)->second]));
                }
                ptr = ptr->fail;
            }
        }

    }
    
    return output;
}



Node::~Node(){
    if(hash != NULL){
        delete hash;
//        delete hash2;
    }
    for(int i = 0; i < 2; i++){
        if(_next[i] != NULL){
            delete _next[i];
        }
    }
}

Node* Node::build(int idx){
//Node* Node::build(char ch){
//    int idx = char_to_idx(ch);
    if(_next[idx] == NULL){
        _next[idx] = new Node;
        _next[idx]->depth = depth + 1;
        ++num_nodes;
    }
    return _next[idx];
}

Node* Node::next(int idx){
//Node* Node::next(char ch){
//    int idx = char_to_idx(ch);
    return _next[idx];
}

