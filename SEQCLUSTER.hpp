#ifndef SEQCLUSTER_H
#define SEQCLUSTER_H
#include "HELPER.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include <map>
using namespace std;
using namespace helper;

/*
 *  Cluster Info:
 *	c_name = cluster name
 */
struct CINFO
{
	CINFO(string name = "", int start = 0) : c_name(name), c_start(start){}
	string c_name;
	int c_start;
	vector <string> tnames;
};

class SEQCLUSTER
{

	private:
		int num_cluster;
		int num_sequence;
		vector<int> sIndex_lookup;	// sequence id starts from 1
		vector<int> cIndex_lookup;	// cluster id starts from 1
		vector<int> start_pos;		// store the starting position of original sequence
		vector<CINFO> cluster_info;
		string text;

	public:
		// constructor
		SEQCLUSTER(const char *filename);
		
		// function
		int get_num_cluster();
		int get_num_sequence();
		int pos_to_sIndex(int char_pos);
		int sIndex_to_cIndex(int sIndex);
		int get_start_pos(int sIndex);
		int get_transcript_length(int sIndex);
		string get_transcript_name(int sIndex);
		int convert_to_forward_pos(int pos, int sigmer_len);
		string getText();
		CINFO get_cluster_info(int cIndex);
};

#endif



