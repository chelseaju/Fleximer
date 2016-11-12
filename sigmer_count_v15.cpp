#include "SEQCLUSTER.h"
#include "HELPER.h"
#include "aho01.hpp"
#include <iostream>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <unordered_map>
using namespace std;
using namespace helper;


/*
 * Read in sigmer file and put it in hash
 */
void load_sigmers(const string file, AC &aho_sigmer){

	// read in sigmer
	ifstream infile(file, ifstream::in);
	if(infile.is_open()){
		string line;
		while(getline(infile, line)){
			// load sequence
			vector <string> sigmer_data;
			split(line, '\t', sigmer_data);
			string seq = sigmer_data[1];

			// extract transcript ids
			vector <string> text_tids;
//			set <int> int_tids;
			vector<int> *int_tids = new vector<int>();

			split(sigmer_data[2], ',', text_tids);

			for(auto t:text_tids){
				if(t != ""){
					int_tids -> push_back(stoi(t));
//					int_tids.insert(stoi(t));
				}
			}

			sort(int_tids->begin(), int_tids->end());
			// put sigmer to tries
			aho_sigmer.add_pattern(seq, int_tids);
			reverse_complement(seq);
			aho_sigmer.add_pattern(seq, int_tids);
			
		}
	}
	infile.close();
	aho_sigmer.construct();
}

/*
 * Go through each read sequence, use sigmer information to determine the origin of each read
 */
void count_reads(const string readfile, AC &aho_sigmer, SEQCLUSTER* sc, vector<unordered_map<string, int>> &equivalent_class){
	ifstream infile(readfile, ifstream::in);
	if(infile.is_open()){
		string line;
		while(getline(infile, line)){
			if(line[0] != '@' && line[0] != '~' && line[0] != '>' && line[0] != '+'){

				set <int> appearance_profile;

				// querying reads
				vector<MatchPair> results = aho_sigmer.query(line);

				for(size_t i = 0; i < results.size(); i++){
					
					if(i == 0)
						appearance_profile.insert(results[i].value -> begin(), results[i].value -> end());
					else{
						set <int> tmp;
						set_intersection(appearance_profile.begin(), appearance_profile.end(), results[i].value -> begin(), results[i].value -> end(), inserter(tmp, tmp.begin()));
						appearance_profile = tmp;
					}
					
					if(appearance_profile.empty())
						break;

				}

				// get cluster size

				int cid = sc -> sIndex_to_cIndex( *(appearance_profile.begin()) * 2);
				CINFO c_info = sc -> get_cluster_info(cid);
				int csize = c_info.tnames.size();

				// make Ymatrix in string format
				string read_profile = "";
				for (int i = 0; i < csize; i++)
					read_profile += "0";


				for (auto t: appearance_profile)
					read_profile[(t*2 - c_info.c_start)/2] = '1';
				

				// collapse read profile into equivalent class, with equivalent class count
				if(equivalent_class[cid].find(read_profile) == equivalent_class[cid].end())
					equivalent_class[cid][read_profile] = 1;
				else
					equivalent_class[cid][read_profile] += 1;

			}
		}
	}
	infile.close();


}

void export_observation(const vector<unordered_map<string, int>> &Y, const string &outfile){
	
	ofstream outfh;
        outfh.open(outfile, ios_base::out);

	for(int c = 0; c < Y.size(); c++){
		for(auto profile:Y[c]){
			outfh << c << "\t" << profile.first << "\t" << profile.second << endl;
		}
	}
	outfh.close();
}

int main(int argc, char* argv[])
{
        if(argc < 4)
        {
                cerr << "Usage: " << argv[0] << " SEQFILE SIGMERFILE OUTFILE READ1 READ2" << endl;
                return 1;
        }

	const char *seqfile = argv[1];
        const char *sigmer_file = argv[2];
	const char *outfile = argv[3];
        const char *read1 = argv[4];

	SEQCLUSTER *sc = new SEQCLUSTER(seqfile);

        echo("Loading sigmers");
	AC aho_sigmer;
	load_sigmers(sigmer_file, aho_sigmer);

	echo("Counting reads");
	vector <unordered_map <string, int>> equivalent_class_count(sc -> get_num_cluster(), unordered_map<string, int>()); // condensed version of Ymatrix, separated by cluster


	count_reads(read1, aho_sigmer, sc, equivalent_class_count);

	echo("Export Observation Matrix");
	export_observation(equivalent_class_count, outfile);
	
	echo("Done");

}
