#include "SEQCLUSTER.hpp"
#include "HELPER.hpp"
#include "aho01new.hpp"
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
			vector <int> *int_tids = new vector<int>();
			split(sigmer_data[2], ',', text_tids);

			for(auto t:text_tids){
				if(t != ""){
					int_tids -> push_back(stoi(t));
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
void count_reads(const string readfile, AC &aho_sigmer, SEQCLUSTER* sc, vector<unordered_map<string, double>> &equivalent_class){
	size_t mismatched = 0;
	size_t conflict = 0;
	size_t unmapped = 0;

	// this set of vector needs to be redistributed later
	vector< unordered_map <int, string >> ambiguous_profile_vector;

	ifstream infile(readfile, ifstream::in);
	if(infile.is_open()){
		string line;
		while(getline(infile, line)){
			if(line[0] != '@' && line[0] != '~' && line[0] != '>' && line[0] != '+'){

	//			cout << line << "::";
				// identify the appearance_profile for different cluster
				unordered_map <int, set<int>> appearance_profile;

				// querying reads
				vector<MatchPair> results = aho_sigmer.query(line);
	//			cout << results.size() << endl;

				for(int i = 0; i < results.size(); i++){

					// identify the cluster first for this sigmer
					int first_tid = results[i].value -> at(0);
					int cid =  sc -> sIndex_to_cIndex(first_tid * 2);


//					for(int j = 0; j < results[i].value -> size(); j++)
//						cout << results[i].value -> at(j) << ",";
//					cout << endl;
	
					// build transcript profile
					if(appearance_profile.find(cid) == appearance_profile.end()){
						set <int> tmp;
						tmp.insert(results[i].value -> begin(), results[i].value -> end());
						appearance_profile[cid] = tmp;
					}
					else{
						set <int> tmp;
						set_intersection(appearance_profile[cid].begin(), appearance_profile[cid].end(), results[i].value->begin(), results[i].value->end(), inserter(tmp, tmp.begin()));
						appearance_profile[cid] = tmp;
					}

				}
	
				if(results.size() == 0){
					unmapped ++;
				}

				// record ambiguous profile
				unordered_map <int, string> ambiguous_profile;

				// iterate through transcript profile
				for(auto p: appearance_profile){

					int cid = p.first;
					CINFO c_info = sc -> get_cluster_info(cid);
					size_t csize = c_info.tnames.size();
			
					if(p.second.size() == 0)	// debug
						conflict ++;
		
					// make Ymatrix in string fromat
					string read_profile = "";
					for (size_t i = 0; i < csize; i++)
						read_profile += "0";
					for (auto t: p.second)
						read_profile[(t*2 - c_info.c_start)/2] = '1';


					// if sigmers of this read are all in the same cluster
					// collapse read profile into equivalent class, with equivalent class count
					if(appearance_profile.size() == 1){
						if(equivalent_class[cid].find(read_profile) == equivalent_class[cid].end())
							equivalent_class[cid][read_profile] = 1;
						else
							equivalent_class[cid][read_profile] += 1;

//						cout  << cid << "\t" << read_profile << endl;
					}
					// if this read is ambiguous, put aside the information
					else{
						ambiguous_profile[cid] = read_profile;
					}	
				}
				if(appearance_profile.size() > 1){
					mismatched++;
					ambiguous_profile_vector.push_back(ambiguous_profile);
				}
			}
		}

		// redistribute ambiguous reads
		for(size_t i = 0; i < ambiguous_profile_vector.size(); i++){
			unordered_map <int, string> ambiguous_profile = ambiguous_profile_vector[i];
			double total_count = 0;
			for(auto p: ambiguous_profile){
				if(equivalent_class[p.first].find(p.second) != equivalent_class[p.first].end())
					total_count += equivalent_class[p.first][p.second];	
			}

			for(auto p:ambiguous_profile){
				if(equivalent_class[p.first].find(p.second) != equivalent_class[p.first].end())
					equivalent_class[p.first][p.second] += equivalent_class[p.first][p.second] / total_count;
			}
			
		}
	}


	infile.close();

	cout << "Unmapped Reads: " << unmapped << endl;
	cout << "Conflict Reads: " << conflict << endl;
	cout << "Mismatched Reads: " << mismatched << endl;
}

void export_observation(const vector<unordered_map<string, double>> &Y, const string &outfile){
	
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
	vector <unordered_map <string, double>> equivalent_class_count(sc -> get_num_cluster(), unordered_map<string, double>()); // condensed version of Ymatrix, separated by cluster


	count_reads(read1, aho_sigmer, sc, equivalent_class_count);

	echo("Export Observation Matrix");
	export_observation(equivalent_class_count, outfile);
	
	echo("Done");

}
