#include "SEQCLUSTER.h"
#include "HELPER.h"
#include <iostream>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <set>
#include <tuple>
#include <unordered_map>

using namespace std;
using namespace helper;

struct sigmer_metadata{
	int sigmer_size;
	set<int> positions;
//	int uniqueness; // record the number of transcripts containing this sigmer
	set<int> tids;
	bool selected = false;
};

bool cmp (const sigmer_metadata *a, const sigmer_metadata *b){

	return a -> tids.size() < b -> tids.size() ||
		(a -> tids.size() == b -> tids.size()  && a -> sigmer_size < b -> sigmer_size);
}



/*
 * Split the sigmer into different files based on their cluster ids
 * Filter the sigmer based on the boundary
 */

void filter_sigmer(ofstream* outarray, const char* sigmerfile, int min_size, int max_size, int core)
{

	//reading sigmerfile
	ifstream infile(sigmerfile, ifstream::in);
	if(infile.is_open())
	{
		string line;
		while(getline(infile, line))
		{
			vector<string> sigmer;
			split(line, '\t', sigmer);
			int sigmer_min_len = atoi(sigmer[2].c_str());
			int sigmer_max_len = atoi(sigmer[3].c_str());
			int rlen = sigmer_min_len;

			// look for minimum length
			if(rlen < min_size)
				rlen = min_size;

                        // the sigmer size has to be within the sigmer len boundary
                        // it also needs to be smaller than the upper bound (which is the read length)
                        if(rlen <= sigmer_max_len && rlen <= max_size)
			{
				int index = atoi(sigmer[0].c_str()) % core;
				outarray[index] << sigmer[0] << "\t" << sigmer[1] << "\t" << rlen << endl;
			}			

		}

	}
	infile.close();

}


/*      Read in all the sigmers, and normalize the sigmer starting positions
 *	Final structure of sigmer_info:
 * 	[ [size, pos1, pos2...], [size, pos1, pos2...] ... ]
 * 
 */

void build_sigmer_info(vector <sigmer_metadata> &sigmer_info, const string sigmerfile, SEQCLUSTER* scptr)
{
	map <string, int> sigmer_id;  // use to record the sigmer id in sigmer_info

	// reading sigmerfile
	ifstream infile(sigmerfile, ifstream::in);
        if(infile.is_open())
        {
		string line;
		int id = 0;
		while(getline(infile, line))
		{
			vector<string> sigmer;
			split(line, '\t', sigmer);
			int rlen = atoi(sigmer[2].c_str());

			vector<string> start_positions;
			split(sigmer[1], ',', start_positions);

			string seq = scptr -> getText().substr(atoi(start_positions[0].c_str()), rlen);

			// convert positions from string to int
			// convert positions from complementary sequence to forward sequence
			set<int> info;
			for(int s = 0; s < start_positions.size() -1; s++){ // ignore last empty number
				int pos = atoi(start_positions[s].c_str());
				pos = scptr -> convert_to_forward_pos(pos, rlen);
				info.insert(pos);
			}

			string rc_seq = seq;
			reverse_complement(rc_seq);

			// check if this sequence has already been processed before
			if( sigmer_id.find(seq) != sigmer_id.end()) // the forward sequence is already in the database
			{
				int s_id = sigmer_id[seq];
				sigmer_info[s_id].positions.insert(info.begin(), info.end());
			}
			else if ( sigmer_id.find(rc_seq) != sigmer_id.end()) // the reverse complement sequence is already in the database
			{
				int s_id = sigmer_id[rc_seq];
				sigmer_info[s_id].positions.insert(info.begin(), info.end());
			}
			else	// not seen before
			{
				sigmer_metadata sm;
				sm.sigmer_size = rlen;
				sm.positions = info;
				sigmer_info.push_back(sm);
				sigmer_id[seq] = id;
				id++;
			}
		}
	}
	
	infile.close();
}

/*
 *	Function : determine the transcripts each sigmer covers
 *
 */
void sigmer_uniqueness(vector <sigmer_metadata> &sigmer_info, SEQCLUSTER *sc)
{
	for(int sigmer_id = 0; sigmer_id < sigmer_info.size(); sigmer_id++)
	{
		set <int> seq_ids;
		for(auto pos : sigmer_info[sigmer_id].positions)
			seq_ids.insert((sc->pos_to_sIndex(pos)) / 2);
		sigmer_info[sigmer_id].tids = seq_ids;
	}
}


/*
 *  Construct the transcript map:
 *  	each key is position, and values are the sigmer pointer starting at this position
 *	{key -> [sid_1, sid_2, sid_3, ...] }
 */

void build_transcript_map( map <int, vector<sigmer_metadata*>> *transcript_map, vector <sigmer_metadata> &sigmer_info)
{
        for (int i = 0; i < sigmer_info.size(); i++)
        {
                set <int> start_positions = sigmer_info[i].positions;
                for(auto pos : start_positions)
                {
                        if( transcript_map -> find(pos) != transcript_map -> end())
                                transcript_map->operator[](pos).push_back(&sigmer_info[i]);
                        else
                        {
                                vector <sigmer_metadata*> tmp_ids;
                                tmp_ids.push_back(&sigmer_info[i]);
                                transcript_map->operator[](pos) = tmp_ids;
                        }
			sort(transcript_map->operator[](pos).begin(), transcript_map -> operator[](pos).end(), cmp);
                }
        }
}

/* 
 * Function: select the most amount of subsets that represent the biggest set
 */
int select_sigmers( vector<sigmer_metadata*> current_sigmer_ptrs, unordered_map <int, bool> &covered_positions){

	int num_selection = 0;
	set<int> max_tids = current_sigmer_ptrs[current_sigmer_ptrs.size() -1] -> tids;

	for(size_t i = 0; i < current_sigmer_ptrs.size() && max_tids.size() > 0; i++){
		sigmer_metadata* examine_sigmer = current_sigmer_ptrs[i];

		for(auto t:examine_sigmer -> tids){
			if( max_tids.find(t) != max_tids.end()){
				examine_sigmer -> selected = true;
				num_selection++;

				// debug
//				cout << "SELECTED " << *(examine_sigmer -> positions.begin()) << "\t" << examine_sigmer -> sigmer_size << ",";
				
				// marked positions with selected sigmers
				for(auto p:examine_sigmer -> positions)
					covered_positions[p] = true;

				set<int> tmp;
				set_difference(max_tids.begin(), max_tids.end(), examine_sigmer -> tids.begin(), examine_sigmer -> tids.end(), inserter(tmp, tmp.end()));
				max_tids = tmp;
				break;
			}
		}
	}
	return num_selection;
}

/*
 * Function: Traverse the transcript map to pick sig-mer
 * At each position, sig-mer is selected based on its uniquness and size
 * i) if the current position is covered, skipped
 * ii) if the current position is not yet covered, and away from gap, go through the selection process
 * return the number of selected sigmer
 */

int traverse_transcript_map(map <int, vector<sigmer_metadata*>> *transcript_map, int gap, SEQCLUSTER* sc)
{
	int num_selected_sigmers = 0;
	int previous_uniqueness = 0;
	int previous_selected = 0;
//	int previous_fork = 0;
        unordered_map <int, bool> covered_positions; // true indicate sigmer has been selected for this position
	
	for(auto traversal_pos = transcript_map -> begin(); traversal_pos != transcript_map -> end(); ++traversal_pos){
		int current_pos = traversal_pos -> first; // get position
		vector <sigmer_metadata*> current_sigmer_ptrs = traversal_pos -> second;
		int max_uniqueness = current_sigmer_ptrs[current_sigmer_ptrs.size()-1] -> tids.size();
		int candidate = (covered_positions.find(current_pos) == covered_positions.end()); // if this pos hasn't been visited

/*		// debugging
		cout << " at position " << current_pos << "::";
                for(int i = 0; i < current_sigmer_ptrs.size(); i++){
                	sigmer_metadata* examine_sigmer = current_sigmer_ptrs[i];
			cout << "(size = " << examine_sigmer -> sigmer_size << ";";
			for(auto p: examine_sigmer -> positions){
				cout << p << ",";
			}
			cout << ")";
		}

		cout << ":::" << sc -> getText().substr(current_pos, current_sigmer_ptrs[0] -> sigmer_size) <<  " "; 

		if(!candidate){
			cout << " COVERED ";
			if(covered_positions[current_pos]){
				cout << " P-SELECTED ";
			}
		}
		// end of debugging
*/
		// mark the position visisted
		if(candidate){
			for(auto p: (current_sigmer_ptrs[current_sigmer_ptrs.size()-1] -> positions)){
				covered_positions[p] = false;
			}
		}


		// add sigmer at fork (remove the start of a transcript))
		if(current_pos != 0 && sc->pos_to_sIndex(current_pos) / 2 == sc -> pos_to_sIndex(current_pos - 1) / 2 && max_uniqueness != previous_uniqueness){
//			cout << " NODE " ;	
			map <int, vector<sigmer_metadata*>>::iterator backtrack_it = traversal_pos;
			for(size_t i = 0; i < 5 && backtrack_it != transcript_map -> begin() ; i++){
				int backtrack_pos = backtrack_it -> first;
				
				// if it hits the start of the position, skip the loop
				if(backtrack_pos ==0 || sc -> pos_to_sIndex(backtrack_pos) / 2 != sc -> pos_to_sIndex(backtrack_pos - 1) / 2)
					break;
				
				// if the position hasn't been selected yet
				else if(covered_positions[backtrack_pos] == false){
					num_selected_sigmers += select_sigmers(backtrack_it -> second, covered_positions);
				}
				backtrack_it--;
			}
		}

		// add sigmer for
		// 1. start position of a transcript
		// 2. supplementary sigmer to fill the gap
		else if( (current_pos == 0 || sc->pos_to_sIndex(current_pos) / 2 != sc -> pos_to_sIndex(current_pos - 1) / 2 || current_pos >= previous_selected + gap) && candidate){
			
			num_selected_sigmers += select_sigmers(current_sigmer_ptrs, covered_positions);
		}

		// update previous_selected
		if(covered_positions[current_pos])
			previous_selected = current_pos;

		previous_uniqueness = max_uniqueness;
//		cout << endl;
	}

	return num_selected_sigmers;
}


void export_sigmer_detail(vector <sigmer_metadata> &sigmer_info, SEQCLUSTER *sc, ofstream& outfh_sigmer){

	// iterate through sigmers
	for(auto s: sigmer_info){
		if(s.selected){
			int sigmer_start = *(s.positions.begin());
			string seq = sc -> getText().substr(sigmer_start, s.sigmer_size);
			int t_id = sc -> pos_to_sIndex(sigmer_start);
			int c_id = sc -> sIndex_to_cIndex(t_id);
			CINFO cluster_info = sc -> get_cluster_info(c_id);
			outfh_sigmer << cluster_info.c_name << "\t" << seq << "\t" ;

			unordered_map <int, int> tids;
			// iterate through all positions
			for (auto start_p : s.positions){
				t_id = sc -> pos_to_sIndex(start_p);
		
				// sanity check
				if( sc -> sIndex_to_cIndex(t_id) != c_id)
					cerr << "Inconsistent Cluster: " << c_id << " vs " << sc -> sIndex_to_cIndex(t_id) << "\t" << seq << endl;
				
				if(tids.find(t_id) == tids.end())
					tids[t_id] = 0;
			}
			for (auto t : tids)
			{
				int transcript_id = t.first;
//				transcript_id = (transcript_id - cluster_info.c_start) / 2;
				outfh_sigmer << transcript_id / 2 << ",";
			}
			outfh_sigmer << endl;

		}
	}
}


int main(int argc, char* argv[])
{
        if(argc < 6)
        {
                cerr << "Usage: " << argv[0] << " SEQFILE SIGMERFILE READLEN OUTPREFIX CORE" << endl;
                return 1;
        }

	const char *seqfile = argv[1];
	const char *sigmerfile = argv[2];
	const int readlen = atoi(argv[3]);
	const char *outdir = argv[4];
	const int core = atoi(argv[5]);
	int min_len = 25;
	int max_len = 0.8*readlen;
	
	SEQCLUSTER *sc = new SEQCLUSTER(seqfile);
	
	// sigmer partition
	ofstream* outarray = new ofstream[core];	
	for(int i = 0; i < core; i++)
	{
		string filename = string(outdir) + "_sigmers_" + to_string(i) + ".txt";
		outarray[i].open(filename, ios_base::out);
	}

	echo("Filter sigmers");
	filter_sigmer(outarray, sigmerfile, min_len, max_len, core);

	for(int i = 0; i < core; i++)
		outarray[i].close();

	// sigmer selection
	int gap = 30;
	int total_sigmer_count = 0;
	int selected_sigmer_count = 0;

	ofstream sigmerfh;
	sigmerfh.open(string(outdir) + "_selected_sigmers.txt", ios_base::out);
	for(int i = 0; i < core; i++)
	{
                echo("Building sigmers " + to_string(i));
		vector <sigmer_metadata> sigmer_info;
		build_sigmer_info(sigmer_info, string(outdir)+"_sigmers_"+to_string(i)+".txt", sc);
		total_sigmer_count += sigmer_info.size();

		echo("Scanning sigmer uniqueness " + to_string(i));
		sigmer_uniqueness(sigmer_info, sc);

	        echo("Building transcript map");
		map <int, vector <sigmer_metadata*> >* transcript_map = new map <int, vector<sigmer_metadata*>>();
	        build_transcript_map(transcript_map, sigmer_info);

		echo("Traversing transcript map");
		int s = traverse_transcript_map(transcript_map, gap, sc);
		selected_sigmer_count += s;
		delete transcript_map;

		echo("Exporting selected sigmers");
		export_sigmer_detail(sigmer_info, sc, sigmerfh);

	}
	sigmerfh.close();

	cout << "Total sigmer count " << total_sigmer_count << endl;
	cout << "Selected sigmer count " << selected_sigmer_count << endl;

	delete sc;

	echo("Done");

}

