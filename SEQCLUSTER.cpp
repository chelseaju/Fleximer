#include "SEQCLUSTER.hpp"

// constructor
SEQCLUSTER::SEQCLUSTER(const char *filename)
{
	string line;
	
	// initialize variable
	this -> num_cluster = 0;
	this -> num_sequence = 0;
	int char_index = 0;

	// reading files
	ifstream infile(filename, ifstream::in);
	if(infile.is_open())
	{
		int current_start = 0;
		int current_end = 0;		
		vector <string> cname_tname;	
		while(getline(infile, line))
		{
			if(line.substr(0,1) == ">")
			{
				cname_tname.clear();
				split(line, '|', cname_tname);
				this -> num_cluster ++;
			}
			else{
				// break the sequnces in a cluster
				// each sequence is separated by "|"
				stringstream ss(line);
				string tok;
				
				while(getline(ss, tok, '|'))
				{
					// store start positions
					this -> start_pos.push_back( this -> sIndex_lookup.size());
			
					// forward sequence
					upper(tok);
					this -> text += tok;
					this -> text += "$";
					this -> cIndex_lookup.push_back(this -> num_cluster-1); // push in cluster id

					// push in sequence id
					for(char_index; char_index < text.size(); char_index++)
						this -> sIndex_lookup.push_back(this -> num_sequence);
					this -> num_sequence++;	

					// reverse sequence
					reverse_complement(tok);
					this -> text += tok;
					this -> text += "$";
					this -> cIndex_lookup.push_back(this -> num_cluster-1); // push in cluster id
					
					// push in sequence id
					for(char_index; char_index < text.size(); char_index++)
						this -> sIndex_lookup.push_back(this -> num_sequence);
					this -> num_sequence++;
				}
				current_end = this -> num_sequence -1;

				// store cluster info
				CINFO cluster_detail = CINFO(cname_tname[0].substr(1), current_start);
				vector<string> tname_vec(cname_tname.begin()+1, cname_tname.end());
				cluster_detail.tnames = tname_vec;
				this -> cluster_info.push_back(cluster_detail);	
				current_start = current_end + 1;
			}
		}
	}
}

// return total number of cluster
int SEQCLUSTER::get_num_cluster(){
	return this -> num_cluster;
}

// return total number of sequence
int SEQCLUSTER::get_num_sequence(){
	return this -> num_sequence;
}

// for a given position in text, return its cluster id
int SEQCLUSTER::pos_to_sIndex(int char_pos){
	if(char_pos >= text.size() || char_pos < 0)
		return -1;
	else
		return this -> sIndex_lookup[char_pos];
}

// for a given sequence id, return its cluster id
// sIndex starts with 1
int SEQCLUSTER::sIndex_to_cIndex(int sIndex){
	if(sIndex >= cIndex_lookup.size() || sIndex < 0)
		return -1;
	else
		return this -> cIndex_lookup[sIndex];
}

// for a given sequence id, return its starting position
int SEQCLUSTER::get_start_pos(int sIndex){
	if( (sIndex / 2) > start_pos.size())
		return -1;
	else
		return this -> start_pos[sIndex/2];
}

// for a given sequence id, return the transcript length
int SEQCLUSTER::get_transcript_length(int sIndex){
	if( sIndex / 2  >= start_pos.size())
		return -1;
	else if (sIndex / 2 == start_pos.size() -1)
		return this -> sIndex_lookup.size() - this -> get_start_pos(sIndex);
	else
		return this -> get_start_pos(sIndex + 2) - get_start_pos(sIndex);
}

// for a given sequence id, return the transcript name
string SEQCLUSTER::get_transcript_name(int sIndex){
	int cIndex = sIndex_to_cIndex(sIndex);
	CINFO cluster_info = this -> get_cluster_info(cIndex);
	int offset = (sIndex  - cluster_info.c_start) / 2;
	return cluster_info.tnames[offset];
}			



// return the corresponding forward position
int SEQCLUSTER::convert_to_forward_pos(int pos, int sigmer_len){

	int sIndex = pos_to_sIndex(pos);	
/*	if(sIndex % 4 == 0)
		return pos;
*/	

	if(sIndex % 2 == 0)
		return pos;
	else
	{
		int t_start = this->get_start_pos(sIndex);
		int t_length = this-> get_transcript_length(sIndex);
		int t_end = t_start + t_length - 1;
	//	cout << endl << t_start << "\t" << t_length << "\t" << t_end << endl;
		return t_end - pos - sigmer_len + t_start;
	}

	return -1;
}

// return the concatenation of string
string SEQCLUSTER::getText(){
	return this -> text;
}

// return the cluster name based on given cluster index
// cIndex starts with 1
CINFO SEQCLUSTER::get_cluster_info(int cIndex){
	return this -> cluster_info[cIndex];
}
