#include "SEQCLUSTER.hpp"
#include "HELPER.hpp"
#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <stack>

using namespace std;
using namespace helper;
using namespace sdsl;


void output_sigmer_v2(vector<int> positions, int length, ofstream& outfh)
{
	outfh << positions[0] << "\t";
	for(int i = 1; i < positions.size(); i++)
		outfh << positions[i] << ",";
	outfh << "\t" << length << endl;
}

void output_sigmer_v3(vector<int> positions, int prefix_length, int length, ofstream& outfh, SEQCLUSTER *scptr)
{
	int start_pos = positions[1];
	int end_pos = start_pos + length - 1;
	int sIndex = scptr -> pos_to_sIndex(start_pos);
	int seq_end_pos = scptr -> get_start_pos(sIndex) + (scptr -> get_transcript_length(sIndex))/2 - 1;
	int rev_seq_end_pos = scptr -> get_start_pos(sIndex) + scptr -> get_transcript_length(sIndex) - 1;
	
	if(start_pos + prefix_length <= seq_end_pos)
		end_pos = min(end_pos, seq_end_pos);
	else
		end_pos = min(end_pos, rev_seq_end_pos);

	if(end_pos == seq_end_pos || end_pos == rev_seq_end_pos)
		end_pos--;

	int adjust_length = end_pos - start_pos + 1;

	if(prefix_length < adjust_length){  // ?? < or <=
		outfh << positions[0] << "\t";
		for(int i = 1; i < positions.size(); i++)
			outfh << positions[i] << ",";	
		outfh << "\t" << prefix_length + 1 << "\t" << adjust_length<< endl;	

	//	cout << "Writing " << positions << "\t" << prefix_length +1 <<"\t" << adjust_length << endl;
	}
	//else
	//	cout << "No writing " << positions << "\t" << prefix_length + 1 << "\t" << adjust_length << endl;
}

void post_order_traversal(SEQCLUSTER *scptr, cst_sct3<> *cst_ptr, const char* outfile)
{
	ofstream outfh;
	outfh.open(outfile, ios_base::out);
	stack<vector<int>> cluster_stack;
	for (auto it=cst_ptr->begin(); it!=cst_ptr->end(); ++it)
	{
		auto  v = *it;
		vector <int> positions; // the first element is cluster_id, the rest is position

		// a leaf must be unique
		if(cst_ptr -> is_leaf(v))
		{
			// push information to stack
			int text_pos = cst_ptr -> csa[cst_ptr -> id(v)];
			int s_id = scptr -> pos_to_sIndex(text_pos);
			positions.push_back(scptr -> sIndex_to_cIndex(s_id));
			positions.push_back(text_pos);
			cluster_stack.push(positions);
		//	cout << "Child:" << text_pos << "(" << positions << ")" << endl;
		}
			
	
		// do something with the parent if it is the second time of visit (post-order)
		else if(it.visit() == 2 && cst_ptr->root() != v)
		{
			// if the parent string (prefix of a suffix) run across multiple sequences, no need to output the subtree
			bool output_subtree = true;
			int parent_length = cst_ptr -> depth(v);
			int parent_start_pos = cluster_stack.top()[1]; // get the position of the first child
			if( scptr -> pos_to_sIndex(parent_start_pos) != scptr -> pos_to_sIndex(parent_start_pos + parent_length))
				output_subtree = false;
		
		//	cout << "Parents:" << endl;
		
			// output children if its cluster is not -1
			// if any of the children is -1, its cluster is automatically set to be -1
			for(int i = cst_ptr->children(v).size()-1; i >= 0; i--)
			{
				auto child = cst_ptr -> children(v)[i];
				int child_length = cst_ptr -> depth(child);
				vector <int> child_positions = cluster_stack.top();	
				cluster_stack.pop();

				// output children only if 1) output_subtree is true (i.e. the prefix doesn't run across two transcripts)
				// 			2) the children is unique (i.e. positions[0] != -1)
				if( output_subtree && child_positions[0] != -1)
				{
					int child_start_pos = child_positions[1];
					output_sigmer_v3(child_positions, parent_length, child_length, outfh, scptr);
				}

				// update uniqueness
				if(child_positions[0] == -1 || (positions.size() > 0 && positions[0] == -1) || (positions.size() > 0 && child_positions[0] != positions[0]))
				{
					positions.resize(2);
					positions[0] = -1;
					positions[1] = child_positions[1];
				}
				else
				{
					if(positions.size() == 0)
						positions.insert(positions.end(), child_positions.begin(), child_positions.end());
					else
						positions.insert(positions.end(), child_positions.begin()+1, child_positions.end());
				}
			}	
	//		cout << "Insert " << positions << endl;
			cluster_stack.push(positions);

		}
	}
	outfh.close();	

}


int main(int argc, char* argv[])
{
        if(argc < 3)
        {
                cerr << "Usage: " << argv[0] << " INFILE OUTFILE" << endl;
                return 1;
        }

        const char *infile = argv[1];
        const char *outfile = argv[2];
	const char *treefile = argv[3];
	
	echo("Reading Sequence Data");
	SEQCLUSTER *sc = new SEQCLUSTER(infile);

	echo("Building Suffixtree");
	cst_sct3 <> cst;
	construct_im(cst, sc -> getText(), 1);

	cout << "total nodes : " << cst.nodes() << endl;
     	cout << "leaf ndoes : " << cst.csa.size() << endl;


	echo("Traversing suffixtree");
	post_order_traversal(sc, &cst, outfile);

	echo("Done");

}

