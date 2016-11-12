#include "SEQCLUSTER.hpp"
#include "HELPER.hpp"
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
using namespace helper;


void printMatrix(vector <vector<double>> v)
{
        for(int i = 0; i < v.size(); i++)
        {
                for(int j = 0; j < v[0].size(); j++)
                        cout << v[i][j] << ",";
                cout << endl;
        }
}


void E_step(vector <vector <double>> & Y_matrix, vector<vector<double>> &Z_matrix, vector<double> &transcript_abundance, SEQCLUSTER *sc, int readlen, int cid)
{
	int cluster_start = sc -> get_cluster_info(cid).c_start / 2;

	// Z_ik = yik * p_k / sum (yik * p_k)
	for(int i = 0; i < Z_matrix.size(); i++){
		// calculating rowsum
		double rowsum = 0;
		for(int k = 0; k < Z_matrix[i].size(); k++)
		{
			double e_length = (sc -> get_transcript_length ((k + cluster_start)*2)) / 2 - readlen;

			if(e_length > 0)
				rowsum += Y_matrix[i][k] * transcript_abundance[k] / e_length;
		}

		// update Z matrix
		if(rowsum > 0){
			for (int k = 0; k < Z_matrix[i].size(); k++)
			{
                	        double e_length = (sc -> get_transcript_length ((k + cluster_start)*2)) / 2 - readlen;
				if(e_length > 0)
					Z_matrix[i][k] = (Y_matrix[i][k] * transcript_abundance[k] / e_length) / rowsum;	

			}	
		}
	}

}

void printAbundance(vector <double> abundance)
{
        for(int i = 0; i < abundance.size(); i++)
                cout << abundance[i] << "\t";
        cout << endl;
}

double M_step(vector <vector<double>> &Z_matrix, vector<double> &transcript_abundance, const vector<double> &C_vector, const double &total_reads)
{
	double difference = 0;

	for(int k = 0; k < transcript_abundance.size(); k++)
	{
		double col_sum = 0;
		for(int i = 0; i < Z_matrix.size(); i++){
			col_sum += Z_matrix[i][k] * (double) C_vector[i];	
		}
	
		difference += fabs(transcript_abundance[k] - (col_sum / total_reads));
		transcript_abundance[k] = col_sum / total_reads;
	}
	return difference;
}


void initialize_ymatrix(const string yfile, vector <vector <vector<double>>> &ymatrix, vector <vector<double>> &cvector, vector <double> &total_reads, SEQCLUSTER* sc){

        //read in ymatrix, update cvector and total_reads
        ifstream infile(yfile, ifstream::in);
	if(infile.is_open()){
		string line;
		while(getline(infile, line)){
			vector<string> file_data;
			split(line, '\t', file_data); // file_data = [cid, profile, count]

			int cid = stoi(file_data[0]);
			string profile = file_data[1];
			double profile_count = stod(file_data[2]);

			vector<double> tmp_yrow(profile.size(), 0);
			for(int i = 0; i < profile.size(); i++){
				if(profile[i] == '1'){
					tmp_yrow[i] = 1;
				}
			}

			ymatrix[cid].push_back(tmp_yrow);
			total_reads[cid] += profile_count;
			cvector[cid].push_back(profile_count);		
		}
	}			

	infile.close();
}

void export_abundance(const string outfile, SEQCLUSTER* sc, const vector <vector<double>> &transcript_abundance, int readlen, vector<double> total_reads ){

	ofstream outfh;
        outfh.open(outfile, ios_base::out);
	outfh << "Name\tLength\tEffectiveLength\tNumReads\tRPKM\tTPM" << endl;

	double total_tpm = 0;
	double readSum = 0;

	for (int c = 0; c < transcript_abundance.size(); c++){ // for each cluster
		int cluster_start = sc -> get_cluster_info(c).c_start / 2;
	
		for(int t = 0; t < transcript_abundance[c].size(); t++){ 	// for each transcript

			double e_length = (sc -> get_transcript_length ( (t + cluster_start) * 2)) / 2 - readlen;
			if(e_length >0){
				double readnums = floor(transcript_abundance[c][t] * total_reads[c] + 0.5);
				total_tpm += (double) readnums / e_length;
			}
		}
		readSum += total_reads[c];
	}	

        for (int c = 0; c < transcript_abundance.size(); c++){  // for each cluster
                int cluster_start = sc -> get_cluster_info(c).c_start / 2;
                for(int t = 0; t < transcript_abundance[c].size(); t++){   //for each transcript
			int t_length = (sc -> get_transcript_length( (t + cluster_start)*2)) / 2;
                        double e_length = t_length - readlen;
                        double readnums = floor(transcript_abundance[c][t] * total_reads[c]);

			double tpm = 0;
			double rpkm = 0;
			if(e_length > 0){
				tpm = ((double) readnums / e_length) / total_tpm * 1000000;
				rpkm = (readnums * 1000000000) / (e_length * readSum) ;
			}
			string t_name = sc -> get_transcript_name( (t + cluster_start)*2);
			outfh << t_name << "\t" << t_length - 1 << "\t" << e_length << "\t" << readnums << "\t" << rpkm << "\t" << tpm << endl;
                }
        }

	outfh.close();
}


int main(int argc, char* argv[]){

       if(argc < 4){
                cerr << "Usage: " << argv[0] << " SEQFILE YMATRIX FRAGLEN OUTFILE" << endl;
                return 1;
        }

        const char *seqfile = argv[1];
	const char *yfile = argv[2];
	const int fraglen = stoi(argv[3]);
	const char *outfile = argv[4];

	SEQCLUSTER *sc = new SEQCLUSTER(seqfile);

	vector <vector <double>> transcript_abundance (sc -> get_num_cluster(), vector<double>()); // k* 1 (for k transcripts)
	vector <vector < vector <double>>> Y_matrix (sc-> get_num_cluster(), vector< vector<double>>()); 
	vector <vector <double>> C_vector (sc -> get_num_cluster(), vector<double>());
	vector <vector < vector <double>>> Z_matrix (sc -> get_num_cluster(), vector< vector<double>>());
	vector <double> total_reads (sc-> get_num_cluster(), 0);

	// initialize Y_matrix
	echo("Loading Observation");
	initialize_ymatrix(yfile, Y_matrix, C_vector, total_reads, sc);


	// initialize transcript abundance and Z_matrix
	echo("Initializaing Transcript Abundance and Hidden Variable");
	for (int c = 0; c < Y_matrix.size(); c++){

		if(Y_matrix[c].size() > 0){
			Z_matrix[c] = vector<vector<double>>( Y_matrix[c].size(), vector<double>(Y_matrix[c][0].size(), 0));
			transcript_abundance[c] = vector<double> (Y_matrix[c][0].size(), 1.0 / (double) Y_matrix[c][0].size());
		}
		else
			transcript_abundance[c] = vector<double> (sc-> get_cluster_info(c).tnames.size(), 0.0);
	}	

	echo("Running EM");
	for(int i = 0; i < sc -> get_num_cluster(); i++){
		int it = 0;
		double diff = 0;
		if(Y_matrix[i].size() > 0){
			diff = 100;
			while(diff > 0.000001 && it < 10000)
			{
				E_step(Y_matrix[i], Z_matrix[i], transcript_abundance[i], sc, fraglen, i);
				diff = M_step(Z_matrix[i], transcript_abundance[i], C_vector[i], total_reads[i]);
				it ++;
			}
		}

	//	cout << "Cluster " << i << ": converge at step " << it << ": " << diff << endl;
	}


	export_abundance(outfile, sc, transcript_abundance, fraglen, total_reads);
	echo("Done");
}
