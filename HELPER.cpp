#include "HELPER.hpp"

namespace helper{

	void echo(const string msg)
	{
		time_t now = time(0);
		struct tm *current = localtime(&now);
		printf("[%i-%i-%i %i:%i:%i] %s \n", (current->tm_year + 1900), (current->tm_mon + 1), current->tm_mday, current->tm_hour, current->tm_min, current->tm_sec, msg.c_str());
	}

	void upper(string &seq)
	{		
        	for(int i = 0; i < seq.size(); i ++)
        	{
                	seq[i] = toupper(seq[i]);
                	if(!(seq[i] == 'A' || seq[i] == 'C' || seq[i] == 'G' || seq[i] == 'T'))
                        	seq[i] = 'N';
        	}
	}	

	void reverse(string &seq)
	{
		reverse(begin(seq), end(seq));	

	}

	void complement(string &seq)
	{
		for(int i = 0; i < seq.size(); i++)
		{
			switch (seq[i])
			{
				case 'A':
					seq[i] = 'T';
					break;
				case 'C':
					seq[i] = 'G';
					break;
				case 'G':
					seq[i] = 'C';
					break;
				case 'T':
					seq[i] = 'A';
					break;
				default:
					break;
			}

		}
	}



	void reverse_complement(string &seq)
	{
        	reverse(begin(seq), end(seq));
        	for(int i = 0; i < seq.size(); i ++)
        	{
                	switch (seq[i])
                	{
                	case 'A':
                        	seq[i] = 'T';
                        	break;
                	case 'C':
                        	seq[i] = 'G';
                        	break;
                	case 'G':
                        	seq[i] = 'C';
                        	break;
                	case 'T':
                        	seq[i] = 'A';
                        	break;
                	default:
                        	break;
                	}		
        	}
	}

	void split(const string& s, char c, vector<string>& v) 
	{
   		string::size_type i = 0;
   		string::size_type j = s.find(c);

   		while (j != string::npos) {
      			v.push_back(s.substr(i, j-i));
      			i = ++j;
      			j = s.find(c, j);

      			if (j == string::npos)
         			v.push_back(s.substr(i, s.length()));
   		}
	}
}
	
