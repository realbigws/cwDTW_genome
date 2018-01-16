#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "proc/proc.h"
#include "util/exception.h"
#include "5mer/5mer_index.h"
#include "wavelib.h"
#include <malloc.h>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace g::proc;

//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}
//---------------------- utility ------//over



//------------ read DNA --------//
void Read_DNA(string &input,string &out,string &first_line)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(input.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",input.c_str());
		exit(-1);
	}
	//skip first
	if(!getline(fin,first_line,'\n'))
	{
		fprintf(stderr,"file %s format bad!\n",input.c_str());
		exit(-1);
	}
	//read content
	out="";
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		out+=buf;
	}
}

//------------ reverse ---------//
char DNA_Reverse(char c,int &bad)
{
	bad=0;
	switch(c)
	{
		case 'A':return 'T';
		case 'T':return 'A';
		case 'C':return 'G';
		case 'G':return 'C';
		case 'N':
		{
			bad=1;
//			fprintf(stderr,"BAD N CODE !!! \r");
			return 'N';
		}
		default:
		{
			bad=1;
//			fprintf(stderr,"BAD CODE !!! %c \n",c);
			return 'N';
		}
	}
}

//------------ reverse complement ---------//
int Reverse_Complement(string &in,string &out)
{
	int i;
	out=in;
	int len=(int)in.length();
	int BAD_CODE=0;
	int fin_bad=0;
	for(i=len-1;i>=0;i--)
	{
		out[len-1-i]=DNA_Reverse(in[i],BAD_CODE);
		if(BAD_CODE==1)fin_bad=1;
	}
	return fin_bad;
}

//-------- calcualte N ratio --------//
double Calc_N_Ratio(string &in)
{
	int i;
	int len=(int)in.length();
	int count=0;
	int BAD_CODE=0;
	for(i=0;i<len;i++)
	{
		if(DNA_Reverse(in[i],BAD_CODE)=='N')count++;
	}
	return 1.0*count/len;
}



//================================ local extension ============================//

//--------- local match -------//
int local_match(string &seq1,string &seq2,long ii,long jj,long wind)
{
	for(long i=-wind;i<=wind;i++)
	{
		if(seq1[ii+i]!=seq2[jj+i])return 0;
	}
	return 1;
}

//--------- local mis-match -------//
long local_mismatch(string &seq1,string &seq2,long ii,long jj,long wind)
{
	long mismatch=0;
	for(long i=-wind;i<=wind;i++)
	{
		if(seq1[ii+i]!=seq2[jj+i])mismatch++;
	}
	return mismatch;
}
long local_mismatch_left(string &seq1,string &seq2,long ii,long jj,long wind)
{
	long mismatch=0;
	for(long i=-wind;i<=0;i++)
	{
		if(seq1[ii+i]!=seq2[jj+i])mismatch++;
	}
	return mismatch;
}
long local_mismatch_right(string &seq1,string &seq2,long ii,long jj,long wind)
{
	long mismatch=0;
	for(long i=0;i<=wind;i++)
	{
		if(seq1[ii+i]!=seq2[jj+i])mismatch++;
	}
	return mismatch;
}

//--------- local search --------//
int Local_Search(string &seq1,string &seq2,long ini_ii,long ini_jj,
	long &ii,long &jj,long left_range, long right_range, long wind)
{
	long i,j;
	long len1=seq1.length();
	long len2=seq2.length();
	ii=-1;
	jj=-1;
	for(i=ini_ii-left_range;i<=ini_ii+right_range;i++)
	{
		if(i<wind)continue;
		if(i>len1-1-wind)break;
		for(j=ini_jj-left_range;j<=ini_jj+right_range;j++)
		{
			if(j<wind)continue;
			if(j>len2-1-wind)break;
			int match=local_match(seq1,seq2,i,j,wind);
			if(match==1)
			{
				ii=i;
				jj=j;
				return 1; //found
			}
		}
	}
	//not found
	return 0;
}



//--------------- gapless extension ----------------//
void Gapless_Extension(string &seq1,string &seq2,
	vector<pair<long, double> > &sigpeak1, 
	vector<pair<long, double> > &sigpeak2, 
	vector<pair<long,long> > &align,
	vector <long> &ali1,vector <long> &ali2)
{
	//init
	long len1=seq1.length();
	long len2=seq2.length();
	ali1.assign(len1,-1);
	ali2.assign(len2,-1);
	//maxcur
	long cur1=0;
	long cur2=0;
	//gapless extension
	for(long i=1;i<align.size()-1;i++)
	{
		//---- current position
		long pos1=sigpeak1[align[i].first].first;
		long pos2=sigpeak2[align[i].second].first;
		if(ali1[pos1]!=-1 || ali2[pos2]!=-1)continue;
		if(pos1<cur1 || pos2<cur2)continue;

		//----- local search
//		if(seq1[pos1]!=seq2[pos2])
		{
/*
			//---> left boundary
			long left_pos1=sigpeak1[align[i-1].first].first;
			long left_pos2=sigpeak2[align[i-1].second].first;
			long left_len1=pos1-left_pos1-1;
			long left_len2=pos2-left_pos2-1;
			long left_len=left_len1<left_len2?left_len1:left_len2;
	
			//---> right boundary
			long right_pos1=sigpeak1[align[i+1].first].first;
			long right_pos2=sigpeak1[align[i+1].second].first;
			long right_len1=right_pos1-pos1-1;
			long right_len2=right_pos2-pos2-1;
			long right_len=right_len1<right_len2?right_len1:right_len2;
*/

			//---> local search
			long pos1_,pos2_;
			long range=10;
			long wind=3;
			int retv=Local_Search(seq1,seq2,pos1,pos2,pos1_,pos2_,range,range,wind);
			if(retv==1)
			{
				pos1=pos1_;
				pos2=pos2_;
			}
			else
			{
				continue;
			}
		}
		//final check
		if(ali1[pos1]!=-1 || ali2[pos2]!=-1)continue;
		if(pos1<cur1 || pos2<cur2)continue;
			
		ali1[pos1]=pos2;
		ali2[pos2]=pos1;
		//---- left extension
		//-> get boundary
		long wind=3;
		long misnum=3;

		//-> left extension
		long mismatches_left=0;
		long maxmismatch_left=(wind+1)*misnum;
		for(long j=1;j<=align.size();j++)
		{
			if(pos1-j<0 || pos2-j<0)break; 
			if(ali1[pos1-j]!=-1 || ali2[pos2-j]!=-1)break;

			//--| check mismatch
			long mis=local_mismatch_left(seq1,seq2,pos1-j,pos2-j,wind);
			if(mis==0)
			{
				ali1[pos1-j]=pos2-j;
				ali2[pos2-j]=pos1-j;
			}
			else
			{
				if(mis>1)break;
				if(mis==1)
				{
					mismatches_left++;
					if(mismatches_left>maxmismatch_left)break;
					ali1[pos1-j]=pos2-j;
					ali2[pos2-j]=pos1-j;
				}
			}
		}
		//---- right extension

		//-> left extension
		long mismatches_right=0;
		long maxmismatch_right=(wind+1)*misnum;
		for(long j=1;j<=align.size();j++)
		{
			if(pos1+1>len1-1 || pos2+j>len2-1)break;
			if(ali1[pos1+j]!=-1 || ali2[pos2+j]!=-1)break;

			//--| check mismatch
			long mis=local_mismatch_right(seq1,seq2,pos1+j,pos2+j,wind);
			if(mis==0)
			{
				ali1[pos1+j]=pos2+j;
				ali2[pos2+j]=pos1+j;
				cur1=pos1+j;
				cur2=pos2+j;
			}
			else
			{
				if(mis>1)break;
				if(mis==1)
				{
					mismatches_right++;
					if(mismatches_right>maxmismatch_right)break;
					ali1[pos1+j]=pos2+j;
					ali2[pos2+j]=pos1+j;
					cur1=pos1+j;
					cur2=pos2+j;
				}
			}
		}
	}
}





//--------- pore_models: from genome to expected signal ----------//
void Genomes2SignalSequence(const std::vector<char> &genomes, 
	std::vector<int>& index, std::vector<double>& signals, 
	int scale, int FIVE_or_SIX, int ZSCO_or_NOT)
{
	long bound;
	if(FIVE_or_SIX==0) //-> 5mer model
	{
		g::Mer2Signal::Genome2Index_5mer(genomes, index);
		bound = genomes.size()-5;//genomes.size()%5;
		signals.assign(bound*scale,0);
		long cur=0;
		for(long i = 0; i < bound; i++){
			double sigval;
			if(index[i]<0)sigval=0;
			else{
				sigval = g::Mer2Signal::AvgSignalAt_5mer(index[i]);
				if(ZSCO_or_NOT==1) //-> transfer to Zsco
				{
					sigval = (sigval-90.208351)/12.832660;
				}
				else               //-> use original int value
				{
					sigval = (int)(5.7*sigval+14);
				}
			}
			for(int c = scale; c--;){
				signals[cur]=sigval;
				cur++;
			}
		}
	}
	else               //-> 6mer model
	{
		g::Mer2Signal::Genome2Index_6mer(genomes, index);
		bound = genomes.size()-6;//genomes.size()%5;
		signals.assign(bound*scale,0);
		long cur=0;
		for(long i = 0; i < bound; i++){
			double sigval;
			if(index[i]<0)sigval=0;
			else{
				sigval = g::Mer2Signal::AvgSignalAt_6mer(index[i]);
				if(ZSCO_or_NOT==1) //-> transfer to Zsco
				{
					sigval = (sigval-90.208199)/12.868652;
				}
				else               //-> use original int value
				{
					sigval = (int)(5.7*sigval+14);
				}
			}
			for(int c = scale; c--;){
				signals[cur]=sigval;
				cur++;
			}
		}
	}

	//---- tail five_mer ------//
	for(long i = bound; i < genomes.size(); i++){
		for(int c = scale; c--;){
			if(ZSCO_or_NOT==1) //-> transfer to Zsco
			{
				signals.push_back(0);
			}
			else
			{
				signals.push_back(100);
			}
		}
	}

}

//--------------- continuous wavelet transform (CWT) analysis -----------------//
/** @scale0: level0 pyramind scale;  @dscale: scale_i = scale0*(2^{i*dsacle} ); @npyr: total number of pyramind*/
void CWTAnalysis(const std::vector<double>& raw, std::vector<std::vector<double> >& output, 
	double scale0, double dscale, int npyr)
{
	const double* sigs = &raw[0];		//sst_nino3.dat
	cwt_object wt;

	size_t N = raw.size();
	double dt = 1;//2;		//sample rate	>  maybe we should use 2?

	wt = cwt_init("dog", 2.0, N, dt, npyr);	//"morlet", "dog", "paul"
	setCWTScales(wt, scale0, dscale, "pow", 2.0);
	cwt(wt, sigs);

	output.resize(npyr);
	for(size_t k = npyr; k--;){
		int idx = npyr-k-1;
		
		output[idx].resize(raw.size());
		size_t offset = k*raw.size();
		for(size_t i = 0; i < output[idx].size(); i++){
			output[idx][i] = wt->output[i+offset].re;
		}
	}

	cwt_free(wt);
}

//----------------- boundary generation for constrained dynamic time warping (cDTW) -------------//
void BoundGeneration(std::vector<std::pair<long,long> >& cosali, 
	long neib, std::vector<std::pair<long,long> >& bound, int mode,
	int RENMIN_or_SHENG=0)
{

//if(mode!=-1) //-> mode = -1 means Renmin mode
vector<pair<long,long> > cosali_=cosali;
if(mode==1) //-> use partial-diagonol alignment
{
	//-------- generate partial-diagonol alignment -------//
	vector<pair<long,long> > cosali2;
	cosali2.push_back(cosali[0]);
	for(long i = 1; i < cosali.size(); i++){
		if(cosali[i].first != cosali2[cosali2.size()-1].first){
			cosali2.push_back(cosali[i]);
		}
		else{
			cosali2[cosali2.size()-1].second = cosali[i].second;
		}
	}

	cosali_.clear();
	for(long i = 1; i < cosali2.size(); i++)
	{
		long pre_idx = cosali2[i-1].first, cur_idx = cosali2[i].first;
		long pre_anchor = cosali2[i-1].second, cur_anchor = cosali2[i].second;
		double anchor_diff = 1.0*(cur_anchor-pre_anchor)/(double)(cur_idx-pre_idx);
		for(long k = pre_idx, count = 0; k < cur_idx; k++, count++)
		{
			long mid = pre_anchor+(long)(count*anchor_diff);  //assign point relationship
			if(mid<pre_anchor)mid=pre_anchor;
			if(mid>cur_anchor)mid=cur_anchor;
			cosali_.push_back(make_pair(k,mid));
		}
	}
	cosali_.push_back(cosali2[cosali2.size()-1]);
}

	//----> output pre-bound alignnment -------
//	{
//		for(int i = 0; i < cosali_.size(); i++)
//		{
//			printf("%d -> %d  %d \n",i,cosali_[i].first,cosali_[i].second);
//		}
//	}


	//---------- use block bound ------------//start
	//-> get signal length
	long moln1=cosali_[cosali_.size()-1].first+1;
	long moln2=cosali_[cosali_.size()-1].second+1;
	//-> renmin align to sheng style
	std::vector<std::pair<long,long> > align_sheng;
	Renmin_To_Sheng_align(moln1,moln2,cosali_,align_sheng);
	//-> get bound in sheng style
	std::vector<std::pair<long,long> > bound_sheng;
	From_Align_Get_Bound(moln1,moln2,align_sheng,bound_sheng,neib);
	//-> transfer bound to renmin style
	Sheng_To_Renmin_bound(moln1,moln2,bound_sheng,bound);
	//----- maybe useful -----//
	bound[0].first = 0;
	bound[bound.size()-1].first = bound[bound.size()-2].first;
	bound[bound.size()-1].second = cosali[cosali.size()-1].second;


	if(RENMIN_or_SHENG==1) //-> use Sheng's bound definition
	{
		bound=bound_sheng;
	}

	//----> output post-bound alignnment -------
//	{
//		for(int i = 0; i < bound.size(); i++)
//		{
//			printf("%d -> %d  %d \n",i,bound[i].first,bound[i].second);
//		}
//		exit(-1);
//	}

	return;
	//---------- use block bound ------------//over


/*
	int radius=neib1;
	std::vector<std::pair<int,int> > cosali2;
	
	cosali2.push_back(cosali[0]);
	for(int i = 1; i < cosali.size(); i++){
		if(cosali[i].first != cosali2[cosali2.size()-1].first){
			cosali2.push_back(cosali[i]);
		}
		else{
			cosali2[cosali2.size()-1].second = cosali[i].second;
		}
	}
	
	cosali = cosali2;
	
	bound.resize(cosali[cosali.size()-1].first+1);
	bound.assign(cosali[cosali.size()-1].first+1, std::make_pair(-1,-1));
	
	for(int i = 1; i < cosali.size(); i++){
		int pre_idx = cosali[i-1].first, cur_idx = cosali[i].first;
		int pre_anchor = cosali[i-1].second, cur_anchor = cosali[i].second;
		
		double anchor_diff;
		
		anchor_diff = (cur_anchor-pre_anchor)/(cur_idx-pre_idx);
		
		int neighbor = int(anchor_diff+0.5)+radius;
		
		for(int i = pre_idx, count = 0; i < cur_idx; i++, count++){
			int mid = pre_anchor+std::round(count*anchor_diff);		//assign point relationship
			bound[i].first = mid-neighbor < 0 ? 0 : mid-neighbor;
			bound[i].second = mid+neighbor > cosali[cosali.size()-1].second ? cosali[cosali.size()-1].second : mid+neighbor;
		}
	}
	
// 	for(int i = 0; i < bound.size(); i++){
// 		std::cout<<bound[i].first<<" "<<bound[i].second<<std::endl;
// 	}
	
	bound[0].first = 0;
	bound[bound.size()-1].first = bound[bound.size()-2].first;
	bound[bound.size()-1].second = cosali[cosali.size()-1].second;
*/
}

//====================== continous wavelet dynamic time warping (cwDTW) ========================//
void MultiLevel_WaveletDTW(std::vector<double>& in1, std::vector<double>& in2,
	std::vector<std::vector<double> >& sig1, std::vector<std::vector<double> >& sig2, 
	std::vector<std::pair<long,long> >& alignment, long radius, int test, int mode, 
	double* totaldiff = 0)
{
	double tdiff;
	std::vector<std::pair<long, double> > sig1peaks, sig2peaks;
	double length1 = sig1[0].size(), length2 = sig2[0].size();

	int tot_size=sig1.size();
	for(long k = 0; k < tot_size; k++)
	{

//printf("layer k = %d \n",k);

		//------ peakpick CWT signal -------------//
		g::proc::PeakPick(sig1[k], sig1peaks);
		g::proc::PeakPick(sig2[k], sig2peaks);
//		std::cout<<sig1peaks.size()<<"\t"<<sig2peaks.size()<<std::endl;
		std::vector<double> peak1(sig1peaks.size());
		std::vector<double> peak2(sig2peaks.size());


		//-------- use peak picking result ---------//
		for(long i = sig1peaks.size(); i--;){
			peak1[i] = sig1peaks[i].second;
		}
		for(long i = sig2peaks.size(); i--;){
			peak2[i] = sig2peaks[i].second;
		}

		//----- apply DTW or cDTW dependent on k-th level -------//
		if(k == 0){
			tdiff = g::proc::DynamicTimeWarping(peak1, peak2, alignment);
		}
		else{
			//----- ReMapIndex_partI (map ground level upon k-th level) -----//
			long c = 0;
			for(long i = 0; i < alignment.size(); i++){
				while(sig1peaks[c].first < alignment[i].first){
					c++;
				}
				alignment[i].first = c;	
			}
			
			long d = 0;
			for(long i = 0; i < alignment.size(); i++){
				while(sig2peaks[d].first < alignment[i].second){
					d++;
				}
				alignment[i].second = d;	
			}
			//----- cDWT (constrained DWT) -------//
			std::vector<std::pair<long,long> > bound;
			int neib=radius*pow(2,tot_size-k); //adaptive radius
			BoundGeneration(alignment, neib, bound, mode);
			tdiff = g::proc::BoundDynamicTimeWarping(peak1, peak2, bound, alignment);
		}
		//----- ReMapIndex_partII (map k-th level back to ground level) -----//
		for(long i = alignment.size(); i--;){
			alignment[i].first = sig1peaks[alignment[i].first].first;
			alignment[i].second = sig2peaks[alignment[i].second].first;
		}

	}
	
	if(totaldiff){
		*totaldiff = tdiff;
	}
}


//================================== DNA_DynaProg_NanoPore_assess_FixFast =============================//

//=================== NanoPore DynaProg ===============//
//-> ATCG to digit
//-> see below link for base pairing:
//       http://www.biology-pages.info/B/BasePairing.html
//-> see below link for molecular weight:
//       http://www.genomics.agilent.com/files/Mobio/Nucleic%20Acids_Sizes_and_Molecular_Weights_2pgs.pdf
int DNA_To_Int(char c)
{
	switch(c)
	{
		case 'A': return -2;
		case 'T': return -1;
		case 'C': return  1;
		case 'G': return  2;
		return 0;
	}
}

//-------- DNA calc ------//
int DNA_Calc(char a,char b)
{
	int mat=1;
	int mis=-2;
	int s1=DNA_To_Int(a);
	int s2=DNA_To_Int(b);
	if(s1==0 || s2==0)return mis;
	if(a==b)return mat;
	else return mis;
}

//---- part 2. generate alignment from a given bound ------//
int Advance_Align_Dyna_Prog_Double_Bound(long n1,long n2,const vector <vector<double> > &score,
	double GAP_OPEN1,double GAP_EXT1,double GAP_OPEN2,double GAP_EXT2,
	double GAP_HEAD1,double GAP_TAIL1,double GAP_HEAD2,double GAP_TAIL2,
	vector<pair<long,long> > &bound, vector<pair<long,long> > & alignment,double &ali_sco)
{
	long i,j,k;
	//input
	long m = n1 + 1;  // +1 to account for the extra row,col in
	long n = n2 + 1;  // the DP matrices corresponding to gaps
	//const value
	const int _H_  = 0;
	const int _S_  = 1;
	const int _V_  = 2;

	//create D and M
	vector <vector <vector <int> > > D;      // the path (directions) matrix
	vector <vector <vector <double> > > M;   // the current scores (values) matrix
	D.resize(3);
	M.resize(3);
	//resize(m,n)
	long wsstart,wsend;
	long prev_start,prev_end;
	for (k = 0; k < 3; k++) 
	{
		D[k].resize(m);
		M[k].resize(m);
		//-> init
		for(i=0; i<m; i++)
		{
			//current bound
			wsstart=bound[i].first;
			wsend=bound[i].second;
			long bound_len=wsend-wsstart+1;
			D[k][i].resize(bound_len,-1);
			M[k][i].resize(bound_len,-1);
		}
	}
	//init()
	double WS_MIN=-1000000;
	D[_S_][0][0] = -1;
	D[_H_][0][0] = -1;
	D[_V_][0][0] = -1;
	M[_S_][0][0] = 0;
	M[_H_][0][0] = WS_MIN;
	M[_V_][0][0] = WS_MIN;
	//first line
	i=0;
	wsstart=bound[i].first;
	wsend=bound[i].second;
	for (j = 1; j <= wsend; j++) 
	{
		D[_S_][i][j-wsstart] = _H_;
		D[_H_][i][j-wsstart] = _H_;
		D[_V_][i][j-wsstart] = _H_;
		M[_S_][i][j-wsstart] = WS_MIN;
		M[_H_][i][j-wsstart] = j*GAP_HEAD2; //-(Params::GAP_OPEN + (j-1)*Params::GAP_EXT);
		M[_V_][i][j-wsstart] = WS_MIN;
	}
	prev_start=wsstart;
	prev_end=wsend;
	//fill(firstSeq, secondSeq, distFunc);
	double gap_open;
	double gap_ext;
	double v1,v2,v3;
	double dist;
	for (i = 1; i < m; i++) 
	{
		//current bound
		wsstart=bound[i].first;
		wsend=bound[i].second;
		//process
		for (j = wsstart; j <= wsend; j++) 
		{
			//------ 1. condition upper -----//
			if(j==n-1)
			{
				gap_open=GAP_TAIL1;
				gap_ext=GAP_TAIL1;
			}
			else
			{
				gap_open=GAP_OPEN1;
				gap_ext=GAP_EXT1;
			}
			//-- boundary --//(start)
			if(j>prev_end)
			{
				M[_V_][i][j-wsstart] = WS_MIN;
				D[_V_][i][j-wsstart] = -1;
			}
			else
			{
				if(j==0)
				{
					gap_ext=GAP_HEAD2;
					gap_open=GAP_HEAD2;
				}
				long curwsstart=bound[i-1].first;
				v1 = M[_V_][(i-1)][j-curwsstart] + gap_ext;
				v2 = M[_S_][(i-1)][j-curwsstart] + gap_open;
				v3 = M[_H_][(i-1)][j-curwsstart] + gap_open;
				M[_V_][i][j-wsstart] = std::max(v1, std::max(v2, v3));
				if (M[_V_][i][j-wsstart] == v1) D[_V_][i][j-wsstart] = _V_;
				else if(M[_V_][i][j-wsstart] == v2) D[_V_][i][j-wsstart] = _S_;
				else D[_V_][i][j-wsstart] = _H_;
			}
			//-- boundary --//(end)

			//------ 2. condition left -----//
			if(i==m-1)
			{
				gap_open=GAP_TAIL2;
				gap_ext=GAP_TAIL2;
			}
			else
			{
				gap_open=GAP_OPEN2;
				gap_ext=GAP_EXT2;
			}
			//-- boundary --//(start)
			if(j==wsstart)
			{
				M[_H_][i][j-wsstart] = WS_MIN;
				D[_H_][i][j-wsstart] = -1;
			}
			else
			{
				v1 = M[_H_][i][j-1-wsstart] + gap_ext;
				v2 = M[_S_][i][j-1-wsstart] + gap_open;
				v3 = M[_V_][i][j-1-wsstart] + gap_open;
				M[_H_][i][j-wsstart] = std::max(v1, std::max(v2, v3));
				if (M[_H_][i][j-wsstart] == v1) D[_H_][i][j-wsstart] = _H_;
				else if(M[_H_][i][j-wsstart] == v2) D[_H_][i][j-wsstart] = _S_;
				else D[_H_][i][j-wsstart] = _V_;
			}
			//-- boundary --//(end)

			//------ 3. condition diag -----//
			//-- boundary --//(start)
			if(j==prev_start || j>prev_end+1)
			{
				M[_S_][i][j-wsstart] = WS_MIN;
				D[_S_][i][j-wsstart] = -1;
			}
			else
			{
				long curwsstart=bound[i-1].first;
				long start=wsstart-1<0?0:wsstart-1;
				dist = score[i-1][j-1-start];
				v1 = M[_V_][(i-1)][j-1-curwsstart] + dist;
				v2 = M[_H_][(i-1)][j-1-curwsstart] + dist;
				v3 = M[_S_][(i-1)][j-1-curwsstart] + dist;
				M[_S_][i][j-wsstart] = std::max(v1, std::max(v2, v3));
				if (M[_S_][i][j-wsstart] == v3) D[_S_][i][j-wsstart] = _S_;
				else if (M[_S_][i][j-wsstart] == v1) D[_S_][i][j-wsstart] = _V_;
				else D[_S_][i][j-wsstart] = _H_;
			}
			//-- boundary --//(end)
		}
		//next bound
		prev_start=wsstart;
		prev_end=wsend;
	}

	//build(ali, firstSeq, secondSeq, distFunc);
	i = m-1;
	wsstart=bound[i].first;
	j = n-1;
	v1=M[_V_][i][j-wsstart];
	v2=M[_H_][i][j-wsstart];
	v3=M[_S_][i][j-wsstart];
	double maximal = std::max(v1, std::max(v2, v3));
	k = -1;
	if(v3==maximal)k = _S_;
	else if(v2==maximal)k = _H_;
	else k = _V_;
	//trace_back
	alignment.clear();
	long count = 0;
	long matches = 0;
	long cur_case=k;
	long pre_case;
	for(;;)
	{
		if(i==0||j==0)break;
		wsstart=bound[i].first;
		wsend=bound[i].second;
		pre_case=D[cur_case][i][j-wsstart];
		switch (cur_case)
		{
			case _S_:
				alignment.push_back(pair<long,long>(i,j)); 
				i--;
				j--;
				++matches;
				break;
			case _V_:
				alignment.push_back(pair<long,long>(i,-j)); 
				i--;
				break;
			case _H_:
				alignment.push_back(pair<long,long>(-i,j)); 
				j--;
				break;
			default:
				cout << "ERROR!! -> advance_global: invalid direction D[" << k << "](" << i << ", " << j << ") = " 
				<< D[k][i][j-wsstart] << endl;
				exit(-1);
		}
		cur_case=pre_case;
		count++;
	}
	while (j> 0) alignment.push_back(pair<long,long>(-i,j)),j--;
	while (i> 0) alignment.push_back(pair<long,long>(i,0)), i--;
	reverse(alignment.begin(), alignment.end());
	ali_sco=maximal;
	return matches;
}


//----------------- align two DNA sequences ---------------------//
//--> generate initial alignment
void Gen_Init_Align(long n1,long n2,vector<pair<long,long> > &align)
{
	align.clear();
	long i;
	double ratio;
	vector <long> ali1(n1,-1);
	if(n1<=n2)
	{
		ratio=1.0*n2/n1;
		for(i=0;i<n1;i++)
		{
			long pos2=(long)(1.0*i*ratio);
			ali1[i]=pos2;
		}
	}
	else
	{
		ratio=1.0*n1/n2;
		for(i=0;i<n2;i++)
		{
			long pos1=(long)(1.0*i*ratio);
			ali1[pos1]=i;
		}
	}
	Ali1_To_AliPair(n1,n2,ali1,align);
}

//--> apply dynamic programing
double process_oriami_record_simp(string &seq_,string &ami_,
	vector<pair<long,long> > &WWW_alignment,long &matchs,
	long neib_min,double neib_ratio)
{
	long i,j;
	long n1,n2;
	//--[0] initial alignment
	n1=seq_.length();
	n2=ami_.length();
	long neib=neib_ratio*(abs(n1-n2)+neib_min);  //-> 2*(abs(n1-n2)+5)
	vector<pair<long,long> > init_ali;
	Gen_Init_Align(n1,n2,init_ali);
	//--[1] bounded alignment	
	vector<pair<long,long> > bound;
	From_Align_Get_Bound(n1,n2,init_ali,bound,neib);
	//--[2] alignment score
	vector <vector <double> > BOUND_score;
	BOUND_score.resize(n1);
	for(i=0;i<n1;i++)
	{
		//current bound
		long wsstart=bound[i+1].first-1;
		long wsend=bound[i+1].second-1;
		//process
		long start=wsstart<0?0:wsstart;
		long end=wsend>=n2?n2-1:wsend;
		for (j = start; j <= end; j++) 
		{
			double score=DNA_Calc(seq_[i],ami_[j]);
			BOUND_score[i].push_back(score);
		}
	}
	//--[3] bounded dyna_prog
	double sco;
	matchs=Advance_Align_Dyna_Prog_Double_Bound(n1,n2,BOUND_score,-2,-1,-2,-1,-1,-1,-1,-1,
		bound,WWW_alignment,sco);
	return sco;
}


//----- assign ali1 and ali2 ------//
void AliPair_To_Ali1_Ali2_Init(long ini_ii,long ini_jj,
	vector<pair<long,long> > & alignment_in,
	vector <long> &ali1,vector <long> &ali2)
{
	//proc
	for(long i=0;i<(long)alignment_in.size();i++)
	{
		long ii=alignment_in[i].first;
		long jj=alignment_in[i].second;
		if(ii>0 && jj>0)
		{
			ali1[ini_ii+ii-1]=ini_jj+jj-1;
			ali2[ini_jj+jj-1]=ini_ii+ii-1;
		}
	}
}



//============= Output Alignment ========//
//-> fasta_output_simp
void FASTA_Output_Simp(FILE *fp,string &nam1,string &nam2,
	const char *ami1,const char *ami2,
	vector<pair<long,long> > &alignment,
	string &assess,
	const char *ami1_rev, vector <long> &ali1_rev)
{
	//alignment->output
	//output
	char c;
	long i;
	long ii,jj;
	long size=alignment.size();
	long n1=(long)strlen(ami1);
	long n2=(long)strlen(ami2);


	for(i=0;i<size;i++)
	{
		//-> output seq1
		ii=alignment[i].first;
		if(ii<=0)fprintf(fp,"-");
		else
		{
			c=ami1[ii-1];
			fprintf(fp,"%c",c);
		}
		//-> output seq2
		jj=alignment[i].second;
		if(jj<=0)fprintf(fp,"-");
		else
		{
			c=ami2[jj-1];
			fprintf(fp,"%c",c);
		}
		//-> output assess
		fprintf(fp,"%c",assess[i]);
		//-> output rev
		if(ii>0)
		{
			fprintf(fp," %c%d\n",ami1_rev[ii-1],ali1_rev[ii-1]);
		}
		else
		{
			fprintf(fp," __\n");
		}
	}

/*
	//--> output seq1
	fprintf(fp,">%s\n",nam1.c_str());
	for(i=0;i<size;i++)
	{
		ii=alignment[i].first;
		if(ii<=0)fprintf(fp,"-");
		else
		{
			c=ami1[ii-1];
			fprintf(fp,"%c",c);
		}
	}
	fprintf(fp,"\n");
	//--> output assess
	if(assess!="")
	{
		fprintf(fp,"%s\n",assess.c_str());
	}
	//--> output seq2
	fprintf(fp,">%s\n",nam2.c_str());
	for(i=0;i<size;i++)
	{
		jj=alignment[i].second;
		if(jj<=0)fprintf(fp,"-");
		else
		{
			c=ami2[jj-1];
			fprintf(fp,"%c",c);
		}
	}
	fprintf(fp,"\n");
*/
}


//------------ cwDTW_align --------------//
double cwDTW_align(string &seqres1,string &seqres2,
	long start1,long start2,long len1,long len2,
	vector<double> &in1,vector<double> &in2,
	vector<pair<long,long> > &WWW_alignment,
	int level,double scale,long radius,long neib)
{
	//-> extract signal
	vector<double> sig1(in1.begin()+start1,in1.begin()+start1+len1);
	vector<double> sig2(in2.begin()+start2,in2.begin()+start2+len2);

	//-> Z-normalization 
	ZScoreNormalize(sig1);
	ZScoreNormalize(sig2);

	//-> continous wavelet transform
	vector<vector<double> > rcwt, pcwt;
	int npyr = level;         // default: 3         
	double scale0 = scale;	  // default: sqrt(2)
	double dscale = 1;        // default: 1
	CWTAnalysis(sig1, rcwt, scale0, dscale, npyr);	
	CWTAnalysis(sig2, pcwt, scale0, dscale, npyr);

	//-> multi-level WaveletDTW 
	double tdiff;
	vector<pair<long,long> > cosali;
	MultiLevel_WaveletDTW(sig1, sig2, rcwt, pcwt, cosali, radius, 0, 0, &tdiff);

	//-> boundary generation
	vector<pair<long,long> > bound;
	BoundGeneration(cosali, neib, bound, 0, 1);

	//-> boundary DynaProg
	long n1=seqres1.length();
	long n2=seqres2.length();
	vector <vector <double> > BOUND_score;
	BOUND_score.resize(n1);
	for(long i=0;i<n1;i++)
	{
		//current bound
		long wsstart=bound[i+1].first-1;
		long wsend=bound[i+1].second-1;
		//process
		long start=wsstart<0?0:wsstart;
		long end=wsend>=n2?n2-1:wsend;
		for (long j = start; j <= end; j++)
		{
			double score=DNA_Calc(seqres1[i],seqres2[j]);
			BOUND_score[i].push_back(score);
		}
	}
	double sco;
	long matchs=Advance_Align_Dyna_Prog_Double_Bound(n1,n2,BOUND_score,-2,-1,-2,-1,-1,-1,-1,-1,
		bound,WWW_alignment,sco);
	//return
	return sco;
}



//----------------- main -------------------//
int main(int argc, char **argv)
{
	struct options opts;
	opts.radius  = 15;
	opts.level   = 4;
	opts.scale0  = sqrt(2);
	opts.neib    = 25;
	opts.verbose = 0;       //-> [0] no verbose; 1 verbose
	opts.test    = 0;       //-> [0] not use test mode; 1 equal_ave, 2 peak_ave, 3 Fast_DTW
	opts.mode    = 0;       //-> [0] block bound; 1 diagonol bound
	opts.kmer    = 0;       //-> [0] to use 5mer; 1 to use 6mer
	

	//----- parse arguments -----//
	if(GetOpts(argc, argv, &opts) < 0){
		EX_TRACE("**WRONG INPUT!**\n");
		return -1;
	}

	std::string input1=opts.input;
	std::string input2=opts.peer;
	std::string output=opts.output;
	if(input1=="" || input2=="")
	{
		fprintf(stderr,"input1 or input2 is NULL \n");
		return -1;
	}


	//======================= START Procedure ===================================//

	//=========================================//
	//----- 1. read genome translated signal -----//
	std::vector<char> genomes1;   //genome sequence
	if(!g::io::ReadATCG(opts.input, genomes1)){
		EX_TRACE("Cannot open %s.\n", opts.input);
		return -1;
	}
	//----- 1.1 pore_model: transform genome sequence to expected signal -------//
	std::vector<int> refer_orig;
	std::vector<double> reference;  //reference: genome signal 1
	Genomes2SignalSequence(genomes1, refer_orig, reference, 1, opts.kmer, 1);
	//----- 1.2 get input name1 --------//
	std::string genom_name_orig=opts.input;
	std::string genom_name;
	getBaseName(genom_name_orig,genom_name,'/','.');


printf("read seq1 \n");


	//========================================//
	//----- 2. read genome translated signal -----//
	std::vector<char> genomes2;   //genome sequence
	if(!g::io::ReadATCG(opts.peer, genomes2)){
		EX_TRACE("Cannot open %s.\n", opts.peer);
		return -1;
	}
	//----- 2.1 pore_model: transform genome sequence to expected signal -------//
	std::vector<int> peer_orig;
	std::vector<double> peer;  //peer: genome signal 2
	Genomes2SignalSequence(genomes2, peer_orig, peer, 1, opts.kmer, 1);
	//----- 2.2 get input name1 --------//
	std::string signal_name_orig=opts.peer;
	std::string signal_name;
	getBaseName(signal_name_orig,signal_name,'/','.');


printf("read seq2 \n");


//exit(-1);    //-> 2.038047 seconds

	//==================================================//
	//------3. process initial input signals ----------//

	//----- 3.1 Zscore normaliza on both signals -----//
//	g::proc::ZScoreNormalize(reference);
//	g::proc::ZScoreNormalize(peer);

	//----- 3.2 calculate length ratio between input signals -----//	
//	double alpha = (double)peer.size()/reference.size();
//	double alpha = 1;


	//----------------- generate coarsening signal in consideration of the sequence length -----------------------//
	vector<pair<long, double> > sigpeaks1;
	vector<pair<long, double> > sigpeaks2;
	double scale0_ini=16*sqrt(2);    //-> compress ~40 times
	//-> process sigpeaks1
	{
		std::vector<std::vector<double> > rcwt;
		long npyr = 1;                   // default: 1
		double dscale = 1;              // default: 1
		CWTAnalysis(reference, rcwt, scale0_ini, dscale, npyr);
		std::vector<std::pair<long, double> > sigpeaks;
		g::proc::PeakPick(rcwt[0], sigpeaks1);
	}

printf("compress seq1 \n");

	//-> process sigpeaks2
	{
		std::vector<std::vector<double> > rcwt;
		long npyr = 1;                   // default: 1
		double dscale = 1;              // default: 1
		CWTAnalysis(peer, rcwt, scale0_ini, dscale, npyr);
		std::vector<std::pair<long, double> > sigpeaks;
		g::proc::PeakPick(rcwt[0], sigpeaks2);
	}

printf("compress seq2 \n");

	//------- perform coarse-grained cwDTW -----------------//
	vector <double> sigpeaks1_sig(sigpeaks1.size(),0);
	vector <double> sigpeaks2_sig(sigpeaks2.size(),0);
	for(long i=0;i<sigpeaks1.size();i++)sigpeaks1_sig[i]=sigpeaks1[i].second;
	for(long i=0;i<sigpeaks2.size();i++)sigpeaks2_sig[i]=sigpeaks2[i].second;

	//---- Z-normalization ----//
	g::proc::ZScoreNormalize(sigpeaks1_sig);
	g::proc::ZScoreNormalize(sigpeaks2_sig);


	//====================================================//
	//----- 4. continous wavelet transform --------------//
	std::vector<std::vector<double> > rcwt, pcwt;

	if(opts.verbose ==1){
		EX_TRACE("CWT Analysis...\n");
	}
	
	int npyr = opts.level;          // default: 3           //---> should be set to 8, or some other value
	double scale0 = opts.scale0;	  // default: sqrt(2)
	double dscale = 1;              // default: 1

	CWTAnalysis(sigpeaks1_sig, rcwt, scale0, dscale, npyr);	
	CWTAnalysis(sigpeaks2_sig, pcwt, scale0, dscale, npyr);


printf("CWTAnalysis done \n");

	//------ 4.1 Zscore normaliza on both CWT signals -----//	
	//if multiscale is used, pyr logical should be added.
//	for(int i = 0; i < npyr; i++){
//		g::proc::ZScoreNormalize(rcwt[i]);
//		g::proc::ZScoreNormalize(pcwt[i]);
//	}

//exit(-1);   //-> 2.162879 seconds

	//============================================//
	//------ 5. multi-level WaveletDTW ----------//
	std::vector<std::pair<long,long> > cosali;
	double tdiff;

	if(opts.verbose ==1){
		EX_TRACE("Coarse Alignment...\n");
	}
	MultiLevel_WaveletDTW(sigpeaks1_sig, sigpeaks2_sig, rcwt, pcwt, cosali, opts.radius, opts.test, opts.mode, &tdiff);
	if(opts.verbose ==1){
		EX_TRACE("Average Deviation (%.1lf/%ld=%.3lf)\n", tdiff, cosali.size(), tdiff/cosali.size());
	}


printf("MultiLevel_WaveletDTW done \n");

	//------ 5.1 generate final boundary -------//
	std::vector<std::pair<long,long> > bound;
	BoundGeneration(cosali, opts.neib, bound, opts.mode, 1);

	//------ 5.2 generate final alignment via cDTW ------//
	std::vector<std::pair<long,long> > alignment;
	tdiff = g::proc::BoundDynamicTimeWarping(sigpeaks1_sig, sigpeaks2_sig, bound, alignment);
	fprintf(stderr,"%s %s %lf %d %lf\n",signal_name.c_str(),genom_name.c_str(),tdiff,alignment.size(),tdiff/alignment.size());

//exit(-1);   //-> 3.386344 seconds

	//=================================================//
	//------ 6. output final alignment to file -------//
//	if(output!="")
//		WriteSequenceAlignment_nano(opts.output, reference, peer, refer_orig, peer_orig, alignment, swap, genome_str);


	//============================== whole genome alignment ==========================//

	//-> 1. read fasta sequence
	string seqres1(genomes1.begin(), genomes1.end() );
	string seqres2(genomes2.begin(), genomes2.end() );
	long n1=seqres1.length();
	long n2=seqres2.length();

	//------------- we set SEQ2 as the REFERENCE sequence -----------------//
	string seqres1_rev;
	string seqres1_revfin=seqres1;
	Reverse_Complement(seqres1,seqres1_rev);
	vector<char> genomes1_rev(seqres1_rev.begin(),seqres1_rev.end());
	vector <long> ali1_rev;  //default: no REV
	ali1_rev.assign(n1,0);
	std::vector<int> refer_orig_rev;
	std::vector<double> reference_rev;  //reference: genome signal 1
	Genomes2SignalSequence(genomes1_rev, refer_orig_rev, reference_rev, 1, opts.kmer, 1);



	//-------- gapless extension -------------//
	vector <long> ali1;
	vector <long> ali2;
	Gapless_Extension(seqres1,seqres2,sigpeaks1,sigpeaks2,alignment,ali1,ali2);

printf("Gapless_Extension done \n");

	//-------- transfer to AliPair -----------//
	vector <vector <long> > AFP_Cor;
	long thres=(long)scale0_ini;
	Ali2_To_AfpCor(n1,n2,ali2,AFP_Cor, thres);

printf("Ali2_To_AfpCor done \n");

	AfpCor_To_Ali1_Ali2(n1,n2,AFP_Cor,ali1,ali2);

printf("AfpCor_To_Ali1_Ali2 done \n");


	//-------- process each "big gap" --------//
	long neib_min=3;
	double neib_ratio=1;
	for(long i=0;i<AFP_Cor.size()-1;i++)
	{
		//--> get position
		long start1=AFP_Cor[i][0]+AFP_Cor[i][2];
		long start2=AFP_Cor[i][1]+AFP_Cor[i][2];
		long end1=AFP_Cor[i+1][0]-1;
		long end2=AFP_Cor[i+1][1]-1;
		//--> extract sequence
		string str1,str2;
		str1=seqres1.substr(start1,end1-start1+1);
		str2=seqres2.substr(start2,end2-start2+1);
		string str1_rev;
		int start1_rev=n1-start1-1;
		int end1_rev=n1-end1-1;
		str1_rev=seqres1_rev.substr(end1_rev,start1_rev-end1_rev+1);
		//--> check length ratio
		long size1=str1.length();
		long size2=str2.length();
		long size_small=size1<size2?size1:size2;
		long size_big=size1>size2?size1:size2;
		if(size_small==0)continue;

		//--> check bad
		{
			double r1=Calc_N_Ratio(str1);
			double r2=Calc_N_Ratio(str2);
			double r1_rev=Calc_N_Ratio(str1_rev);

if(r1!=r1_rev)
{
	fprintf(stderr,"r1 %lf not equal to r1_rev %lf \n",r1,r1_rev);
	exit(-1);
}

			if(r1>0.8 || r2>0.8)continue;
		}

		//--> perform alignment
		double sco_for,sco_rev;
		int rev_or_not=0;  //default: no REV
		vector<pair<long,long> > subali;
		vector<pair<long,long> > subali_rev;
		if(size_big>1000000)
		{
			printf("size_big %d , size_small %d > 1000000 \n",size_big,size_small);
			if(1.0*size_big/size_small>1.2)continue;
			else
			{
				sco_for=cwDTW_align(str1,str2,start1,start2,end1-start1+1,end2-start2+1,
					reference,peer,subali,16,sqrt(2),5,1);
				sco_rev=sco_for-9999;
//				if(sco_for<0 && 1.0*size_big/size_small<1.2)
				if(sco_for<0)
				{
					sco_rev=cwDTW_align(str1_rev,str2,end1_rev,start2,start1_rev-end1_rev+1,end2-start2+1,
						reference_rev,peer,subali_rev,16,sqrt(2),5,1);
				}
				goto next;
			}
		}

//		if(size_big>2*thres && 1.0*size_big/size_small>2.0)continue;
		//--> perform DTW
		if(size_big>1000 && size_small>400)
		{
			sco_for=cwDTW_align(str1,str2,start1,start2,end1-start1+1,end2-start2+1,
				reference,peer,subali,opts.level,opts.scale0,opts.radius,opts.neib);
			sco_rev=sco_for-9999;
//			if(sco_for<0 && 1.0*size_big/size_small<1.2)
			if(sco_for<0)
			{
//				printf("%d / %d -> %lf %lf %d %d -> begin\n",
//					i,AFP_Cor.size(),sco_for,sco_rev,size_small,size_big);
				sco_rev=cwDTW_align(str1_rev,str2,end1_rev,start2,start1_rev-end1_rev+1,end2-start2+1,
					reference_rev,peer,subali_rev,opts.level,opts.scale0,opts.radius,opts.neib);
//				printf("%d / %d -> %lf %lf %d %d -> end\n",
//					i,AFP_Cor.size(),sco_for,sco_rev,size_small,size_big);
			}
		}
		else if(size_big<1000)
		{
			long matchs;
			sco_for=process_oriami_record_simp(str1,str2,subali,matchs,neib_min,neib_ratio);
			sco_rev=sco_for-9999;
//			if(sco_for<0 && 1.0*size_big/size_small<1.2 && size_small>50)
			if(sco_for<0 && size_small>50)
			{
				sco_rev=process_oriami_record_simp(str1_rev,str2,subali_rev,matchs,neib_min,neib_ratio);
			}
		}
		else continue;
		//--> remap index
next:

		//check
		if(size_small>10000)
		{
			printf("%d / %d -> %lf %lf %d %d \n",
				i,AFP_Cor.size(),sco_for,sco_rev,size_small,size_big);
		}
		if(sco_for<sco_rev && sco_rev>0)
		{
			rev_or_not=1;
			subali=subali_rev;
		}
		//rev
		if(rev_or_not==1)
		{
			long count=0;
			for(long k=start1;k<=end1;k++)
			{
				ali1_rev[k]=1;
				seqres1_revfin[k]=str1_rev[count];
				count++;
			}
			printf("REV FOUND %d / %d -> %lf %lf %d %d !!!  ->  %d %d\n",
				i,AFP_Cor.size(),sco_for,sco_rev,size_small,size_big,start1,end1-start1+1);
		}
		AliPair_To_Ali1_Ali2_Init(start1,start2,subali,ali1,ali2);

//if(i%100==0)printf("now %d / %d \r",i,AFP_Cor.size());

	}

printf("Filling Gap done \n");

	//----- final alignment -----//
	vector<pair<long,long> > WWW_alignment;
	Ali1_To_AliPair(n1,n2,ali1,WWW_alignment);

printf("Ali1_To_AliPair done \n");

	//======================== DNA_DynaProg_NanoPore_assess_FixFast ======================//
	{
/*
		//-> 1. read fasta sequence
		std::string seqres1(genomes1.begin(), genomes1.end() );
		std::string seqres2(genomes2.begin(), genomes2.end() );
		long n1=(int)seqres1.length();
		long n2=(int)seqres2.length();

		//-> 2. dynamic programming alignment
		vector <vector <double> > BOUND_score;
		BOUND_score.resize(n1);
		for(long i=0;i<n1;i++)
		{
			//current bound
			long wsstart=bound[i+1].first-1;
			long wsend=bound[i+1].second-1;
			//process
			long start=wsstart<0?0:wsstart;
			long end=wsend>=n2?n2-1:wsend;
			for (long j = start; j <= end; j++) 
			{
				double score=DNA_Calc(seqres1[i],seqres2[j]);
				BOUND_score[i].push_back(score);
			}
		}
		double sco;
		vector<pair<long,long> > WWW_alignment;
		int matchs=Advance_Align_Dyna_Prog_Double_Bound(n1,n2,BOUND_score,-2,-1,-2,-1,0,0,0,0,
			bound,WWW_alignment,sco);
*/

		//-- remove head and tail gap --//
		long iden=0;
		long igap=0;
		long ilen=0;
		long ogap=0;
		long mism=0;
		int first=1;
		string assess="";
		{
			//--> determine start
			long start=-1;
			for(long i=0;i<(long)WWW_alignment.size();i++)
			{
				long ii=WWW_alignment[i].first;
				long jj=WWW_alignment[i].second;
				if(ii>0 && jj>0)
				{
					start=i;
					break;
				}
			}
			//--> determine end
			long end=-1;
			for(long i=(long)WWW_alignment.size()-1;i>=0;i--)
			{
				long ii=WWW_alignment[i].first;
				long jj=WWW_alignment[i].second;
				if(ii>0 && jj>0)
				{
					end=i;
					break;
				}
			}
			//--> final check
			if(start==-1 || end==-1)
			{
				fprintf(stderr,"BAD HERE !!! Null Position Aligned !! \n");
				exit(-1);
			}
			//--> determine
			for(long i=0;i<start;i++)assess.push_back('.');
			for(long i=start;i<=end;i++)
			{
				long ii=WWW_alignment[i].first;
				long jj=WWW_alignment[i].second;
				if(ii>0 && jj>0)
				{
					if(seqres1_revfin[ii-1]==seqres2[jj-1])
					{
						assess.push_back('|');
						iden++;
					}
					else
					{
						assess.push_back('X');
						mism++;
					}
					//-- count gap open --//
					first=1;
				}
				else
				{
					assess.push_back('-');
					igap++;
					//-- count gap open --//
					if(first==1)
					{
						first=0;
						ogap++;
					}
				}
				ilen++;
			}
			for(long i=end+1;i<(long)WWW_alignment.size();i++)assess.push_back('.');
		}

		//-> 3. output alignment
		string nam1=genom_name;
		string nam2=signal_name;
		if(output!="")
		{
			string ali_out=output;
			FILE *fp=fopen(ali_out.c_str(),"wb");
			FASTA_Output_Simp(fp,nam1,nam2,seqres1.c_str(),seqres2.c_str(),WWW_alignment,assess,seqres1_revfin.c_str(),ali1_rev);
			fclose(fp);
		}
		printf("%s %s -> %d -> %d %d %d -> %lf %lf -> %d/%d(%lf) | %d/%d(%lf) -> %d/%d(%lf) | %d/%d(%lf) \n",
			nam1.c_str(),nam2.c_str(),0,
			iden,seqres1.length(),seqres2.length(),
			1.0*iden/seqres1.length(),1.0*iden/seqres2.length(),
			iden,ilen,1.0*iden/ilen,mism,ilen,1.0*mism/ilen,
			igap,ilen,1.0*igap/ilen,ogap,ilen,1.0*ogap/ilen);
	}

	//----- exit -----//	
	return 0;
}

