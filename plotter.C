#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <time.h>
#include <math.h>
#include <cmath>
#include "TRandom3.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include <TLegend.h>
#include "TColor.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TClonesArray.h"

#include <sys/types.h>
#include <dirent.h>

#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <TVector3.h>
#include <TMultiGraph.h>
#include <TH2Poly.h>
#include <TLine.h>

using namespace std;using namespace ROOT::Math;

struct runProperties {
    int runno;
    int runname;
    int datetime;
    float ratioValue;
    float ratioValue1;
    float ratioValue3;
    float ratioError;
    float ratioNoCuts;
    float pedestalValue1;
    float pedestalValue2;
};
runProperties singleRun;
vector <runProperties> allRuns;

struct channel {
    int index;
    vector<runProperties> runForOneChannel;
    TGraphErrors* ratioPlot;
    TGraph* ratioPlot_1;
    TGraph* ratioPlot_2;
    TGraph* ratioPlot_3;
    TGraph* ratioCut;
    TGraph* ratioNoCut;
    TGraph* pedestalGraph1;
    TGraph* pedestalGraph2;
    TGraph* pedestalNormalized1;
    TGraph* pedestalNormalized2;
    TGraph* cutRelErrors;
    TGraph* dtRelErrors;
    TGraph* day0Errors;

    float sysTDCerror;
    float sysDTerror;

    TGraphErrors* ratioEtas;
    TGraph2D* ratio2D;
};

struct edata {
	vector <Int_t> *ieta;
	vector <Int_t> *iphi;
	vector <Int_t> *depth;
	vector <vector <Int_t>> *pulse;
	vector <vector <Int_t>> *tdc;
	vector <vector <Int_t>> *capid;
};
edata ed;

struct semapin {
	int crate;
	int slot;
	int fiber;
	int channel;
	string PMTname;
	string Channelname;
	string winchester;
	int qiecard;
	string qiecrate;
	int qieslot;
	int PMTid;
	int BBid;
	int ieta;
	int iphi;
	int depth;
	int box;
	string boxname;
	string boxbarcode;
	int VA;
	int VB;
	int VC;

	float SX5gain;
	float SX5gainerr;
	float B904pedestal;
	float B904pedestalsigma;
};
semapin semin;

struct semap {
	int crate;
	int slot;
	int fiber;
	int channel;
	string PMTname;
	string Channelname;
	string winchester;
	int qiecard;
	int ieta;
	int iphi;
	int depth;
	int box;
	string boxname;

	float SX5gain;
	float SX5gainerr;
	float B904pedestal;
};
semap sem;

//array of semap structures
vector <semap> SEM;

struct ledlist {
	int mapind;
	float Qmean;
	float Qsigma;
	float npemean;
	float npeerr;
	float gain;
	float gainerr;
	float nbadcapid;
	float qf[32];//??
	float ng[32];//??

    /*TH1F == 1D histogram with floats*/

    TH1F* p;//pulse
	TH1F* t;//tdc
    TH1F* Q1;//Q1 histogram
    TH1F* Q2;//Q2 histogram
    TH1F* Qtot;//Qtot
    TH1F* Qpedest; //Qpedestal (average of charges in TS0 and TS1) histogram
    TH1F* TDC_TS2;//TDC 1 histogram
    TH1F* TDC_TS3;//TDC 2 histogram
    TH2F* Q1vsQ2; //Charges of ts2 and ts3
    TH2F* Q1vsTDC; //Charge - TDC histogram for Q1
    TH2F* Q2vsTDC; //Charge - TDC histogram for Q2
    TH2F* Q2_Q1vsTDC1; // Q2/Q1 vs TDC TS2 2D histogram
    TH2F* Q2_Q1vsTDC2; // Q2/Q1 vs TDC TS3 2D histogram
    TH2F* QtotvsTDC1; // Qtot vs TDC TS2 2D histogram
    TH2F* QtotvsTDC2; // Qtot vs TDC TS3 2D histogram
    TH2F* TDC1vsTDC2; // TDC1 vs TDC2 2D histogram
    TH2F* Q2_Q1vsQ1; // Q2/Q1 vs Q1 2D histogram
    TH2F* Q2_Q1vsQ2; // Q2/Q1 vs Q2 2D histogram
    TH2F* Q2_Q1vsQtot; // Q2/Q1 vs Qtotal 2D histogram
    TH2F* Q2_Q1vsQsum; // Q2/Q1 vs Qtotal 2D histogram

    TH2F* tdc1_1vstdc1_2_first;
    TH2F* tdc2_1vstdc2_2_first;
    TH2F* tdc1_1vstdc1_2;
    TH2F* tdc2_1vstdc2_2;
    TH2F* tdc1_1vstdc1_2_third;
    TH2F* tdc2_1vstdc2_2_third;
    TH2F* q1_1vsq1_2; //Comparison of 2 PMT Channels charge values for TS2
    TH2F* q2_1vsq2_2; //Comparison of 2 PMT Channels charge values for TS3

    TH1F* Q2overQ1_first;//Q2/Q1 for exact delay time values in selectedTSs array (DELAY TIME ~3)
    TH1F* Q2overQ1;//Q2/Q1 for latter delay time values in selectedTSs array (DELAY TIME ~4)
    TH1F* Q2overQ1_third;//delay time 5 for HFP and majority of HFM, delay time 3 for some HFM (DELAY TIME ~5)
    TH1F* Q2overQ1_noCuts;
    TGraph* Q1_Q2overQtotal;//(Q1+Q2)/Qtot Graph

    Double_t q2overq1;
    Double_t q2overq1_nocuts;
    Double_t q2overq1_fi;
    Double_t q2overq1_th;
    Double_t q2overq1_error;
    Double_t qped1 = 0;
    Double_t qped2 = 0;
    int runnumber;

    TH1F* pn;//pulse norm
	TH1F* Q[6][32];//Q
	TGraph* QP[6];//Q profile
    double PT[10][3];//pulse, TDC, norm

};

//from JM 070516
//Converts analogue signal to digital values
float adc2fC_QIE10[256]={

  // =========== RANGE 0 ===========

  // --------- subrange 1 ---------
  -14.45,-11.35,-8.25,-5.15,-2.05,1.05,4.15,7.25,10.35,13.45,16.55,19.65,22.75,25.85,28.95,32.05,
  // --------- subrange 2 ---------
  36.7,42.9,49.1,55.3,61.5,67.7,73.9,80.1,86.3,92.5,98.7,104.9,111.1,117.3,123.5,129.7,135.9,142.1,148.3,154.5,
  // --------- subrange 3 ---------
  163.8,176.2,188.6,201.0,213.4,225.8,238.2,250.6,263.0,275.4,287.8,300.2,312.6,325.0,337.4,349.8,362.2,374.6,387.0,399.4,411.8,
  // --------- subrange 4 ---------
  430.4,455.2,480.0,504.8,529.6,554.4,579.2,
  // =========== RANGE 1 ===========

  // --------- subrange 1 ---------
  529.4,554.2,579.0,603.8,628.6,653.4,678.2,703.0,727.8,752.6,777.4,802.2,827.0,851.8,876.6,901.4,
  // --------- subrange 2 ---------
  938.6,988.2,1037.8,1087.4,1137.0,1186.6,1236.2,1285.8,1335.4,1385.0,1434.6,1484.2,1533.8,1583.4,1633.0,1682.6,1732.2,1781.8,1831.4,1881.0,
  // --------- subrange 3 ---------
  1955.4,2054.6,2153.8,2253.0,2352.2,2451.4,2550.6,2649.8,2749.0,2848.2,2947.4,3046.6,3145.8,3245.0,3344.2,3443.4,3542.6,3641.8,3741.0,3840.2,3939.4,
  // --------- subrange 4 ---------
  4088.2,4286.6,4485.0,4683.4,4881.8,5080.2,5278.6,
  // =========== RANGE 2 ===========

  // --------- subrange 1 ---------
  4879.2,5077.6,5276.0,5474.4,5672.8,5871.2,6069.6,6268.0,6466.4,6664.8,6863.2,7061.6,7260.0,7458.4,7656.8,7855.2,
  // --------- subrange 2 ---------
  8152.8,8549.6,8946.4,9343.2,9740.0,10136.8,10533.6,10930.4,11327.2,11724.0,12120.8,12517.6,12914.4,13311.2,13708.0,14104.8,14501.6,14898.4,15295.2,15692.0,
  // --------- subrange 3 ---------
  16287.2,17080.8,17874.4,18668.0,19461.6,20255.2,21048.8,21842.4,22636.0,23429.6,24223.2,25016.8,25810.4,26604.0,27397.6,28191.2,28984.8,29778.4,30572.0,31365.6,32159.2,
  // --------- subrange 4 ---------
  33349.6,34936.8,36524.0,38111.2,39698.4,41285.6,42872.8,
  // =========== RANGE 3 ===========

  // --------- subrange 1 ---------
  39693.5,41280.5,42867.5,44454.5,46041.5,47628.5,49215.5,50802.5,52389.5,53976.5,55563.5,57150.5,58737.5,60324.5,61911.5,63498.5,
  // --------- subrange 2 ---------
  65879.0,69053.0,72227.0,75401.0,78575.0,81749.0,84923.0,88097.0,91271.0,94445.0,97619.0,100793.0,103967.0,107141.0,110315.0,113489.0,116663.0,119837.0,123011.0,126185.0,
  // --------- subrange 3 ---------
  130946.0,137294.0,143642.0,149990.0,156338.0,162686.0,169034.0,175382.0,181730.0,188078.0,194426.0,200774.0,207122.0,213470.0,219818.0,226166.0,232514.0,238862.0,245210.0,251558.0,257906.0,
  // --------- subrange 4 ---------
  267428.0,280124.0,292820.0,305516.0,318212.0,330908.0,343604.0

};

int HFMBoxMap[37]={0,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19};

//Initialization of run number, will be defined in main
int RunNo=0;int RunType=0;

int plindMap[2][50][75]={{{0}}};


//Structure of a raddam channel
typedef struct {
    int ieta;
    int iphi;
    int depth;
} Raddam_ch;

//Array that contains 56 raddam channels ieta, iphi and idepth
static const Raddam_ch RADDAM_CH[56]= {{30,1,2},{30,21,1},{30,37,2},{30,57,1},{32,1,1},{32,21,2},{32,37,1},{32,57,2},{34,1,2},{34,21,1},{34,37,2},{34,57,1},{36,1,1},{36,21,2 },{36,37,1},{36,57,2},{38,1,2},{38,21,1},{38,37,2},{38,57,1},{40,19,2},{40,35,1},{40,55,2},{40,71,1},{41,19,1},{41,35,2},{41,55,1},{41,71,2},{-30,15,2},{-30,35,1},{-30,51,2},{-30,71,1},{-32,15,1},{-32,35,2},{-32,51,1},{-32,71,2},{-34,15,2},{-34,35,1},{-34,51,2},{-34,71,1},{-36,15,1},{-36,35,2},{-36,51,1},{-36,71,2},{-38,15,2},{-38,35,1},{-38,51,2},{-38,71,1}, {-40,15,1},{-40,35,2},{-40,51,1},{-40,71,2},{-41,15,2},{-41,35,1},{-41,51,2},{-41,71,1}};

//TDC2 cuts for delay time 5 for HFP and majority of HFM, delay time 3 for some HFM  (DELAY TIME ~5)
static const int tdcCuts3_TS3[56][2]= {{22,25},{60,64},{16,20},{60,64},     {18,22},{60,64},{20,26},{13,18},    {19,22},{60,64},{60,64},{12,14},    {22,28},{14,17},{60,64},{60,64},    {20,27},{60,64},{60,64},{60,64},     {60,64},{20,25},{60,64},{60,64},    {60,64},{20,30},{60,64},{60,64},    {60,64},{60,64},{60,64},{60,64},    {60,64},{60,64},{60,64},{60,64},{60,64},{60,64},{60,64},{60,64},    {15,19},{60,64},{60,64},{60,64},    {60,64},{60,64},{60,64},{60,64},    {60,64},{60,64},{60,64},{60,64},    {60,64},{60,64},{60,64},{60,64}};

//TDC2 cuts for latter of delay time values in selectedTSs array  (DELAY TIME ~4)
static const int tdcCuts2_TS3[56][2]= {{12,16},{60,64},{6,10},{0,8},    {10,12},{60,64},{12,16},{5,7},  {10,12},{60,64},{13,15},{2,5},    {13,17},{5,7},{13,15},{60,64},  {13,16},{60,64},{60,64},{10,12},  {12,14},{12,14},{11,13},{16,18},    {12,13},{15,17},{60,64},{16,22},    {11,13},{60,64},{60,64},{60,64},    {60,64},{60,64},{60,64},{60,64},    {60,64},{60,64},{60,64},{60,64},    {6,9},{0,64},{60,64},{60,64},    {13,15},{60,64},{60,64},{60,64},    {14,16},{60,64},{60,64},{60,64},    {13,15},{60,64},{60,64},{60,64}};

//TDC1 cuts for exact latter of time values in selectedTSs array  (DELAY TIME ~4)
static const int tdcCuts2_TS2[56][2]= {{10,16},{0,64},{0,64},{0,64},    {0,64},{0,64},{0,64},{0,64},  {0,64},{0,64},{0,64},{0,64},   {0,64},{0,64},{0,64},{0,64},  {0,64},{0,64},{0,64},{0,64},  {0,64},{0,64},{0,64},{0,64},    {0,64},{14,18},{0,64},{0,64},  {0,64},{0,64},{0,64},{0,64},   {0,64},{0,64},{0,64},{0,64},  {0,64},{0,64},{0,64},{0,64},    {0,64},{0,64},{0,64},{0,64},    {0,64},{0,64},{0,64},{0,64},  {0,64},{0,64},{0,64},{0,64},   {0,64},{0,64},{0,64},{0,64}};

//TDC2 cuts for exact delay time values in selectedTSs array  (DELAY TIME ~3)
static const int tdcCuts1_TS3[56][2]= {{4,6},{60,64},{60,64},{60,64},    {0,3},{60,64},{4,6},{60,64},  {0,3},{60,64},{4,7},{60,64},    {5,7},{60,64},{2,6},{60,64},  {4,6},{60,64},{60,64},{0,2},  {1,5},{3,5},{0,3},{6,10},    {0,5},{4,8},{60,64},{7,14},    {0,6},{0,64},{60,64},{60,64},    {60,64},{60,64},{60,64},{60,64},    {60,64},{60,64},{60,64},{60,64},    {60,64},{60,64},{60,64},{60,64},    {4,6},{60,64},{60,64},{60,64},    {5,6},{60,64},{60,64},{60,64},    {4,7},{60,64},{60,64},{60,64}};

//TDC1 cuts for exact delay time values in selectedTSs array  (DELAY TIME ~3)
static const int tdcCuts1_TS2[56][2]= {{2,6},{0,4},{60,64},{0,3},    {0,2},{60,64},{4,7},{60,64},  {0,3},{6,8},{4,7},{60,64},    {4,6},{60,64},{3,5},{2,4},  {3,6},{4,6},{60,64},{0,4},  {2,6},{4,6},{0,4},{0,64},    {2,6},{4,8},{0,4},{6,10},    {0,7},{0,64},{0,2},{4,8},    {0,10},{0,10},{0,10},{0,10},    {0,10},{0,64},{0,10},{0,10},    {60,64},{0,10},{60,64},{0,10},    {0,10},{0,10},{0,10},{0,10},    {0,10},{0,10},{0,10},{0,64},    {0,10},{0,10},{0,10},{0,10}};

static const int selectedTSs[56]= {3,3,3,3, 3,3,4,3, 3,4,4,3, 4,3,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,2, 3,3,3,2, 3,4,4,2, 4,4,4,2, 3,4,3,2, 4,4,4,2, 4,4,4,2, 4,4,4,4};

int getmap() {

    //reads the file semapex_v3.txt
	ifstream semapinf("semapex_v3.txt");
	while(!semapinf.eof()) {
        //Reads formatted input from the file, and initializes an object semin
				semapinf>>semin.crate>>semin.slot>>semin.fiber>>semin.channel>>semin.PMTname>>semin.Channelname>>semin.winchester>>semin.qiecard>>semin.qiecrate>>semin.qieslot>>semin.PMTid>>semin.BBid>>semin.ieta>>semin.iphi>>semin.depth>>semin.box>>semin.boxname>>semin.boxbarcode>>semin.VA>>semin.VB>>semin.VC>>semin.SX5gain>>semin.SX5gainerr>>semin.B904pedestal>>semin.B904pedestalsigma;

		sem.crate=semin.crate;
		sem.slot=semin.slot;
		sem.fiber=semin.fiber;
		sem.channel=semin.channel;
		sem.PMTname=semin.PMTname;
		sem.Channelname=semin.Channelname;
		sem.winchester=semin.winchester;
		sem.qiecard=semin.qiecard;
		sem.ieta=semin.ieta;
		sem.iphi=semin.iphi;
		sem.depth=semin.depth;
		sem.box=semin.box;
		sem.boxname=semin.boxname;
		sem.SX5gain=semin.SX5gain;
		sem.SX5gainerr=semin.SX5gainerr;
		sem.B904pedestal=semin.B904pedestal;
		SEM.push_back(sem);
	}
	semapinf.close();
}

/*
 * Reads the NTuple file, get events and channel info, initializes histograms,
 * select events and fills histograms, writes the results to a file
 */
int plotleds() {

	char hname[500];char hname2[500];
	sprintf(hname,"../NTuples/N_%d.root",RunNo); //Input file that is created by ntupler, i.e. config file
	TFile* inroot=new TFile(hname); //Root binary file that includes events
	TTree *tree = (TTree*)inroot->Get("Events"); //Reads file and gets events
	tree->SetBranchAddress("ieta",&ed.ieta);
	tree->SetBranchAddress("iphi",&ed.iphi);
	tree->SetBranchAddress("depth",&ed.depth);
	tree->SetBranchAddress("pulse",&ed.pulse);
	tree->SetBranchAddress("tdc",&ed.tdc);
	tree->SetBranchAddress("capid",&ed.capid);

	//int ev1=e1;int ev2=e2;
	//if(ev1==-1) {
    int ev1=0;
    int ev2=tree->GetEntries();
    //}

	ledlist pl; //Contains all histograms in ledlist structure for a single channel
    vector <ledlist> PL; // list of pl structure of all channels
    int plind=-1;

	tree->GetEntry(0);
	int nTS=ed.pulse->at(0).size(); // number of time slices

	int qpcols[6]={9,1,2,3,4,6};//???

    //HISTO INITIALIZATIONS
	for(int i1=0;i1<ed.ieta->size();i1++) {
		plind=-1;

		for(int i2=0;i2<PL.size();i2++) {

			if(SEM[PL[i2].mapind].ieta==ed.ieta->at(i1) && SEM[PL[i2].mapind].iphi==ed.iphi->at(i1)){
				plind=i2;
                break;
			}
		}
		if(plind==-1) {

			pl.mapind=-1;
			for(int i2=0;i2<SEM.size();i2++) {

				if(SEM[i2].ieta==ed.ieta->at(i1) && SEM[i2].iphi==ed.iphi->at(i1)) {
					pl.mapind=i2;
                    break;
				}
			}
			pl.Qmean=0.;
			pl.Qsigma=0.;
			pl.npemean=0.;
			pl.npeerr=0.;
			pl.gain=0.;
			pl.gainerr=0.;
			pl.nbadcapid=0.;

            pl.runnumber = RunNo; //Run number

			for(int ik1=0;ik1<32;ik1++) {
				pl.qf[ik1]=0.;
				pl.ng[ik1]=0.;
			}

			/* 2D Histogram of Pulse Shapes of Q1 and Q2 */
            sprintf(hname, "Q1 vs Q2 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q1vsQ2 = new TH2F(hname,hname,200,0.,2000.,200,0.,2000.);
            pl.Q1vsQ2->SetCanExtend(TH2::kAllAxes);
            pl.Q1vsQ2->GetXaxis()->SetTitle("Pulse Q1");pl.Q1vsQ2->GetXaxis()->CenterTitle();
            pl.Q1vsQ2->GetYaxis()->SetTitle("Pulse Q2");pl.Q1vsQ2->GetYaxis()->CenterTitle();

            /* 2D Histogram of ADC vs TDC for Q1 */
            sprintf(hname, "Q1 vs TDC1 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q1vsTDC = new TH2F(hname,hname,200,0,200,130,0,129);
            pl.Q1vsTDC->SetCanExtend(TH2::kAllAxes);
            pl.Q1vsTDC->GetXaxis()->SetTitle("Q1 (fC)");pl.Q1vsTDC->GetXaxis()->CenterTitle();
            pl.Q1vsTDC->GetYaxis()->SetTitle("TDC1");pl.Q1vsTDC->GetYaxis()->CenterTitle();

            /* 2D Histogram of ADC vs TDC for Q2 */
            sprintf(hname, "Q2 vs TDC2 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q2vsTDC = new TH2F(hname,hname,200,0,200,130,0,129);
            pl.Q2vsTDC->SetCanExtend(TH2::kAllAxes);
            pl.Q2vsTDC->GetXaxis()->SetTitle("Q2 (fC)");pl.Q2vsTDC->GetXaxis()->CenterTitle();
            pl.Q2vsTDC->GetYaxis()->SetTitle("TDC2");pl.Q2vsTDC->GetYaxis()->CenterTitle();

            /* 2D Histogram of Q2/Q1 vs TDC TS2 */
            sprintf(hname, "Q2/Q1 vs TDC1 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q2_Q1vsTDC1 = new TH2F(hname,hname,30,-3,3.,130,0,129);
            pl.Q2_Q1vsTDC1->SetCanExtend(TH2::kAllAxes);
            pl.Q2_Q1vsTDC1->GetXaxis()->SetTitle("Q2/Q1");pl.Q2_Q1vsTDC1->GetXaxis()->CenterTitle();
            pl.Q2_Q1vsTDC1->GetYaxis()->SetTitle("TDC1");pl.Q2_Q1vsTDC1->GetYaxis()->CenterTitle();

            /* 2D Histogram of Q2/Q1 vs TDC TS3 */
            sprintf(hname, "Q2/Q1 vs TDC2 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q2_Q1vsTDC2 = new TH2F(hname,hname,30,-3.,3.,130,0,129);
            pl.Q2_Q1vsTDC2->SetCanExtend(TH2::kAllAxes);
            pl.Q2_Q1vsTDC2->GetXaxis()->SetTitle("Q2/Q1");pl.Q2_Q1vsTDC2->GetXaxis()->CenterTitle();
            pl.Q2_Q1vsTDC2->GetYaxis()->SetTitle("TDC2");pl.Q2_Q1vsTDC2->GetYaxis()->CenterTitle();

            /* 2D Histogram of Q2/Q1 vs Q1 */
            sprintf(hname, "Q2/Q1 vs Q1 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q2_Q1vsQ1 = new TH2F(hname,hname,30,-3,3.,10,0,9);
            pl.Q2_Q1vsQ1->SetCanExtend(TH2::kAllAxes);
            pl.Q2_Q1vsQ1->GetXaxis()->SetTitle("Q2/Q1");pl.Q2_Q1vsQ1->GetXaxis()->CenterTitle();
            pl.Q2_Q1vsQ1->GetYaxis()->SetTitle("Q1 (fC)");pl.Q2_Q1vsQ1->GetYaxis()->CenterTitle();

            /* 2D Histogram of Q2/Q1 vs Q2 */
            sprintf(hname, "Q2/Q1 vs Q2 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q2_Q1vsQ2 = new TH2F(hname,hname,30,-3.,3.,10,0,9);
            pl.Q2_Q1vsQ2->SetCanExtend(TH2::kAllAxes);
            pl.Q2_Q1vsQ2->GetXaxis()->SetTitle("Q2/Q1");pl.Q2_Q1vsQ2->GetXaxis()->CenterTitle();
            pl.Q2_Q1vsQ2->GetYaxis()->SetTitle("Q2 (fC)");pl.Q2_Q1vsQ2->GetYaxis()->CenterTitle();

            /* 2D Histogram of Q2/Q1 vs Qtotal */
            sprintf(hname, "Q2/Q1 vs Qtotal (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q2_Q1vsQtot = new TH2F(hname,hname,30,-3.,3.,10,0,9);
            pl.Q2_Q1vsQtot->SetCanExtend(TH2::kAllAxes);
            pl.Q2_Q1vsQtot->GetXaxis()->SetTitle("Q2/Q1");pl.Q2_Q1vsQtot->GetXaxis()->CenterTitle();
            pl.Q2_Q1vsQtot->GetYaxis()->SetTitle("Qt (fC)");pl.Q2_Q1vsQtot->GetYaxis()->CenterTitle();

            /* 2D Histogram of Q2/Q1 vs Qtotal */
            sprintf(hname, "Q2/Q1 vs (Q1+Q2) (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q2_Q1vsQsum = new TH2F(hname,hname,30,-3.,3.,10,0,9);
            pl.Q2_Q1vsQsum->SetCanExtend(TH2::kAllAxes);
            pl.Q2_Q1vsQsum->GetXaxis()->SetTitle("Q_{2}/Q_{1}");pl.Q2_Q1vsQsum->GetXaxis()->CenterTitle();
            pl.Q2_Q1vsQsum->GetYaxis()->SetTitle("Q1+Q2 (fC)");pl.Q2_Q1vsQsum->GetYaxis()->CenterTitle();

            /* 2D Histogram of Qtot vs TDC TS2 */
            sprintf(hname, "Qt vs TDC1 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.QtotvsTDC1 = new TH2F(hname,hname,21,-10.,10.,130,0,129);
            pl.QtotvsTDC1->SetCanExtend(TH2::kAllAxes);
            pl.QtotvsTDC1->GetXaxis()->SetTitle("Qt (fC)");pl.QtotvsTDC1->GetXaxis()->CenterTitle();
            pl.QtotvsTDC1->GetYaxis()->SetTitle("TDC1");pl.QtotvsTDC1->GetYaxis()->CenterTitle();

            /* 2D Histogram of Qtot vs TDC TS3 */
            sprintf(hname, "Qt vs TDC2 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.QtotvsTDC2 = new TH2F(hname,hname,21,-10.,10.,130,0,129);
            pl.QtotvsTDC2->SetCanExtend(TH2::kAllAxes);
            pl.QtotvsTDC2->GetXaxis()->SetTitle("Qt (fC)");pl.QtotvsTDC2->GetXaxis()->CenterTitle();
            pl.QtotvsTDC2->GetYaxis()->SetTitle("TDC2");pl.QtotvsTDC2->GetYaxis()->CenterTitle();

            /* 2D Histogram of TDC1 vs TDC2 */
            sprintf(hname, "TDC1 vs TDC2 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.TDC1vsTDC2 = new TH2F(hname,hname,130,0,129,130,0,129);
            pl.TDC1vsTDC2->GetXaxis()->SetTitle("TDC1");pl.TDC1vsTDC2->GetXaxis()->CenterTitle();
            pl.TDC1vsTDC2->GetYaxis()->SetTitle("TDC2");pl.TDC1vsTDC2->GetYaxis()->CenterTitle();

            /* Q1 1D histogram */
            sprintf(hname,"Q1 %d %d",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q1=new TH1F(hname,hname,100,0.,1000.);
            pl.Q1->SetCanExtend(TH1::kAllAxes);
            pl.Q1->GetXaxis()->SetTitle("Q_{1}(fC)");pl.Q1->GetXaxis()->CenterTitle();
            pl.Q1->GetYaxis()->SetTitle("Events");pl.Q1->GetYaxis()->CenterTitle();
            pl.Q1->SetLineColor(3); pl.Q1->SetFillColor(3);

            /* Q2 1D histogram */
            sprintf(hname,"Q2 %d %d",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q2=new TH1F(hname,hname,100,0.,1000.);
            pl.Q2->SetCanExtend(TH1::kAllAxes);
            pl.Q2->GetXaxis()->SetTitle("Q_{2}(fC)");pl.Q2->GetXaxis()->CenterTitle();
            pl.Q2->GetYaxis()->SetTitle("Events");pl.Q2->GetYaxis()->CenterTitle();
            pl.Q2->SetLineColor(7); pl.Q2->SetFillColor(7);

            /* Qpedest 1D histogram */
            sprintf(hname,"Qpedestal %d %d",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Qpedest=new TH1F(hname,hname,100,-5.,5.);
            pl.Qpedest->SetCanExtend(TH1::kAllAxes);
            pl.Qpedest->GetXaxis()->SetTitle("Q_{p}(fC)");pl.Qpedest->GetXaxis()->CenterTitle();
            pl.Qpedest->GetYaxis()->SetTitle("Events");pl.Qpedest->GetYaxis()->CenterTitle();
            pl.Qpedest->SetLineColor(51); pl.Qpedest->SetFillColor(51);

            /* TDC 1 1D histogram */
            sprintf(hname,"TDC Q1 %d %d",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.TDC_TS2=new TH1F(hname,hname,130,0,129);
            pl.TDC_TS2->GetXaxis()->SetTitle("TDC_{1}");pl.TDC_TS2->GetXaxis()->CenterTitle();
            pl.TDC_TS2->GetYaxis()->SetTitle("Events");pl.TDC_TS2->GetYaxis()->CenterTitle();
            pl.TDC_TS2->SetLineColor(6);

            /* TDC 2 1D histogram */
            sprintf(hname,"TDC Q2 %d %d",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.TDC_TS3=new TH1F(hname,hname,130,0,129);
            pl.TDC_TS3->GetXaxis()->SetTitle("TDC_{2}");pl.TDC_TS3->GetXaxis()->CenterTitle();
            pl.TDC_TS3->GetYaxis()->SetTitle("Events");pl.TDC_TS3->GetYaxis()->CenterTitle();
            pl.TDC_TS3->SetLineColor(8);

            /* s1/s2 1d histogram */
			sprintf(hname2,"Q_{2}/Q_{1} %d %d",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
			pl.Q2overQ1=new TH1F(hname2,"",100,0.,5.);
            pl.Q2overQ1->SetCanExtend(TH1::kAllAxes);
			pl.Q2overQ1->GetXaxis()->SetTitle("Q_{2}/Q_{1}");pl.Q2overQ1->GetXaxis()->CenterTitle();
			pl.Q2overQ1->GetYaxis()->SetTitle("Events");pl.Q2overQ1->GetYaxis()->CenterTitle();
			pl.Q2overQ1->SetLineColor(30);

            /* s1/s2 1d histogram */
            sprintf(hname2,"Q_{2}/Q_{1} Former %d %d",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q2overQ1_first=new TH1F(hname2,hname2,100,0.,10.);
            pl.Q2overQ1_first->SetCanExtend(TH1::kAllAxes);
            pl.Q2overQ1_first->GetXaxis()->SetTitle("Q_{2}/Q_{1} 1st");pl.Q2overQ1_first->GetXaxis()->CenterTitle();
            pl.Q2overQ1_first->GetYaxis()->SetTitle("Events");pl.Q2overQ1_first->GetYaxis()->CenterTitle();
            pl.Q2overQ1_first->SetLineColor(25);

            /* s1/s2 1d histogram */
            sprintf(hname2,"Q_{2}/Q_{1} Latter %d %d",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q2overQ1_third=new TH1F(hname2,hname2,100,0.,10.);
            pl.Q2overQ1_third->SetCanExtend(TH1::kAllAxes);
            pl.Q2overQ1_third->GetXaxis()->SetTitle("Q_{2}/Q_{1} 3rd");pl.Q2overQ1_third->GetXaxis()->CenterTitle();
            pl.Q2overQ1_third->GetYaxis()->SetTitle("Events");pl.Q2overQ1_third->GetYaxis()->CenterTitle();
            pl.Q2overQ1_third->SetLineColor(35);

            /* s1/s2 1d histogram */
            sprintf(hname2,"Q2/Q1 No TDC Cuts %d %d",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.Q2overQ1_noCuts=new TH1F(hname2,hname2,100,0.,10.);
            pl.Q2overQ1_noCuts->SetCanExtend(TH1::kAllAxes);
            pl.Q2overQ1_noCuts->GetXaxis()->SetTitle("Q_{2}/Q_{1}");pl.Q2overQ1_noCuts->GetXaxis()->CenterTitle();
            pl.Q2overQ1_noCuts->GetYaxis()->SetTitle("Events");pl.Q2overQ1_noCuts->GetYaxis()->CenterTitle();
            pl.Q2overQ1_noCuts->SetLineColor(5);

            /* Qtot 1d histogram */
			sprintf(hname2,"Qt %d %d",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
			pl.Qtot=new TH1F(hname2,hname2,1000,0.,1000.);
            pl.Qtot->SetCanExtend(TH1::kAllAxes);
			pl.Qtot->GetXaxis()->SetTitle("Qt (fC)"); pl.Qtot->GetXaxis()->CenterTitle();
			pl.Qtot->GetYaxis()->SetTitle("Events"); pl.Qtot->GetYaxis()->CenterTitle();
			pl.Qtot->SetLineColor(8);

            /* 2D Histogram of TDC1-1 vs TDC1-2 */
            sprintf(hname, "TDC1-1 vs TDC1-2 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.tdc1_1vstdc1_2 = new TH2F(hname,hname,67,0,66,67,0,66);
            pl.tdc1_1vstdc1_2->GetXaxis()->SetTitle("TDC1-1");pl.tdc1_1vstdc1_2->GetXaxis()->CenterTitle();
            pl.tdc1_1vstdc1_2->GetYaxis()->SetTitle("TDC1-2");pl.tdc1_1vstdc1_2->GetYaxis()->CenterTitle();

            /* 2D Histogram of TDC2-1 vs TDC2-2 */
            sprintf(hname, "TDC2-1 vs TDC2-2 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.tdc2_1vstdc2_2 = new TH2F(hname,hname,67,0,66,67,0,66);
            pl.tdc2_1vstdc2_2->GetXaxis()->SetTitle("TDC2-1");pl.tdc2_1vstdc2_2->GetXaxis()->CenterTitle();
            pl.tdc2_1vstdc2_2->GetYaxis()->SetTitle("TDC2-2");pl.tdc2_1vstdc2_2->GetYaxis()->CenterTitle();

            /* 2D Histogram of TDC1-1 vs TDC1-2 FORMER*/
            sprintf(hname, "TDC1-1 vs TDC1-2 Former (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.tdc1_1vstdc1_2_first = new TH2F(hname,hname,67,0,66,67,0,66);
            pl.tdc1_1vstdc1_2_first->GetXaxis()->SetTitle("TDC1-1");pl.tdc1_1vstdc1_2_first->GetXaxis()->CenterTitle();
            pl.tdc1_1vstdc1_2_first->GetYaxis()->SetTitle("TDC1-2");pl.tdc1_1vstdc1_2_first->GetYaxis()->CenterTitle();

            /* 2D Histogram of TDC2-1 vs TDC2-2 FORMER*/
            sprintf(hname, "TDC2-1 vs TDC2-2 Former (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.tdc2_1vstdc2_2_first = new TH2F(hname,hname,67,0,66,67,0,66);
            pl.tdc2_1vstdc2_2_first->GetXaxis()->SetTitle("TDC2-1");pl.tdc2_1vstdc2_2_first->GetXaxis()->CenterTitle();
            pl.tdc2_1vstdc2_2_first->GetYaxis()->SetTitle("TDC2-2");pl.tdc2_1vstdc2_2_first->GetYaxis()->CenterTitle();

            /* 2D Histogram of TDC1-1 vs TDC1-2 LATTER */
            sprintf(hname, "TDC1-1 vs TDC1-2 Latter (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.tdc1_1vstdc1_2_third = new TH2F(hname,hname,67,0,66,67,0,66);
            pl.tdc1_1vstdc1_2_third->GetXaxis()->SetTitle("TDC1-1");pl.tdc1_1vstdc1_2_third->GetXaxis()->CenterTitle();
            pl.tdc1_1vstdc1_2_third->GetYaxis()->SetTitle("TDC1-2");pl.tdc1_1vstdc1_2_third->GetYaxis()->CenterTitle();

            /* 2D Histogram of TDC2-1 vs TDC2-2 */
            sprintf(hname, "TDC2-1 vs TDC2-2 Latter (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.tdc2_1vstdc2_2_third = new TH2F(hname,hname,67,0,66,67,0,66);
            pl.tdc2_1vstdc2_2_third->GetXaxis()->SetTitle("TDC2-1");pl.tdc2_1vstdc2_2_third->GetXaxis()->CenterTitle();
            pl.tdc2_1vstdc2_2_third->GetYaxis()->SetTitle("TDC2-2");pl.tdc2_1vstdc2_2_third->GetYaxis()->CenterTitle();

            //2D Histogram of Q1-1 vs Q1-2
            sprintf(hname, "Q1-1 vs Q1-2 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.q1_1vsq1_2 = new TH2F(hname,hname,50,0,100,50,0,100);
            pl.q1_1vsq1_2->SetCanExtend(TH2::kAllAxes);
            pl.q1_1vsq1_2->GetXaxis()->SetTitle("Q1-1 (fC)");pl.q1_1vsq1_2->GetXaxis()->CenterTitle();
            pl.q1_1vsq1_2->GetYaxis()->SetTitle("Q1-2 (fC)");pl.q1_1vsq1_2->GetYaxis()->CenterTitle();

            //2D Histogram of Q2-1 vs Q2-2
            sprintf(hname, "Q2-1 vs Q2-2 (%d %d)",SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            pl.q2_1vsq2_2 = new TH2F(hname,hname,100,0,100,100,0,100);
            pl.q2_1vsq2_2->SetCanExtend(TH2::kAllAxes);
            pl.q2_1vsq2_2->GetXaxis()->SetTitle("Q2-1 (fC)");pl.q2_1vsq2_2->GetXaxis()->CenterTitle();
            pl.q2_1vsq2_2->GetYaxis()->SetTitle("Q2-2 (fC)");pl.q2_1vsq2_2->GetYaxis()->CenterTitle();


			for(int is1=0;is1<10;is1++) {
				for(int is2=0;is2<3;is2++) {
					pl.PT[is1][is2]=0.;
				}
			}

            /* Pulse Shapes 1D Histogram */
			sprintf(hname,"#splitline{Pulse %d %s %s}{      (%d, %d)}", SEM[pl.mapind].box,SEM[pl.mapind].PMTname.c_str(),SEM[pl.mapind].Channelname.c_str(),SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
            //Name of the histogram contains box number, PMT name and Channel name sequentially
			sprintf(hname2,"Pulse %d %s %s",SEM[pl.mapind].box,SEM[pl.mapind].PMTname.c_str(),SEM[pl.mapind].Channelname.c_str());
			pl.p=new TH1F(hname2,hname,nTS,-0.5,((float)nTS)-0.5); //Divided the axis to the TS number
			pl.p->GetXaxis()->SetTitle("TS (x25 ns)");pl.p->GetXaxis()->CenterTitle();
			pl.p->GetYaxis()->SetTitle("Mean Charge per TS (fC)");pl.p->GetYaxis()->CenterTitle();
			pl.p->SetFillColor(4);pl.p->SetLineColor(4);

            /* TDC Shapes 1D Histogram */
			sprintf(hname,"#splitline{TDC %d %s %s}{      (%d, %d)}", SEM[pl.mapind].box,SEM[pl.mapind].PMTname.c_str(),SEM[pl.mapind].Channelname.c_str(),SEM[pl.mapind].ieta,SEM[pl.mapind].iphi);
			sprintf(hname2,"TDC %d %s %s",SEM[pl.mapind].box,SEM[pl.mapind].PMTname.c_str(),SEM[pl.mapind].Channelname.c_str());
			pl.t=new TH1F(hname2,hname,nTS,-0.5,((float)nTS)-0.5);
			pl.t->GetXaxis()->SetTitle("TS (x25 ns)");pl.t->GetXaxis()->CenterTitle();
			pl.t->GetYaxis()->SetTitle("Mean TDC per TS");pl.t->GetYaxis()->CenterTitle();
			pl.t->SetLineColor(46);

            //(Q1+Q2)/Qtot graph
            pl.Q1_Q2overQtotal = new TGraph();

			PL.push_back(pl);
		}
	}
    //END OF HISTO INITIALIZATIONS

	int MI=-1;
	for(int i2=0;i2<PL.size();i2++) {

		MI=PL[i2].mapind;
        //HFM plindMap[0] and HFP plindMap[1] = i2
        //plindMap maps the PL vector content to channels (ieta,iphi,idepth)
		plindMap[(SEM[MI].ieta<0)?0:1][abs(SEM[MI].ieta)][SEM[MI].iphi]=i2;
	}

	double ped=0.;double sig=0.;bool capidOK=true;int cap0=0;

    //EVENT SELECTION AND FILLING THE HISTOGRAMS
    for(int i=ev1;i<ev2;i++) {

        tree->GetEntry(i);// gets i'th event

        float tdc1values[PL.size()], tdc2values[PL.size()], tdc1values_1[PL.size()], tdc2values_1[PL.size()], tdc1values_3[PL.size()], tdc2values_3[PL.size()], q1s_1[PL.size()], q2s_1[PL.size()], q1s_3[PL.size()], q2s_3[PL.size()],q1s[PL.size()], q2s[PL.size()],q1s_nocuts[PL.size()], q2s_nocuts[PL.size()], qtots[PL.size()], qpeds[PL.size()] ;
        float NQ1, NQ2, NQp;

        for(int i1=0;i1<ed.ieta->size();i1++){

            plind=plindMap[(ed.ieta->at(i1)<0)?0:1][abs(ed.ieta->at(i1))][ed.iphi->at(i1)];
            ped=0.;sig=0.;capidOK=true;
            float QQ[10]={0.};float Qtot=0.;float qf=0.; //float TDC[10]={0.};

            for(int i2=0;i2<nTS;i2++) {
                QQ[i2]=adc2fC_QIE10[ed.pulse->at(i1)[i2]];
                Qtot += adc2fC_QIE10[ed.pulse->at(i1)[i2]];
            }
            cout << i << " " << (QQ[2]+QQ[3])/Qtot << endl;
            PL[plind].Q1_Q2overQtotal->SetPoint(PL[plind].Q1_Q2overQtotal->GetN(), i, (QQ[2]+QQ[3])/Qtot);

            bool e1, e2, e3;

            for( int chan = 0 ; chan < 56 ; chan++) {
                if ((ed.ieta->at(i1)==RADDAM_CH[chan].ieta) && (ed.iphi->at(i1)==RADDAM_CH[chan].iphi)){
                    if (RADDAM_CH[chan].ieta < 0 && RADDAM_CH[chan].ieta != -41 && RADDAM_CH[chan].iphi == 35 && RunNo > 322851){
                        e1 = (i >= 5000) && i < 10000;
                        e2 = (i >= 10000) && i < 15000;
                        e3 = (i >= 15000) && i < 20000;
                        break;
                    }
                    e1 = (i >= (selectedTSs[chan]-1)*5000 && i < selectedTSs[chan]*5000);
                    e2 = (i >= selectedTSs[chan]*5000 && i < (selectedTSs[chan]+1)*5000);
                    e3 = (i >= (selectedTSs[chan]+1)*5000 && i < (selectedTSs[chan]+2)*5000);
                    break;
                }
            }

            if(e1){ //THE EARLIER ONE

                //Cuts for invalid TDC values:
                bool f1, f2, f3, f ;
                for( int chan = 0 ; chan < 56 ; chan++) {
                    if ((ed.ieta->at(i1)==RADDAM_CH[chan].ieta) && (ed.iphi->at(i1)==RADDAM_CH[chan].iphi)){
                        f1 = ((ed.tdc->at(i1)[3] > tdcCuts1_TS3[chan][1]) || (ed.tdc->at(i1)[3] < tdcCuts1_TS3[chan][0]));
                        f2 = ((ed.tdc->at(i1)[2] > tdcCuts1_TS2[chan][1]) || (ed.tdc->at(i1)[2] < tdcCuts1_TS2[chan][0]));
                        f3 = (QQ[2] < 10) || (QQ[3] < 10);
                        f = f1 || f2 || f3;
                        break;
                    }
                }
                if(!f) {
                    if(i1 < ed.ieta->size()/2){

                        tdc1values_1[i1] = ed.tdc->at(i1)[2];
                        tdc2values_1[i1] = ed.tdc->at(i1)[3];
                        q1s_1[i1] = QQ[2] - (QQ[0]+QQ[1])/2 ;
                        q2s_1[i1] = QQ[3] - (QQ[0]+QQ[1])/2 ;

                    } else if (q1s_1[i1-ed.ieta->size()/2] != -1) {

                        //Normalized Q1 and Q2 values for 2 channels
                        NQ1 = QQ[2] - (QQ[0]+QQ[1])/2 + q1s_1[i1-ed.ieta->size()/2];
                        NQ2 = QQ[3] - (QQ[0]+QQ[1])/2 + q2s_1[i1-ed.ieta->size()/2];

                        PL[plind].Q2overQ1_first->Fill(NQ2/NQ1);
                        PL[plind].tdc1_1vstdc1_2_first->Fill(tdc1values_1[i1-ed.ieta->size()/2], ed.tdc->at(i1)[2]);
                        PL[plind].tdc2_1vstdc2_2_first->Fill(tdc2values_1[i1-ed.ieta->size()/2],ed.tdc->at(i1)[3]);
                    }
                } else {
                    q1s_1[i1] = -1;
                }

            }

            if(e2){ //THE MAIN ONE

                if((QQ[2] > 10) && QQ[3] > 10) {
                    if(i1 < ed.ieta->size()/2){
                        PL[plind].qped1 += (QQ[0]+QQ[1])/2; //cout << plind << " birinci " << ed.ieta->at(i1) << " " << ed.iphi->at(i1) << " " << ed.depth->at(i1) << endl;
                        q1s_nocuts[i1] = QQ[2] - (QQ[0]+QQ[1])/2 ;
                        q2s_nocuts[i1] = QQ[3] - (QQ[0]+QQ[1])/2 ;
                    } else {
                        PL[plind].qped2 += (QQ[0]+QQ[1])/2; //cout << plind << " ikinci " << ed.ieta->at(i1) << " " << ed.iphi->at(i1) << " " << ed.depth->at(i1) << endl;
                        NQ1 = QQ[2] - (QQ[0]+QQ[1])/2 + q1s_nocuts[i1-ed.ieta->size()/2];
                        NQ2 = QQ[3] - (QQ[0]+QQ[1])/2 + q2s_nocuts[i1-ed.ieta->size()/2];
                        PL[plind].Q2overQ1_noCuts->Fill(NQ2/NQ1);
                    }
                }

                for(int i2=0;i2<nTS;i2++) {
                    //PL->PT is a [10][3] array whic corresponds to values at 10 TS of pulse tcd and norm
                    PL[plind].PT[i2][0]+=adc2fC_QIE10[ed.pulse->at(i1)[i2]]; //Pulse
                    PL[plind].PT[i2][1]+=ed.tdc->at(i1)[i2]; // Tdc
                    PL[plind].PT[i2][2]+=1.; //Norm
                    Qtot+=QQ[i2];//Total charge of this channel this event

                    if(i2 > 0) {
                        if(!((ed.capid->at(i1)[i2]-ed.capid->at(i1)[i2-1])==1 || (ed.capid->at(i1)[i2]-ed.capid->at(i1)[i2-1])==-3)) {
                            capidOK=false;
                        }
                    }
                }

                //Cuts for invalid TDC values:
                bool f1, f2, f3, f ;
                for( int chan = 0 ; chan < 56 ; chan++) {
                    if ((ed.ieta->at(i1)==RADDAM_CH[chan].ieta) && (ed.iphi->at(i1)==RADDAM_CH[chan].iphi)){
                        f1 = ((ed.tdc->at(i1)[3] > tdcCuts2_TS3[chan][1]) || (ed.tdc->at(i1)[3] < tdcCuts2_TS3[chan][0]));
                        f2 = ((ed.tdc->at(i1)[2] > tdcCuts2_TS2[chan][1]) || (ed.tdc->at(i1)[2] < tdcCuts2_TS2[chan][0]));
                        f3 = (QQ[2] < 10) || (QQ[3] < 10);
                        f = f1 || f2 || f3;
                        break;
                    }
                }

                if(!capidOK) PL[plind].nbadcapid+=1.;

                if(!f) {

                    if(i1 < ed.ieta->size()/2){

                        tdc1values[i1] = ed.tdc->at(i1)[2];
                        tdc2values[i1] = ed.tdc->at(i1)[3];
                        q1s[i1] = QQ[2] - (QQ[0]+QQ[1])/2 ;
                        q2s[i1] = QQ[3] - (QQ[0]+QQ[1])/2 ;
                        qtots[i1] = Qtot ;
                        qpeds[i1] = (QQ[0]+QQ[1])/2 ;

                    } else if (q1s[i1-ed.ieta->size()/2] != -1) {

                        //Normalized Q1 and Q2 values for 2 channels
                        NQ1 = QQ[2] - (QQ[0]+QQ[1])/2 + q1s[i1-ed.ieta->size()/2];
                        NQ2 = QQ[3] - (QQ[0]+QQ[1])/2 + q2s[i1-ed.ieta->size()/2];
                        NQp = (QQ[0]+QQ[1])/2 + qpeds[i1-ed.ieta->size()/2];

                        //Charge histograms
                        PL[plind].Q2overQ1->Fill(NQ2/NQ1);
                        PL[plind].Q1vsQ2->Fill(NQ1, NQ2); //Fills Q1 vs Q2 histogram
                        PL[plind].Q1->Fill(NQ1);
                        PL[plind].Q2->Fill(NQ2);
                        PL[plind].Qtot->Fill(Qtot+ qtots[i1-ed.ieta->size()/2]);
                        PL[plind].Qpedest->Fill(NQp);
                        PL[plind].Q2_Q1vsQ1->Fill(NQ2/NQ1,NQ1); //Fills Q2/Q1 vs Q1
                        PL[plind].Q2_Q1vsQ2->Fill(NQ2/NQ1,NQ2); //Fills Q2/Q1 vs Q2
                        PL[plind].Q2_Q1vsQtot->Fill(NQ2/NQ1,Qtot+ qtots[i1-ed.ieta->size()/2]); //Fills Q2/Q1 vs Qtotal
                        PL[plind].Q2_Q1vsQsum->Fill(NQ2/NQ1,NQ1+NQ2); //Fills Q2/Q1 vs Q1+Q2

                        //TDC Histograms
                        PL[plind].TDC_TS2->Fill(ed.tdc->at(i1)[2] + tdc1values[i1-ed.ieta->size()/2] );
                        PL[plind].TDC_TS3->Fill(ed.tdc->at(i1)[3] + tdc2values[i1-ed.ieta->size()/2]);
                        PL[plind].QtotvsTDC1->Fill(Qtot+ qtots[i1-ed.ieta->size()/2],ed.tdc->at(i1)[2] + tdc1values[i1-ed.ieta->size()/2]); //Fills Qtot vs TDC TS2
                        PL[plind].QtotvsTDC2->Fill(Qtot+ qtots[i1-ed.ieta->size()/2],ed.tdc->at(i1)[3] + tdc2values[i1-ed.ieta->size()/2]); //Fills Qtot vs TDC TS3
                        PL[plind].TDC1vsTDC2->Fill(ed.tdc->at(i1)[2] + tdc1values[i1-ed.ieta->size()/2] ,ed.tdc->at(i1)[3] + tdc2values[i1-ed.ieta->size()/2] ); //Fills TDC1 vs TDC2
                        PL[plind].Q1vsTDC->Fill(NQ1,ed.tdc->at(i1)[2] + tdc1values[i1-ed.ieta->size()/2]); //Fills ADC vs TDC histogram for Q1
                        PL[plind].Q2vsTDC->Fill(NQ2,ed.tdc->at(i1)[3] + tdc2values[i1-ed.ieta->size()/2]); //Fills ADC vs TDC histogram for Q2
                        PL[plind].Q2_Q1vsTDC1->Fill(NQ2/NQ1,ed.tdc->at(i1)[2]+ tdc1values[i1-ed.ieta->size()/2]); //Fills Q2/Q1 vs TDC TS2
                        PL[plind].Q2_Q1vsTDC2->Fill(NQ2/NQ1,ed.tdc->at(i1)[3]+ tdc2values[i1-ed.ieta->size()/2]); //Fills Q2/Q1 vs TDC TS3

                        //PMT Channels comparison histograms
                        PL[plind].q1_1vsq1_2->Fill(q1s[i1-ed.ieta->size()/2], QQ[2]);
                        PL[plind].q2_1vsq2_2->Fill(q2s[i1-ed.ieta->size()/2], QQ[3]);
                        PL[plind].tdc1_1vstdc1_2->Fill(tdc1values[i1-ed.ieta->size()/2], ed.tdc->at(i1)[2]);
                        PL[plind].tdc2_1vstdc2_2->Fill(tdc2values[i1-ed.ieta->size()/2],ed.tdc->at(i1)[3]);


                    }
                } else {
                    q1s[i1] = -1;
                }
            }

            if(e3){ //THE LATTER ONE

                //Cuts for invalid TDC values:
                bool f1, f2, f3, f ;
                for( int chan = 0 ; chan < 56 ; chan++) {
                    if ((ed.ieta->at(i1)==RADDAM_CH[chan].ieta) && (ed.iphi->at(i1)==RADDAM_CH[chan].iphi)){
                        f1 = ((ed.tdc->at(i1)[3] > tdcCuts3_TS3[chan][1]) || (ed.tdc->at(i1)[3] < tdcCuts3_TS3[chan][0]));
                        if(RADDAM_CH[chan].ieta == 41 && RADDAM_CH[chan].iphi == 35) {
                            f2 = (ed.tdc->at(i1)[2] > 28) || (ed.tdc->at(i1)[2] < 22);
                        } else if(RADDAM_CH[chan].ieta == 30 && RADDAM_CH[chan].iphi == 1){
                            f2 = (ed.tdc->at(i1)[2] > 26) || (ed.tdc->at(i1)[2] < 20);
                        } else {
                            f2 = false;
                        }
                        f3 = (QQ[2] < 10) || (QQ[3] < 10);
                        f = f1 || f2 || f3;
                        break;
                    }
                }

                if(!f) {
                    if(i1 < ed.ieta->size()/2){

                        tdc1values_3[i1] = ed.tdc->at(i1)[2];
                        tdc2values_3[i1] = ed.tdc->at(i1)[3];
                        q1s_3[i1] = QQ[2] - (QQ[0]+QQ[1])/2 ;
                        q2s_3[i1] = QQ[3] - (QQ[0]+QQ[1])/2 ;

                    } else if (q1s_3[i1-ed.ieta->size()/2] != -1) {

                        //Normalized Q1 and Q2 values for 2 channels
                        NQ1 = QQ[2] - (QQ[0]+QQ[1])/2 + q1s_3[i1-ed.ieta->size()/2];
                        NQ2 = QQ[3] - (QQ[0]+QQ[1])/2 + q2s_3[i1-ed.ieta->size()/2];

                        PL[plind].Q2overQ1_third->Fill(NQ2/NQ1);
                        PL[plind].tdc1_1vstdc1_2_third->Fill(tdc1values_3[i1-ed.ieta->size()/2], ed.tdc->at(i1)[2]);
                        PL[plind].tdc2_1vstdc2_2_third->Fill(tdc2values_3[i1-ed.ieta->size()/2],ed.tdc->at(i1)[3]);
                    }
                }  else {
                    q1s_3[i1] = -1;
                }

            }
        }
    }
    //END OF EVENT SELECTION AND FILLING THE HISTOGRAMS

	int nn=0;
	int iside=0;int iquad=0;
	float fnevt=((float)(ev2-ev1));

	TFile* outfile=new TFile("RadDam.root","recreate");

	for(int i2=0;i2<PL.size();i2++) {
		for(int is1=0;is1<nTS;is1++) {
			if(PL[i2].PT[is1][2]>0) {
                PL[i2].PT[is1][0]/=PL[i2].PT[is1][2]; //Pulse is divided by norm
                PL[i2].PT[is1][1]/=PL[i2].PT[is1][2]; //Tdc is divided by norm
            }
			PL[i2].p->SetBinContent(is1+1,PL[i2].PT[is1][0]); //Pulse Shapes is initiated
			PL[i2].t->SetBinContent(is1+1,PL[i2].PT[is1][1]); //TDC Shapes is initiated
		}
		PL[i2].nbadcapid/=fnevt;
	}

	{
        //HISTOGRAM FIT, DRAW AND WRITE
		string pmtnames[24]={"A2","A4","A6","A8","A1","A3","A5","A7","B2","B4","B6","B8","B1","B3","B5","B7","C2","C4","C6","C8","C1","C3","C5","C7"};
		string chnames[2]={"1+2","3+4"};
		string side[2]={"P","M"};

        int dividexaxes = 1;
        int divideyaxes = 2;
        int pagenumber = 28;

        string histonames[]= {"PulseShapes.pdf","TDCShapes.pdf","Q2_Q1vsTDC1.pdf","Q2_Q1vsTDC2.pdf","Q1vsQ2.pdf","Q1vsTDC1.pdf","Q2vsTDC2.pdf","QTotal.pdf","Q2overQ1.pdf","Q1.pdf","Q2.pdf","TDC1.pdf","TDC2.pdf","QtotvsTDC1.pdf","QtotvsTDC2.pdf", "Qpedestal.pdf","TDC1vsTDC2.pdf", "TDC1-1vsTDC1-2.pdf", "TDC2-1vsTDC2-2.pdf", "Q2_Q1vsQ1.pdf" , "Q2_Q1vsQ2.pdf" , "Q2_Q1vsQtot.pdf", "Q2_Q1vsQsum.pdf", "Q2overQ1first.pdf" , "Q2overQ1third.pdf", "TDC1comp_first.pdf", "TDC2comp_first.pdf", "TDC1comp_third.pdf", "TDC2comp_third.pdf" , "Q2overQ1_noCuts.pdf", "Q1-1vsQ1-2.pdf", "Q2-1vsQ2-2.pdf", "EventvsDelay.pdf" };

        int histonumber =  sizeof(histonames)/sizeof(histonames[0]);

        for( int n=0 ; n<histonumber ; n++){

            int l = histonames[n].length();

            char histoname1[20];
            char histoname[20] ;
            histonames[n].copy(histoname,l,0);

            for(int foo=l; foo<20;foo++){
                histoname[foo]= '\0';
            }

            TCanvas* canvas[pagenumber];

            for(int x=0;x<pagenumber;x++){
                char name[100];
                sprintf(name, "cc%d%d", n,x);
                canvas[x] = new TCanvas(name,name,5000,6000);
                gStyle->SetOptStat(0);
                gStyle->SetTitleFontSize(0.1);
                canvas[x]->Divide(dividexaxes,divideyaxes,0.001,0.001);
            }

            int currentpage = 0; int currentcanv = 1;

            for(int i1=0;i1<56;i1++) {

                    plind=-1;
                    for(int ik1=0;ik1<PL.size();ik1++)
                    {
                        MI=PL[ik1].mapind;
                        if(SEM[MI].ieta==RADDAM_CH[i1].ieta && SEM[MI].iphi==RADDAM_CH[i1].iphi)
                        {
                            plind=ik1;
                            break;
                        }
                    }

                    canvas[currentpage]->cd(currentcanv);
                    currentcanv++;
                    if(currentcanv> divideyaxes*dividexaxes ) {
                        currentpage++;
                        currentcanv = 1;
                    }
                    gPad-> SetLogy(0);
                    gStyle->SetOptStat(1111);

                    if(n==0){ //Pulse Shapes Histogram

                        PL[plind].p->Draw("hist");
                        PL[plind].p->GetYaxis()->SetRangeUser(-10.,(PL[plind].p->GetBinContent(PL[plind].p->GetMaximumBin()))*1.5);
                        PL[plind].p->Write();

                    } else if(n==1){ //Tdc histogram

                        PL[plind].t->Draw();
                        PL[plind].t->GetYaxis()->SetRangeUser(0.,75.);
                        PL[plind].t->Write();

                    } else if(n==2){ //Q2/Q1 vs TDC TS2 histogram

                        PL[plind].Q2_Q1vsTDC1->Draw("COLZ");
                        PL[plind].Q2_Q1vsTDC1->Write();

                    } else if(n==3){ //Q2/Q1 vs TDC TS3 histogram

                        PL[plind].Q2_Q1vsTDC2->Draw("COLZ");
                        PL[plind].Q2_Q1vsTDC2->Write();

                    } else if(n==4){ //Q1 vs Q2 histogram

                        PL[plind].Q1vsQ2->Draw("COLZ");
                        PL[plind].Q1vsQ2->Write();

                    } else if(n==5){ //Q1 vs TDC 1 histogram

                        PL[plind].Q1vsTDC->Draw("COLZ");
                        PL[plind].Q1vsTDC->Write();

                    } else if(n==6){ //Q2 vs TDC 2 histogram

                        PL[plind].Q2vsTDC->Draw("COLZ");
                        PL[plind].Q2vsTDC->Write();

                    } else if(n==7){ //Qtot histogram

                        gPad-> SetLogy(1);
                        PL[plind].Qtot->GetYaxis()->SetRangeUser(1,(PL[plind].Qtot->GetBinContent(PL[plind].Qtot->GetMaximumBin()))*1.5);
                        PL[plind].Qtot->Draw("hist");
                        PL[plind].Qtot->Write();

                    } else if(n==8){ //Q2/Q1 histogram

                        PL[plind].Q2overQ1->GetYaxis()->SetRangeUser(-1.,(PL[plind].Q2overQ1->GetBinContent(PL[plind].Q2overQ1->GetMaximumBin()))*1.5);
                        PL[plind].Q2overQ1->Fit("gaus");
                        if ( PL[plind].Q2overQ1->GetEntries() != 0) {
                            PL[plind].q2overq1 = PL[plind].Q2overQ1->GetFunction("gaus")->GetParameter(1);
                            PL[plind].q2overq1_error = PL[plind].Q2overQ1->GetFunction("gaus")->GetParError(1);
                            //cout << "ieta: " << RADDAM_CH[plind].ieta << " iphi: "<< RADDAM_CH[plind].iphi<< " " <<PL[plind].q2overq1 << " +-" << PL[plind].q2overq1_error << endl;
                        } else {
                            PL[plind].q2overq1 = 0;
                            PL[plind].q2overq1_error = 0;
                        }
                        gStyle->SetOptFit(0111);
                        PL[plind].Q2overQ1->Draw();
                        PL[plind].Q2overQ1->Write();

                    } else if(n==9){ //Q1 histogram

                        gPad-> SetLogy(1);
                        PL[plind].Q1->GetYaxis()->SetRangeUser(1,(PL[plind].Q1->GetBinContent(PL[plind].Q1->GetMaximumBin()))*1.5);
                        PL[plind].Q1->Draw();
                        PL[plind].Q1->Write();

                    } else if(n==10){ //Q2 histogram

                        gPad-> SetLogy(1);
                        PL[plind].Q2->GetYaxis()->SetRangeUser(1,(PL[plind].Q2->GetBinContent(PL[plind].Q2->GetMaximumBin()))*1.5);
                        PL[plind].Q2->Draw();
                        PL[plind].Q2->Write();

                    } else if(n==11){ //TDC1 histogram

                        PL[plind].TDC_TS2->Draw("hist");
                        PL[plind].TDC_TS2->Write();

                    } else if(n==12){ //TDC2 histogram

                        PL[plind].TDC_TS3->Draw("hist");
                        PL[plind].TDC_TS3->Write();

                    } else if(n==13){ //Qtot vs TDC TS2 histogram

                        PL[plind].QtotvsTDC1->Draw("COLZ");
                        PL[plind].QtotvsTDC1->Write();

                    } else if(n==14){ //Qtot vs TDC TS3 histogram
                        PL[plind].QtotvsTDC2->Draw("COLZ");
                        PL[plind].QtotvsTDC2->Write();

                    }  else if(n==15){ // Qpedestal histogram

                        gPad-> SetLogy(1);
                        PL[plind].Qpedest->GetYaxis()->SetRangeUser(1,(PL[plind].Qpedest->GetBinContent(PL[plind].Qpedest->GetMaximumBin()))*1.5);
                        PL[plind].Qpedest->Draw();
                        PL[plind].Qpedest->Write();

                    } else if(n == 16){ //TDC1 vs TDC2 Histogram

                        PL[plind].TDC1vsTDC2->Draw("COLZ");
                        PL[plind].TDC1vsTDC2->Write();

                    } else if(n == 17){ //TDC1-1 vs TDC1-2 Histogram

                        PL[plind].tdc1_1vstdc1_2->Draw("COLZ");
                        PL[plind].tdc1_1vstdc1_2->Write();

                    } else if(n == 18){ //TDC2-1 vs TDC2-2 Histogram

                        PL[plind].tdc2_1vstdc2_2->Draw("COLZ");
                        PL[plind].tdc2_1vstdc2_2->Write();

                    } else if(n == 19){ //Q2/Q1 vs Q1 histogram

                        PL[plind].Q2_Q1vsQ1->Draw("COLZ");
                        PL[plind].Q2_Q1vsQ1->Write();

                    } else if(n == 20){ //Q2/Q1 vs Q2 histogram

                        PL[plind].Q2_Q1vsQ2->Draw("COLZ");
                        PL[plind].Q2_Q1vsQ2->Write();

                    } else if(n == 21){ //Q2/Q1 vs Qtot histogram

                        PL[plind].Q2_Q1vsQtot->Draw("COLZ");
                        PL[plind].Q2_Q1vsQtot->Write();

                    } else if(n == 22){ //Q2/Q1 vs Qsum histogram

                        PL[plind].Q2_Q1vsQsum->Draw("COLZ");
                        PL[plind].Q2_Q1vsQsum->Write();

                    }  else if(n == 23){ //Q2/Q1  former

                        PL[plind].Q2overQ1_first->GetYaxis()->SetRangeUser(-1.,(PL[plind].Q2overQ1_first->GetBinContent(PL[plind].Q2overQ1_first->GetMaximumBin()))*1.5);
                        PL[plind].Q2overQ1_first->Fit("gaus");
                        if ( PL[plind].Q2overQ1_first->GetEntries() != 0) {
                            //cout << "plind: "<<plind<< " ,number of entries of 1st Q2/Q1: " << PL[plind].Q2overQ1_first->GetEntries() << endl;
                            PL[plind].q2overq1_fi = PL[plind].Q2overQ1_first->GetFunction("gaus")->GetParameter(1);//PROBLEM HERE(?)
                        } else {
                            PL[plind].q2overq1_fi = 0;
                        }
                        gStyle->SetOptFit(0111);
                        PL[plind].Q2overQ1_first->Draw();
                        PL[plind].Q2overQ1_first->Write();

                    } else if(n == 24){ //Q2/Q1 latter

                        PL[plind].Q2overQ1_third->GetYaxis()->SetRangeUser(-1.,(PL[plind].Q2overQ1_third->GetBinContent(PL[plind].Q2overQ1_third->GetMaximumBin()))*1.5);
                        PL[plind].Q2overQ1_third->Fit("gaus");
                        if ( PL[plind].Q2overQ1_third->GetEntries() != 0) {
                            PL[plind].q2overq1_th = PL[plind].Q2overQ1_third->GetFunction("gaus")->GetParameter(1);
                        } else {
                            PL[plind].q2overq1_th = 0;
                        }
                        gStyle->SetOptFit(0111);
                        PL[plind].Q2overQ1_third->Draw();
                        PL[plind].Q2overQ1_third->Write();

                    }  else if(n == 25){ //TDC1-1 vs TDC1-2 Histogram former

                        PL[plind].tdc1_1vstdc1_2_first->Draw("COLZ");
                        PL[plind].tdc1_1vstdc1_2_first->Write();

                    } else if(n == 26){ //TDC2-1 vs TDC2-2 Histogram former

                        PL[plind].tdc2_1vstdc2_2_first->Draw("COLZ");
                        PL[plind].tdc2_1vstdc2_2_first->Write();

                    }  else if(n == 27){ //TDC1-1 vs TDC1-2 Histogram latter

                        PL[plind].tdc1_1vstdc1_2_third->Draw("COLZ");
                        PL[plind].tdc1_1vstdc1_2_third->Write();

                    } else if(n == 28){ //TDC2-1 vs TDC2-2 Histogram latter

                        PL[plind].tdc2_1vstdc2_2_third->Draw("COLZ");
                        PL[plind].tdc2_1vstdc2_2_third->Write();

                    }  else if(n == 29){ // Q2/Q1 histogram with no TDC cuts

                        PL[plind].Q2overQ1_noCuts->GetYaxis()->SetRangeUser(-1.,(PL[plind].Q2overQ1_noCuts->GetBinContent(PL[plind].Q2overQ1_noCuts->GetMaximumBin()))*1.5);
                        PL[plind].Q2overQ1_noCuts->Fit("gaus");
                        if ( PL[plind].Q2overQ1_noCuts->GetEntries() != 0) {
                            PL[plind].q2overq1_nocuts = PL[plind].Q2overQ1_noCuts->GetFunction("gaus")->GetParameter(1);
                        } else {
                            PL[plind].q2overq1_nocuts = 0;
                        }
                        gStyle->SetOptFit(0111);
                        PL[plind].Q2overQ1_noCuts->Draw();
                        PL[plind].Q2overQ1_noCuts->Write();

                    } else if(n == 30){ //Q1-1 vs Q1-2 Histogram

                        PL[plind].q1_1vsq1_2->Draw("COLZ");
                        PL[plind].q1_1vsq1_2->Write();

                    } else if(n == 31){ //Q2-1 vs Q2-2 Histogram

                        PL[plind].q2_1vsq2_2->Draw("COLZ");
                        PL[plind].q2_1vsq2_2->Write();

                    } else if ( n == 32 ){
                        sprintf(hname, "(Q1+Q2)/Qt");
                        PL[plind].Q1_Q2overQtotal->SetName(hname);PL[plind].Q1_Q2overQtotal->SetTitle("");
                        PL[plind].Q1_Q2overQtotal->SetMarkerStyle(1);PL[plind].Q1_Q2overQtotal->SetMarkerColor(9);
                        PL[plind].Q1_Q2overQtotal->GetYaxis()->SetRangeUser(0,1.2);  PL[plind].Q1_Q2overQtotal->GetXaxis()->SetRangeUser(0,35000);
                        PL[plind].Q1_Q2overQtotal->GetYaxis()->SetTitle("(Q1+Q2)/Qt"); PL[plind].Q1_Q2overQtotal->GetXaxis()->SetTitle("Events");
                        PL[plind].Q1_Q2overQtotal->Draw();
                    }
            }
            for(int x=0;x<pagenumber; x++){

                char dummy[10];

                if( x==0 ) sprintf(dummy,"("); //First page of pdf file
                else if(x==pagenumber-1) sprintf(dummy,")"); //Last page of pdf file
                else sprintf(dummy,"");

                sprintf(histoname1,"%s%s",histoname,dummy);
                canvas[x]->SaveAs(histoname1);
                delete canvas[x];
            }
            sprintf(hname,"mv %s ../Plots/%d",histoname,RunNo);system(hname);

        }
        //END OF HISTOGRAM FIT, DRAW AND WRITE

        for(int ik1=0;ik1<PL.size();ik1++) {
            int max=0;float qfmax=0.;
            for(int ik2=0;ik2<32;ik2++) {
                if(PL[ik1].qf[ik2]>qfmax){
                    qfmax=PL[ik1].qf[ik2];
                    max=ik2;
                }
            }
        }

        for(int ik1=0;ik1<PL.size();ik1++) {
            PL[ik1].p=0;PL[ik1].pn=0;PL[ik1].t=0;
        }
    }

    //Q2/Q1 FILE WRITE
    char filename[100]; char command[100];

    fstream results;
    results.open("results.txt",fstream::out);

    if(!results.is_open()) {
        cout << "error"<<endl;
        return 0;
    }

    for(int i1=0;i1<56;i1++) {
            int index =-1;
            for(int ik1=0;ik1<PL.size();ik1++) {
                MI=PL[ik1].mapind;
                if(SEM[MI].ieta==RADDAM_CH[i1].ieta && SEM[MI].iphi==RADDAM_CH[i1].iphi) {
                    index=ik1;
                    break;
                }
            }
            if(i1 == 55) results << SEM[MI].ieta << " " <<SEM[MI].iphi <<" "<< PL[index].q2overq1_fi <<" "<< PL[index].q2overq1_th << " "<< PL[index].q2overq1<< " "<< PL[index].q2overq1_error<< " " << PL[index].q2overq1_nocuts << " " << PL[index].qped1<< " " << PL[index].qped2 ;
            else results << SEM[MI].ieta << " " <<SEM[MI].iphi <<" "<< PL[index].q2overq1_fi <<" "<< PL[index].q2overq1_th<< " " << PL[index].q2overq1<< " "<< PL[index].q2overq1_error<< " " << PL[index].q2overq1_nocuts << " " << PL[index].qped1<< " " << PL[index].qped2 << endl;
    }
    sprintf(command, "mv results.txt ../Data/results_%d.txt", RunNo); system(command);
    results.close();
    //END

    sprintf(hname,"mv RadDam.root ../Histos/RadDam_%d.root",RunNo);system(hname);
    outfile->Close();
}

/*
 Analyzes all runs taken so far, draws Q2/Q1 plots for each channel, saves it into Ratios.pdf
*/
void plotratios(){

    char hname[500]; char graphname[100];

    vector<channel> ChannelList;
    int numberofchannels= 56;

    //Initializes channel list array
    for(int i=0; i<numberofchannels; i++ ){
        channel ch;
        ch.index = i; // Corresponds to the channel RADDAM_CH[i]

        char title[100]; char title2[100]; char title3[100]; char title4[100];
        sprintf(title, "(%d %d)", RADDAM_CH[i].ieta, RADDAM_CH[i].iphi);
        sprintf(title2, "(%d %d %d)", RADDAM_CH[i].ieta, RADDAM_CH[i].iphi, RADDAM_CH[i].depth);
        sprintf(title3, "(%d %d %d)", RADDAM_CH[i].ieta, RADDAM_CH[i].iphi, RADDAM_CH[i].depth+2);
        sprintf(title4, "(%d)", RADDAM_CH[i].ieta);

        ch.ratioPlot = new TGraphErrors(); ch.ratioPlot->SetTitle("");
        ch.ratioPlot_1 = new TGraph(); ch.ratioPlot_1->SetTitle(title);
        ch.ratioPlot_2 = new TGraph();
        ch.ratioPlot_3 = new TGraph();
        ch.ratioCut = new TGraph(); ch.ratioCut->SetTitle(title);
        ch.ratioNoCut = new TGraph();
        ch.pedestalGraph1 = new TGraph(); ch.pedestalGraph1->SetTitle("");
        ch.pedestalGraph2 = new TGraph(); ch.pedestalGraph2->SetTitle(title3);
        ch.pedestalNormalized1 = new TGraph(); ch.pedestalNormalized1->SetTitle(title2);
        ch.pedestalNormalized2 = new TGraph(); ch.pedestalNormalized2->SetTitle(title3);
        ch.cutRelErrors = new TGraph(); ch.cutRelErrors->SetTitle(title);
        ch.dtRelErrors = new TGraph(); ch.dtRelErrors->SetTitle(title);
        ch.day0Errors = new TGraph(); ch.day0Errors->SetTitle(title);
        ch.ratioEtas = new TGraphErrors(); ch.ratioEtas->SetTitle("");

        ChannelList.push_back(ch);
    }

    TGraphErrors* colPlot_1 = new TGraphErrors(); colPlot_1->SetTitle("HF+");
    TGraphErrors* colPlot_2 = new TGraphErrors();TGraphErrors* colPlot_14 = new TGraphErrors();
    TGraphErrors* colPlot_3 = new TGraphErrors();TGraphErrors* colPlot_15 = new TGraphErrors();
    TGraphErrors* colPlot_4 = new TGraphErrors();TGraphErrors* colPlot_16 = new TGraphErrors();
    TGraphErrors* colPlot_5 = new TGraphErrors();TGraphErrors* colPlot_17 = new TGraphErrors();
    TGraphErrors* colPlot_6 = new TGraphErrors();TGraphErrors* colPlot_18 = new TGraphErrors();
    TGraphErrors* colPlot_7 = new TGraphErrors();TGraphErrors* colPlot_19 = new TGraphErrors();
    TGraphErrors* colPlot_8 = new TGraphErrors();TGraphErrors* colPlot_20 = new TGraphErrors();
    TGraphErrors* colPlot_9 = new TGraphErrors();TGraphErrors* colPlot_21 = new TGraphErrors();
    TGraphErrors* colPlot_10 = new TGraphErrors();TGraphErrors* colPlot_22 = new TGraphErrors();
    TGraphErrors* colPlot_11 = new TGraphErrors();TGraphErrors* colPlot_23 = new TGraphErrors();
    TGraphErrors* colPlot_12 = new TGraphErrors();TGraphErrors* colPlot_24 = new TGraphErrors();
    TGraphErrors* colPlot_13 = new TGraphErrors();

    //Opens run list file
    ifstream runlist("RunList.csv");
    if (!runlist.is_open()){
        cout << "Run list cannot be openned!" << endl;
        return;
    }

    int runno, runname; float datetime;
    string a, b, c, d, e, f; //Consecutively: Run no(1,2,3,...) , run number(321446,...), date, time(hour), date+time, date and time converted to a number

    //Reads the run list and initializes run data, stores in allRuns code
    while (runlist.good()) {

        getline(runlist,a, ','); getline(runlist, b,','); getline(runlist, c,',');
        getline(runlist, d,','); getline(runlist, e,','); getline(runlist, f,'\n');
        if (a.empty()) break;

        runno = stoi(a); runname =stoi(b); datetime=stoi(f);

        singleRun.runno = runno;
        singleRun.runname = runname;
        singleRun.datetime = datetime;

        allRuns.push_back(singleRun);
    }
    runlist.close();

    int runs = allRuns.size();
    float refPdstl[numberofchannels][2] = {{0}};//Use averages of all pedestals instead of first one
    float pdstl[runs][numberofchannels][2];


    //Pedestal values of all channels, normalized to the average of all, in a single graph
    TGraph* normPed = new TGraph();
    normPed->SetName("Normalized Pedestals of All Channels");normPed->SetTitle("Normalized Pedestals of All Channels");

    //Pedestal values of all channels in a single graph
    TGraph* pedestalAll = new TGraph();
    pedestalAll->SetName("Pedestals of All Channels");pedestalAll->SetTitle("Pedestals of All Channels");

    //Traces all runs, opens the data file, takes channel info and ratio value for this run
    for(int i=0; i<runs; i++) {

        char filelocation[100];
        sprintf(filelocation, "../Data/results_%d.txt",allRuns[i].runname);
        ifstream run(filelocation);
        if(run.is_open()){

            int ieta, iphi; float ratio1, ratio2, ratio3, rat2error, ratiowithnocuts, pedes1, pedes2;

            while(!run.eof()){
                run >> ieta >> iphi >> ratio1 >> ratio3 >> ratio2 >> rat2error >> ratiowithnocuts >> pedes1 >> pedes2;

                pedes1 /= 5000; pedes2 /=5000;

                //Finds the correct channel
                for(int j=0;j<numberofchannels;j++) {

                    if (ieta == RADDAM_CH[j].ieta && iphi == RADDAM_CH[j].iphi) {

                        runProperties r;
                        r.ratioValue = ratio2;
                        r.ratioValue1 = ratio1;
                        r.ratioValue3 = ratio3;
                        r.ratioError = rat2error;
                        r.ratioNoCuts = ratiowithnocuts;
                        r.pedestalValue1 = pedes1;
                        r.pedestalValue2 = pedes2;
                        r.datetime = allRuns[i].datetime;

                        ChannelList[j].runForOneChannel.push_back(r);

                        refPdstl[j][0] += pedes1; refPdstl[j][1] += pedes2;
                        pdstl[i][j][0] = pedes1; pdstl[i][j][1] = pedes2;

                        pedestalAll->SetPoint(pedestalAll->GetN(), allRuns[i].datetime, pedes1);
                        pedestalAll->SetPoint(pedestalAll->GetN(), allRuns[i].datetime, pedes2);
                    }
                }
            }
            run.close();
        } else {
            cout << "Run number " << allRuns[i].runname << " file does not exists!" << endl;
            break;
        }
    }

    for(int m=0; m < runs ; m++){
        for(int n=0; n<numberofchannels; n++){
            if( m == 0 ){
                refPdstl[n][0] /= runs; refPdstl[n][1] /= runs;
            }
            float l1 = pdstl[m][n][0]/refPdstl[n][0]; float l2 = pdstl[m][n][1]/refPdstl[n][1];

            normPed->SetPoint(normPed->GetN(),allRuns[m].datetime, l1);
            normPed->SetPoint(normPed->GetN(),allRuns[m].datetime, l2);
        }
    }

    //INTEGRATED LUMI
    ifstream lumiday("lumiday.csv");
    if (!lumiday.is_open()){
        cout << "Lumi day by day cannot be openned!" << endl;
        return;
    }

    float lumi[105], lumi2[105], days[105]; int ind=0; float integratedlumi = 0;
    string g, h; //Consecutively: Day number , luminosity

    while (lumiday.good()) {

        getline(lumiday,g, ','); getline(lumiday, h,'\n');
        if (g.empty()) break;

        lumi[ind] = stoi(h)*0.0000000005 +0.8; days[ind] =stoi(g);
        lumi2[ind] = stoi(h)*0.000001;
        integratedlumi += stoi(h);
        ind++;
    }
    lumiday.close();
    //INTEGRATED LUMI END

    //Creates the graph and canvas to display it
    TFile* outfile=new TFile("Ratios.root","recreate");

    TGraph* lumin = new TGraph(ind,days, lumi2); lumin->SetFillColor(46);
    lumin->SetTitle("Day by Day Integrated Luminosity");
    lumin->GetXaxis()->SetTitle("Days"); lumin->GetYaxis()->SetTitle("Luminosity (1/pb)");
    TCanvas* ccc = new TCanvas();
    lumin->Draw("AB");
    ccc->SaveAs("Luminosity.pdf");
    delete ccc;
    cout << "integrated lumin for data taking: " << integratedlumi << endl;

    //ALL PEDESTALS GRAPH
    TCanvas* canv = new TCanvas("page","page",5000,6000);
    canv->Divide(1,2,0.001,0.001);

    canv->cd(1);
    pedestalAll->SetMarkerStyle(21); pedestalAll->SetMarkerColor(3);
    pedestalAll->GetYaxis()->SetRangeUser(-5,5);  //pedestalAll->GetXaxis()->SetRangeUser(0,150);
    pedestalAll->GetYaxis()->SetTitle("Pedestals (fC)"); pedestalAll->GetXaxis()->SetTitle("Days");
    pedestalAll->Draw("AP");

    canv->cd(2);
    normPed->GetYaxis()->SetRangeUser(0.5,1.5);  normPed->GetXaxis()->SetRangeUser(0,150);
    normPed->SetMarkerStyle(19); normPed->SetMarkerColor(4);
    normPed->GetYaxis()->SetTitle("Normalized Pedestals"); normPed->GetXaxis()->SetTitle("Days");
    normPed->Draw("AP");

    canv->SaveAs("AllPedestals.pdf");
    delete canv;
    //END OF ALL PEDESTALS GRAPH

    //RATIOS HISTOGRAMS DRAW
    int dividexaxes = 1;
    int divideyaxes = 2;
    int pagenumber = 28;
    string histonames[]= {"Ratios.pdf","Ratios_DT.pdf","Ratios_TDCCuts.pdf","Pedestals1.pdf", /*"Pedestals2.pdf", "NormPedestals1.pdf", "NormPedestals2.pdf",*/ "TDCErrors.pdf", "DTErrors.pdf"};

    int histonumber =  sizeof(histonames)/sizeof(histonames[0]);
    float reference1[ChannelList.size()] = {0}; float reference2[ChannelList.size()] = {0};
    float etaAverage = 0, count; float arrayHFP[7]; float arrayHFM[7]; float arrayHF[7];
    Double_t max_min[56]; Double_t max_minError[56]; Double_t chans[56]; float wa_minmax, wa_minmax_denom = 0;

    fstream onyuzmilyonbaloncuk;
    onyuzmilyonbaloncuk.open("RatiosExcel.csv",fstream::out);

    if(!onyuzmilyonbaloncuk.is_open()) {
        cout << "error in raddam csv" <<endl;
    }
    onyuzmilyonbaloncuk << "iEta;iPhi;Day;Ratio;Ratio Error;Luminosity(fb-1)" <<endl;

    for( int n=0 ; n<histonumber ; n++){

        int l = histonames[n].length();

        char histoname1[20]; char histoname[20] ;
        histonames[n].copy(histoname,l,0);

        for(int foo=l; foo<20;foo++){
            histoname[foo]= '\0';
        }

        TCanvas* canvas[pagenumber];

        for(int x=0;x<pagenumber;x++){
            char name[100];
            sprintf(name, "page%d%d",x,n);
            canvas[x] = new TCanvas(name,name,5000,6000);
            gStyle->SetOptStat(0);
            gStyle->SetTitleFontSize(0.1);
            canvas[x]->Divide(dividexaxes,divideyaxes,0.001,0.001);
        }
        int currentpage = 0; int currentcanv = 1;

        int numberofruns = ChannelList[0].runForOneChannel.size();
        float errAv = 0;//Average statistical error
        float iEta30[8][2][numberofruns];

        for(int k = 0; k<ChannelList.size(); k++){

            canvas[currentpage]->cd(currentcanv);
            currentcanv++;
            if(currentcanv> divideyaxes*dividexaxes ) {
                currentpage++;
                currentcanv = 1;
            }
            if(n==0) { //Ratios Graph


                float A = ChannelList[k].runForOneChannel[0].ratioValue;
                float DA = ChannelList[k].runForOneChannel[0].ratioError;

                float errAv2 = 0;
                for(int l =0; l<numberofruns; l++){

                    float B = ChannelList[k].runForOneChannel[l].ratioValue;
                    float DB = ChannelList[k].runForOneChannel[l].ratioError;

                    float day = ChannelList[k].runForOneChannel[l].datetime;

                    float rat = B/A, er;
                    er = (B/A)*sqrt((DA/A)*(DA/A) + (DB/B)*(DB/B));

                    if((l >= 0 && l < 10 && l != 4 && l!=1) || l == 11) {
                        onyuzmilyonbaloncuk << RADDAM_CH[k].ieta << ";" << RADDAM_CH[k].iphi << ";" << day << ";" << rat <<";"<<er << "; " << "0"<<endl;
                        /*if(RADDAM_CH[k].ieta == 30 && RADDAM_CH[k].iphi == 1){
                            iEta30[0][0][l] = rat;
                            iEta30[0][1][l] = er;
                        } else if(RADDAM_CH[k].ieta == 30 && RADDAM_CH[k].iphi == 21){
                            iEta30[1][0][l] = rat;
                            iEta30[1][1][l] = er;
                        } else if(RADDAM_CH[k].ieta == 30 && RADDAM_CH[k].iphi == 37){
                            iEta30[2][0][l] = rat;
                            iEta30[2][1][l] = er;
                        } else if(RADDAM_CH[k].ieta == 30 && RADDAM_CH[k].iphi == 57){
                            iEta30[3][0][l] = B/A;
                            iEta30[3][1][l] = (B/A)*sqrt((DA/A)*(DA/A) + (DB/B)*(DB/B));
                        } else if(RADDAM_CH[k].ieta == -30 && RADDAM_CH[k].iphi == 15){
                            iEta30[4][0][l] = B/A;
                            iEta30[4][1][l] = (B/A)*sqrt((DA/A)*(DA/A) + (DB/B)*(DB/B));
                        } else if(RADDAM_CH[k].ieta == -30 && RADDAM_CH[k].iphi == 35){
                            iEta30[5][0][l] = B/A;
                            iEta30[5][1][l] = (B/A)*sqrt((DA/A)*(DA/A) + (DB/B)*(DB/B));
                        } else if(RADDAM_CH[k].ieta == -30 && RADDAM_CH[k].iphi == 51){
                            iEta30[6][0][l] = B/A;
                            iEta30[6][1][l] = (B/A)*sqrt((DA/A)*(DA/A) + (DB/B)*(DB/B));
                        } else if(RADDAM_CH[k].ieta == -30 && RADDAM_CH[k].iphi == 71){
                            iEta30[7][0][l] = B/A;
                            iEta30[7][1][l] = (B/A)*sqrt((DA/A)*(DA/A) + (DB/B)*(DB/B));
                        } else if(RADDAM_CH[k].ieta > 0  && (RADDAM_CH[k].iphi == 1 || RADDAM_CH[k].iphi == 71)){


                            //er = (rat/iEta30[0][0][l])*sqrt((iEta30[0][1][l]/iEta30[0][0][l])*(iEta30[0][1][l]/iEta30[0][0][l]) + (er/rat)*(er/rat));
                            rat /= iEta30[0][0][l];
                            //cout << "Channel ieta=" << RADDAM_CH[k].ieta << ", iphi=" << RADDAM_CH[k].iphi << " ratio over ratio= " << rat <<", error=  " << er << endl;

                            if(RADDAM_CH[k].ieta == 32) {
                                int poni = colPlot_1->GetN();
                                colPlot_1->SetPoint(poni,day,rat);
                                colPlot_1->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 34) {
                                int poni = colPlot_2->GetN();
                                colPlot_2->SetPoint(poni,day,rat);
                                colPlot_2->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 36) {
                                int poni = colPlot_3->GetN();
                                colPlot_3->SetPoint(poni,day,rat);
                                colPlot_3->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 38) {
                                int poni = colPlot_4->GetN();
                                colPlot_4->SetPoint(poni,day,rat);
                                colPlot_4->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 40) {
                                int poni = colPlot_5->GetN();
                                colPlot_5->SetPoint(poni,day,rat);
                                colPlot_5->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 41) {
                                int poni = colPlot_6->GetN();
                                colPlot_6->SetPoint(poni,day,rat);
                                colPlot_6->SetPointError(poni,0,er);
                            }
                        } else if(RADDAM_CH[k].ieta > 0  && (RADDAM_CH[k].iphi == 21 || RADDAM_CH[k].iphi == 19)){

                            //er = (rat/iEta30[1][0][l])*sqrt((iEta30[1][1][l]/iEta30[1][0][l])*(iEta30[1][1][l]/iEta30[1][0][l]) + (er/rat)*(er/rat));
                            rat /= iEta30[1][0][l];

                            if(RADDAM_CH[k].ieta == 32) {
                                int poni = colPlot_7->GetN();
                                colPlot_7->SetPoint(poni,day,rat);
                                colPlot_7->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 34) {
                                int poni = colPlot_8->GetN();
                                colPlot_8->SetPoint(poni,day,rat);
                                colPlot_8->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 36) {
                                int poni = colPlot_9->GetN();
                                colPlot_9->SetPoint(poni,day,rat);
                                colPlot_9->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 38) {
                                int poni = colPlot_10->GetN();
                                colPlot_10->SetPoint(poni,day,rat);
                                colPlot_10->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 40) {
                                int poni = colPlot_11->GetN();
                                colPlot_11->SetPoint(poni,day,rat);
                                colPlot_11->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 41) {
                                int poni = colPlot_12->GetN();
                                colPlot_12->SetPoint(poni,day,rat);
                                colPlot_12->SetPointError(poni,0,er);
                            }
                        } else if(RADDAM_CH[k].ieta > 0  && (RADDAM_CH[k].iphi == 37 || RADDAM_CH[k].iphi == 35)){

                            //er = (rat/iEta30[2][0][l])*sqrt((iEta30[2][1][l]/iEta30[2][0][l])*(iEta30[2][1][l]/iEta30[2][0][l]) + (er/rat)*(er/rat));
                            rat /= iEta30[2][0][l];

                            if(RADDAM_CH[k].ieta == 32) {
                                int poni = colPlot_13->GetN();
                                colPlot_13->SetPoint(poni,day,rat);
                                colPlot_13->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 34) {
                                int poni = colPlot_14->GetN();
                                colPlot_14->SetPoint(poni,day,rat);
                                colPlot_14->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 36) {
                                int poni = colPlot_15->GetN();
                                colPlot_15->SetPoint(poni,day,rat);
                                colPlot_15->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 38) {
                                int poni = colPlot_16->GetN();
                                colPlot_16->SetPoint(poni,day,rat);
                                colPlot_16->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 40) {
                                int poni = colPlot_17->GetN();
                                colPlot_17->SetPoint(poni,day,rat);
                                colPlot_17->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 41) {
                                int poni = colPlot_18->GetN();
                                colPlot_18->SetPoint(poni,day,rat);
                                colPlot_18->SetPointError(poni,0,er);
                            }
                        } else if(RADDAM_CH[k].ieta > 0  && (RADDAM_CH[k].iphi == 57 || RADDAM_CH[k].iphi == 55)){

                            //er = (rat/iEta30[3][0][l])*sqrt((iEta30[3][1][l]/iEta30[3][0][l])*(iEta30[3][1][l]/iEta30[3][0][l]) + (er/rat)*(er/rat));
                            rat /= iEta30[3][0][l];

                            if(RADDAM_CH[k].ieta == 32) {
                                int poni = colPlot_19->GetN();
                                colPlot_19->SetPoint(poni,day,rat);
                                colPlot_19->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 34) {
                                int poni = colPlot_20->GetN();
                                colPlot_20->SetPoint(poni,day,rat);
                                colPlot_20->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 36) {
                                int poni = colPlot_21->GetN();
                                colPlot_21->SetPoint(poni,day,rat);
                                colPlot_21->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 38) {
                                int poni = colPlot_22->GetN();
                                colPlot_22->SetPoint(poni,day,rat);
                                colPlot_22->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 40) {
                                int poni = colPlot_23->GetN();
                                colPlot_23->SetPoint(poni,day,rat);
                                colPlot_23->SetPointError(poni,0,er);
                            } else if (RADDAM_CH[k].ieta == 41) {
                                int poni = colPlot_24->GetN();
                                colPlot_24->SetPoint(poni,day,rat);
                                colPlot_24->SetPointError(poni,0,er);
                            }
                        }*/
                    }

                    int pointi = ChannelList[k].ratioPlot->GetN();
                    if ( k != 41)
                        errAv += er;
                        errAv2 += er;
                    ChannelList[k].ratioPlot->SetPoint(pointi,ChannelList[k].runForOneChannel[l].datetime,B/A);
                    ChannelList[k].ratioPlot->SetPointError(pointi,0,er);

                    //DENEYSEL
                    //Tek bir kanal, tek bir run
                    if( l == 0 ) {
                        float norm = 0, mean = 0, var = 0, sdev, nerror;
                        for(int dene = 0; dene < numberofruns; dene++){
                            //Her bir day 0 icin normalize ediyoruz
                            float C = ChannelList[k].runForOneChannel[dene].ratioValue;
                            mean += B/C;
                            norm += 1;
                        }
                        mean /= norm;

                        for(int dene= 0; dene < numberofruns; dene++){
                            //Her bir day 0 icin normalize ediyoruz
                            float C = ChannelList[k].runForOneChannel[dene].ratioValue;
                            var += (B/C-mean)*(B/C-mean);
                        }
                        var /= (norm-1);
                        sdev = sqrt(var);
                        nerror = sdev/mean;

                        if(k != 41){
                            etaAverage += nerror; count++;
                            //cout << "CHANNEL: " << k << " error: " << nerror << " eta av: " << etaAverage << " count: " << count<<  endl;
                        }

                        if( k % 4 == 3 ){
                            etaAverage /= count;
                            //cout << " eta average: " << etaAverage << endl;
                            if(k < 28)
                                arrayHFP[k/4] = etaAverage;
                            else
                                arrayHF[k/4-7] = (arrayHFP[k/4-7]+etaAverage)/2;
                                arrayHFM[k/4-7] = etaAverage;
                            etaAverage = 0; count = 0;
                        }
                        //ChannelList[k].day0Errors->SetPoint(ChannelList[k].day0Errors->GetN(),ChannelList[k].runForOneChannel[l].datetime,nerror);
                    }
                    //END OF DENEYSEL
                }

                //cout << errAv2/numberofruns << endl;
                ChannelList[k].ratioPlot->GetYaxis()->SetRangeUser(0.8,1.2);  ChannelList[k].ratioPlot->GetXaxis()->SetRangeUser(0,150);
                ChannelList[k].ratioPlot->GetYaxis()->SetTitle("Normalized Ratio (Q2/Q1)"); ChannelList[k].ratioPlot->GetXaxis()->SetTitle("Days");
                ChannelList[k].ratioPlot->SetMarkerStyle(8); ChannelList[k].ratioPlot->SetMarkerSize(5);
                ChannelList[k].ratioPlot->Draw("A*"); ChannelList[k].ratioPlot->Write();

                lumin->SetMarkerStyle(1); lumin->SetMarkerColor(9); lumin->SetFillColorAlpha(9,0.35); lumin->Draw("PB");

                //MAX AND MIN VALUE OPERATIONS
                float maxValue = TMath::MaxElement(numberofruns,ChannelList[k].ratioPlot->GetY());
                float minValue = TMath::MinElement(numberofruns,ChannelList[k].ratioPlot->GetY());
                float maxindex = TMath::LocMax(numberofruns,ChannelList[k].ratioPlot->GetY());
                float minindex = TMath::LocMin(numberofruns,ChannelList[k].ratioPlot->GetY());
                float maxerror = ChannelList[k].ratioPlot->GetErrorY(maxindex);
                float minerror = ChannelList[k].ratioPlot->GetErrorY(minindex);
                //cout << "Channel: " << k << " max=" << maxValue << " & min=" << minValue << endl;
                max_min[k] = (1-minValue)/(maxValue-minValue);
                //max_minError[k] = max_min[k]*sqrt((minerror/(1-minValue))*(minerror/(1-minValue)) + ((maxerror*maxerror+minerror*minerror)/((maxValue-minValue)*(maxValue-minValue))));
                max_minError[k] = (sqrt((1-maxValue)*(1-maxValue)*minerror*minerror + (1-minValue)*(1-minValue)*maxerror*maxerror))/((maxValue-minValue)*(maxValue-minValue));
                chans[k] = k;
                float dum = max_minError[k];
                if(!isnan(dum)){
                    //cout << k << "oran = " << max_min[k] << " ,hata = " << max_minError[k] << endl;
                    float weight = 1/(max_minError[k]*max_minError[k]);
                    wa_minmax += max_min[k]*weight;
                    wa_minmax_denom += weight;
                }

            } else if(n==1) { // Ratios of former, main and latter delay time settings, superimposed

                float percentageSTD[numberofruns][2]; float average, std, waverage;

                float A1 = ChannelList[k].runForOneChannel[0].ratioValue1;
                float A2 = ChannelList[k].runForOneChannel[0].ratioValue;
                float A3 = ChannelList[k].runForOneChannel[0].ratioValue3;

                for(int l =0; l<numberofruns; l++){

                    float B1 = ChannelList[k].runForOneChannel[l].ratioValue1;
                    float B2 = ChannelList[k].runForOneChannel[l].ratioValue;
                    float B3 = ChannelList[k].runForOneChannel[l].ratioValue3;

                    ChannelList[k].ratioPlot_1->SetPoint(ChannelList[k].ratioPlot_1->GetN(), ChannelList[k].runForOneChannel[l].datetime, B1/A1);
                    ChannelList[k].ratioPlot_2->SetPoint(ChannelList[k].ratioPlot_2->GetN(), ChannelList[k].runForOneChannel[l].datetime, B2/A2);
                    ChannelList[k].ratioPlot_3->SetPoint(ChannelList[k].ratioPlot_3->GetN(), ChannelList[k].runForOneChannel[l].datetime, B3/A3);

                    average = (B1/A1 + B2/A2 + B3/A3)/3;
                    std = sqrt(((B1/A1-average)*(B1/A1-average) + (B2/A2-average)*(B2/A2-average) + (B3/A3-average)*(B3/A3-average))/2);
                    percentageSTD[l][0] = std/average;
                    percentageSTD[l][1] = 1/(std*std);

                    ChannelList[k].dtRelErrors->SetPoint(ChannelList[k].dtRelErrors->GetN(),ChannelList[k].runForOneChannel[l].datetime, std/average);
                    //cout << "Channel " << RADDAM_CH[k].ieta << " " << RADDAM_CH[k].iphi << " average: " << average << " std: " << std << endl;

                }

                ChannelList[k].ratioPlot_1->SetMarkerColor(7);
                ChannelList[k].ratioPlot_1->GetYaxis()->SetRangeUser(0.8,1.2);  ChannelList[k].ratioPlot_1->GetXaxis()->SetRangeUser(0,150);
                ChannelList[k].ratioPlot_1->GetYaxis()->SetTitle("Normalized Ratio (Q2/Q1)"); ChannelList[k].ratioPlot_1->GetXaxis()->SetTitle("Days");
                ChannelList[k].ratioPlot_1->Draw("A*");

                ChannelList[k].ratioPlot->SetLineWidth(3);
                ChannelList[k].ratioPlot_2->GetYaxis()->SetRangeUser(0.8,1.2);  ChannelList[k].ratioPlot_2->GetXaxis()->SetRangeUser(0,150);
                ChannelList[k].ratioPlot_2->SetMarkerStyle(21); ChannelList[k].ratioPlot_2->SetMarkerColor(6);
                ChannelList[k].ratioPlot_2->Draw("P");

                ChannelList[k].ratioPlot_3->SetLineWidth(3);
                ChannelList[k].ratioPlot_3->GetYaxis()->SetRangeUser(0.8,1.2);  ChannelList[k].ratioPlot_3->GetXaxis()->SetRangeUser(0,150);
                ChannelList[k].ratioPlot_3->SetMarkerStyle(10); ChannelList[k].ratioPlot_3->SetMarkerColor(4);
                ChannelList[k].ratioPlot_3->Draw("P");

                float wi = 0, miwi = 0;
                for(int l =1; l<numberofruns; l++){
                    if(percentageSTD[l][0] != 0) {
                        miwi += percentageSTD[l][0]*percentageSTD[l][1];
                        wi += percentageSTD[l][1];
                    }
                }
                waverage = miwi/wi;
                ChannelList[k].sysDTerror = waverage;
                //cout << "Channel " << RADDAM_CH[k].ieta << " " << RADDAM_CH[k].iphi << " systematic delay time error " << waverage << endl;

            } else if(n==2) { // Ratios of with tdc cuts and without, superimposed

                float percentageSTD[numberofruns][2]; float average, std, waverage;

                float A1 = ChannelList[k].runForOneChannel[0].ratioValue;
                float A2 = ChannelList[k].runForOneChannel[0].ratioNoCuts;

                for(int l =0; l<numberofruns; l++){

                    float B1 = ChannelList[k].runForOneChannel[l].ratioValue;
                    float B2 = ChannelList[k].runForOneChannel[l].ratioNoCuts;

                    ChannelList[k].ratioCut->SetPoint(ChannelList[k].ratioCut->GetN(), ChannelList[k].runForOneChannel[l].datetime, B1/A1);
                    ChannelList[k].ratioNoCut->SetPoint(ChannelList[k].ratioNoCut->GetN(), ChannelList[k].runForOneChannel[l].datetime, B2/A2);

                    average = (B1/A1 + B2/A2)/2;
                    std = sqrt((B1/A1-average)*(B1/A1-average) + (B2/A2-average)*(B2/A2-average));
                    percentageSTD[l][0] = std/average;
                    percentageSTD[l][1] = 1/(std*std);

                    ChannelList[k].cutRelErrors->SetPoint(ChannelList[k].cutRelErrors->GetN(),ChannelList[k].runForOneChannel[l].datetime, std/average);
                    //cout << "Channel " << RADDAM_CH[k].ieta << " " << RADDAM_CH[k].iphi << " average: " << average << " std: " << std << endl;

                }

                ChannelList[k].ratioCut->SetMarkerColor(8);
                ChannelList[k].ratioCut->GetYaxis()->SetRangeUser(0.8,1.2);  ChannelList[k].ratioCut->GetXaxis()->SetRangeUser(0,150);
                ChannelList[k].ratioCut->GetYaxis()->SetTitle("Normalized Ratio (Q2/Q1)"); ChannelList[k].ratioCut->GetXaxis()->SetTitle("Days");
                ChannelList[k].ratioCut->Draw("A*");

                ChannelList[k].ratioNoCut->SetLineWidth(3);
                ChannelList[k].ratioNoCut->GetYaxis()->SetRangeUser(0.8,1.2);  ChannelList[k].ratioNoCut->GetXaxis()->SetRangeUser(0,150);
                ChannelList[k].ratioNoCut->SetMarkerStyle(21); ChannelList[k].ratioNoCut->SetMarkerColor(9);
                ChannelList[k].ratioNoCut->Draw("P");

                float wi = 0, miwi = 0;
                for(int l =1; l<numberofruns; l++){
                    if(percentageSTD[l][0] != 0) {
                        miwi += percentageSTD[l][0]*percentageSTD[l][1];
                        wi += percentageSTD[l][1];
                    }
                }
                waverage = miwi/wi;
                ChannelList[k].sysTDCerror = waverage;
                //cout << "Channel " << RADDAM_CH[k].ieta << " " << RADDAM_CH[k].iphi << " systematic TDC error " << waverage << endl;

            } else if(n==3) { //Pedestal values of depth <= 2 channels

                for(int l =0; l<numberofruns; l++){
                    ChannelList[k].pedestalGraph1->SetPoint(ChannelList[k].pedestalGraph1->GetN(),ChannelList[k].runForOneChannel[l].datetime, ChannelList[k].runForOneChannel[l].pedestalValue1);
                    reference1[k] += ChannelList[k].runForOneChannel[l].pedestalValue1;
                }

                gStyle->SetOptStat(1101);
                ChannelList[k].pedestalGraph1->GetYaxis()->SetRangeUser(-4,4);  ChannelList[k].pedestalGraph1->GetXaxis()->SetRangeUser(0,150);
                ChannelList[k].pedestalGraph1->GetYaxis()->SetTitle("Pedestals (fC)"); ChannelList[k].pedestalGraph1->GetXaxis()->SetTitle("Days");
                ChannelList[k].pedestalGraph1->SetMarkerStyle(3); ChannelList[k].pedestalGraph1->SetMarkerSize(5);
                ChannelList[k].pedestalGraph1->Draw("A*");
                ChannelList[k].pedestalGraph1->Write();

            } /*else if(n==4) { //Pedestal values of depth > 2 channels

                for(int l =0; l<numberofruns; l++){
                    ChannelList[k].pedestalGraph2->SetPoint(ChannelList[k].pedestalGraph2->GetN(),ChannelList[k].runForOneChannel[l].datetime, ChannelList[k].runForOneChannel[l].pedestalValue2);
                    reference2[k] += ChannelList[k].runForOneChannel[l].pedestalValue2;
                }

                gStyle->SetOptStat(1101);
                ChannelList[k].pedestalGraph2->GetYaxis()->SetRangeUser(-4,4);  ChannelList[k].pedestalGraph2->GetXaxis()->SetRangeUser(0,150);
                ChannelList[k].pedestalGraph2->GetYaxis()->SetTitle("Pedestals (fC)"); ChannelList[k].pedestalGraph2->GetXaxis()->SetTitle("Days");
                ChannelList[k].pedestalGraph2->Draw("A*");
                ChannelList[k].pedestalGraph2->Write();

            } else if (n == 5){ //Normalized pedestals of depth <= 2 channels

                reference1[k] /= numberofruns;

                for(int l =0; l<numberofruns; l++){
                    ChannelList[k].pedestalNormalized1->SetPoint(ChannelList[k].pedestalNormalized1->GetN(),ChannelList[k].runForOneChannel[l].datetime, (ChannelList[k].runForOneChannel[l].pedestalValue1)/reference1[k]);
                }

                ChannelList[k].pedestalNormalized1->GetYaxis()->SetRangeUser(0,2);  ChannelList[k].pedestalNormalized1->GetXaxis()->SetRangeUser(0,150);
                ChannelList[k].pedestalNormalized1->GetYaxis()->SetTitle("Normalized Pedestals"); ChannelList[k].pedestalNormalized1->GetXaxis()->SetTitle("Days");
                ChannelList[k].pedestalNormalized1->Draw("A*");
                ChannelList[k].pedestalNormalized1->Write();

            } else if (n == 6) { //Normalized pedestals of depth > 2 channels
                reference2[k] /= numberofruns;

                for(int l =0; l<numberofruns; l++){
                    ChannelList[k].pedestalNormalized2->SetPoint(ChannelList[k].pedestalNormalized2->GetN(),ChannelList[k].runForOneChannel[l].datetime, (ChannelList[k].runForOneChannel[l].pedestalValue2)/reference2[k]);
                }

                ChannelList[k].pedestalNormalized2->GetYaxis()->SetRangeUser(0,2);  ChannelList[k].pedestalNormalized2->GetXaxis()->SetRangeUser(0,150);
                ChannelList[k].pedestalNormalized2->GetYaxis()->SetTitle("Normalized Pedestals"); ChannelList[k].pedestalNormalized2->GetXaxis()->SetTitle("Days");
                ChannelList[k].pedestalNormalized2->Draw("A*");
                ChannelList[k].pedestalNormalized2->Write();
            } */ else if(n == 4) {

                gStyle->SetOptStat(1101);
                ChannelList[k].cutRelErrors->SetMarkerColor(4);
                ChannelList[k].cutRelErrors->GetYaxis()->SetRangeUser(0.,0.04);  ChannelList[k].cutRelErrors->GetXaxis()->SetRangeUser(0,150);
                ChannelList[k].cutRelErrors->Draw("A*");

            } else if (n==5){

                gStyle->SetOptStat(1101);
                ChannelList[k].dtRelErrors->SetMarkerColor(8);
                ChannelList[k].dtRelErrors->GetYaxis()->SetRangeUser(0.,0.04);  ChannelList[k].dtRelErrors->GetXaxis()->SetRangeUser(0,150);
                ChannelList[k].dtRelErrors->Draw("A*");

            } else if (n==6){

                gStyle->SetOptStat(1101);
                ChannelList[k].day0Errors->SetMarkerColor(7);
                ChannelList[k].day0Errors->GetYaxis()->SetRangeUser(0.,0.04);  ChannelList[k].day0Errors->GetXaxis()->SetRangeUser(0,150);
                ChannelList[k].day0Errors->Draw("A*");

            } else if (n==7){

                float A = ChannelList[k].runForOneChannel[0].ratioValue;
                float DA = ChannelList[k].runForOneChannel[0].ratioError;

                for(int l =0; l<numberofruns; l++){

                    float B = ChannelList[k].runForOneChannel[l].ratioValue;
                    float DB = ChannelList[k].runForOneChannel[l].ratioError;

                }
            }
        }
        onyuzmilyonbaloncuk.close();

        errAv /= ChannelList[0].runForOneChannel.size()*55;
        //cout << "Error = " << errAv << endl;

        for(int x=0;x<pagenumber; x++){

            char dummy[10];

            if( x==0 ) sprintf(dummy,"("); //First page of pdf file
            else if(x==pagenumber-1) sprintf(dummy,")"); //Last page of pdf file
            else sprintf(dummy,"");

            sprintf(histoname1,"%s%s",histoname,dummy);
            canvas[x]->SaveAs(histoname1);
            delete canvas[x];
        }
    }
    //END OF RATIOS HISTOGRAMS DRAW

    TCanvas* hfpcanvas = new TCanvas();
    hfpcanvas->cd();
    colPlot_1->SetMarkerStyle(2); colPlot_1->SetMarkerColor(2);
    colPlot_1->GetYaxis()->SetRangeUser(0.96,1.04);  colPlot_1->GetXaxis()->SetRangeUser(0,150);
    colPlot_1->GetYaxis()->SetTitle("(Q2/Q1)/(Q2_iEta30/Q1_iEta30)"); colPlot_1->GetXaxis()->SetTitle("Days");
    colPlot_1->Draw("AP");

    colPlot_2->SetLineWidth(2);
    colPlot_2->SetMarkerStyle(2); colPlot_2->SetMarkerColor(2);
    colPlot_2->Draw("P");

    colPlot_3->SetLineWidth(2);
    colPlot_3->SetMarkerStyle(2); colPlot_3->SetMarkerColor(2);
    colPlot_3->Draw("P");

    colPlot_4->SetLineWidth(2);
    colPlot_4->SetMarkerStyle(2); colPlot_4->SetMarkerColor(2);
    colPlot_4->Draw("P");

    colPlot_5->SetLineWidth(2);
    colPlot_5->SetMarkerStyle(2); colPlot_5->SetMarkerColor(6);
    colPlot_5->Draw("P");

    colPlot_6->SetLineWidth(2);
    colPlot_6->SetMarkerStyle(2); colPlot_6->SetMarkerColor(6);
    colPlot_6->Draw("P");

    colPlot_7->SetLineWidth(2);
    colPlot_7->SetMarkerStyle(2); colPlot_7->SetMarkerColor(6);
    colPlot_7->Draw("P");

    colPlot_8->SetLineWidth(2);
    colPlot_8->SetMarkerStyle(2); colPlot_8->SetMarkerColor(6);
    colPlot_8->Draw("P");

    colPlot_9->SetLineWidth(2); colPlot_9->SetMarkerStyle(2);
    colPlot_9->SetMarkerColor(7);
    colPlot_9->Draw("P");

    colPlot_10->SetLineWidth(2); colPlot_10->SetMarkerStyle(2);
    colPlot_10->SetMarkerColor(7);
    colPlot_10->Draw("P");

    colPlot_11->SetLineWidth(2);
    colPlot_11->SetMarkerStyle(2); colPlot_11->SetMarkerColor(7);
    colPlot_11->Draw("P");

    colPlot_12->SetLineWidth(2);
    colPlot_12->SetMarkerStyle(2); colPlot_12->SetMarkerColor(7);
    colPlot_12->Draw("P");

    colPlot_13->SetLineWidth(2);
    colPlot_13->SetMarkerStyle(2); colPlot_13->SetMarkerColor(8);
    colPlot_13->Draw("P");

    colPlot_14->SetLineWidth(2);
    colPlot_14->SetMarkerStyle(2); colPlot_14->SetMarkerColor(8);
    colPlot_14->Draw("P");

    colPlot_15->SetLineWidth(2);
    colPlot_15->SetMarkerStyle(2); colPlot_15->SetMarkerColor(8);
    colPlot_15->Draw("P");

    colPlot_16->SetLineWidth(2);
    colPlot_16->SetMarkerStyle(2); colPlot_16->SetMarkerColor(8);
    colPlot_16->Draw("P");

    colPlot_17->SetLineWidth(2);
    colPlot_17->SetMarkerStyle(2); colPlot_17->SetMarkerColor(9);
    colPlot_17->Draw("P");

    colPlot_18->SetLineWidth(2);
    colPlot_18->SetMarkerStyle(2); colPlot_18->SetMarkerColor(9);
    colPlot_18->Draw("P");

    colPlot_19->SetLineWidth(2);
    colPlot_19->SetMarkerStyle(2); colPlot_19->SetMarkerColor(9);
    colPlot_19->Draw("P");

    colPlot_20->SetLineWidth(2);
    colPlot_20->SetMarkerStyle(2); colPlot_20->SetMarkerColor(9);
    colPlot_20->Draw("P");

    colPlot_21->SetLineWidth(2);
    colPlot_21->SetMarkerStyle(2); colPlot_21->SetMarkerColor(11);
    colPlot_21->Draw("P");

    colPlot_22->SetLineWidth(2);
    colPlot_22->SetMarkerStyle(2); colPlot_22->SetMarkerColor(11);
    colPlot_22->Draw("P");

    colPlot_23->SetLineWidth(2);
    colPlot_23->SetMarkerStyle(2); colPlot_23->SetMarkerColor(11);
    colPlot_23->Draw("P");

    colPlot_24->SetLineWidth(2);
    colPlot_24->SetMarkerStyle(2); colPlot_24->SetMarkerColor(11);
    colPlot_24->Draw("P");

    TLegend* legend = new TLegend(0.2,0.2,0.4,0.4);
    legend->SetHeader("Channels");
    legend->AddEntry(colPlot_1,"(32)","p");
    //legend->AddEntry(colPlot_2,"(32 21)","p");
    //legend->AddEntry(colPlot_3,"(32 37)","p");
    //legend->AddEntry(colPlot_4,"(32 57)","p");
    legend->AddEntry(colPlot_5,"(34)","p");
    //legend->AddEntry(colPlot_6,"(34 21)","p");
    //legend->AddEntry(colPlot_7,"(34 37)","p");
    //legend->AddEntry(colPlot_8,"(34 57)","p");
    legend->AddEntry(colPlot_9,"(36)","p");
    legend->AddEntry(colPlot_13,"(38)","p");
    legend->AddEntry(colPlot_17,"(40)","p");
    legend->AddEntry(colPlot_21,"(41)","p");
    legend->Draw();

    hfpcanvas->SaveAs("Onyuzmil.pdf");
    delete hfpcanvas;


    //RADIATION DAMAGE CALCULATIONS
    float ilumi2018 = 6367, ilumiAugDec2018 = 2393, lumiRatio;
    lumiRatio = ilumiAugDec2018/ilumi2018;
    //cout << "lumi ratio = " << lumiRatio << endl;

    fstream raddamCsv;
    raddamCsv.open("RadDam.csv",fstream::out);

    if(!raddamCsv.is_open()) {
        cout << "error in raddam csv" <<endl;
    }

    TH1F* raddamHist = new TH1F("RadDam Histogram", "RadDam Histogram", 50, 0.0, 4.0);

    TGraphErrors* rad_dam = new TGraphErrors(56,chans,max_min,0,max_minError);
    rad_dam->SetName("Radiation Damage Estimation Min/Delta");rad_dam->SetTitle("Radiation Damage Estimation Min/Delta");

    TCanvas* canv2 = new TCanvas("page","page",5000,6000);
    canv2->Divide(1,2,0.001,0.001);

    canv2->cd(1);
    rad_dam->SetMarkerStyle(7); rad_dam->SetMarkerColor(3); rad_dam->SetMarkerSize(9);
    //rad_dam->GetYaxis()->SetRangeUser(0,3);
    rad_dam->GetYaxis()->SetTitle("(1-Min)/(Max-min)"); rad_dam->GetXaxis()->SetTitle("Channels");
    rad_dam->Draw("AP");

    for(int i = 0; i < 56; i++){
        max_min[i] /= lumiRatio;
        raddamHist->Fill(max_min[i]);

        if(i == 55) raddamCsv << RADDAM_CH[i].ieta <<";" << RADDAM_CH[i].iphi <<";" << max_min[i] <<";" << max_minError[i];
        else raddamCsv << RADDAM_CH[i].ieta <<";" << RADDAM_CH[i].iphi <<";" << max_min[i] <<";" << max_minError[i] << endl;
    }
    wa_minmax /= lumiRatio;
    float wa = wa_minmax/wa_minmax_denom;
    cout << "WEIGHTED AVERAGE IS " << wa_minmax << " / " << wa_minmax_denom << " = "<<wa << endl;

    TGraphErrors* rad_dam2 = new TGraphErrors(56,chans,max_min,0,max_minError);
    rad_dam2->SetName("Radiation Damage Estimation");rad_dam2->SetTitle("Radiation Damage Estimation");
    canv2->cd(2);
    rad_dam2->GetYaxis()->SetRangeUser(0,3);
    rad_dam2->SetMarkerStyle(7); rad_dam2->SetMarkerColor(4); rad_dam->SetMarkerSize(9);
    rad_dam2->GetYaxis()->SetTitle("Ratio"); rad_dam2->GetXaxis()->SetTitle("Channels");
    rad_dam2->Draw("AP");

    canv2->SaveAs("Raddam.pdf");
    delete canv2;

    TCanvas* canv3 = new TCanvas("page","page",5000,6000);
    canv3->cd();
    raddamHist->Draw();
    canv3->SaveAs("Raddam2.pdf");
    delete canv3;

    raddamCsv.close();

    //END OF RADIATION DAMAGE CALCULATIONS

    //SYSTEMATIC ERRORS FILE WRITE
    TH1F* tdcErr = new TH1F("TDC Errors", "TDC Errors", 100, 0., 0.005);
    tdcErr->SetCanExtend(TH1::kAllAxes);
    tdcErr->GetXaxis()->SetTitle("Systematic Error Rates"); tdcErr->GetXaxis()->CenterTitle();
    tdcErr->GetYaxis()->SetTitle("Channels"); tdcErr->GetYaxis()->CenterTitle();
    tdcErr->SetLineColor(46); tdcErr->SetFillColor(46);

    TH1F* dtError = new TH1F("Delay Time Errors", "DT Errors", 100, 0., 0.005);
    dtError->SetCanExtend(TH1::kAllAxes);
    dtError->GetXaxis()->SetTitle("Systematic Error Rates"); dtError->GetXaxis()->CenterTitle();
    dtError->GetYaxis()->SetTitle("Channels"); dtError->GetYaxis()->CenterTitle();
    dtError->SetLineColor(40); dtError->SetFillColor(40);


    char filename[100]; char command[100];

    fstream sysErrors;
    sysErrors.open("SystematicErrors.txt",fstream::out);

    fstream sysErrorsWA;
    sysErrorsWA.open("waverage.txt",fstream::out);

    if(!sysErrors.is_open()) {
        cout << "error" <<endl;
    }
    if(!sysErrorsWA.is_open()) {
        cout << "error2" <<endl;
    }

    /*float wi = 0, miwi = 0;
    for(int l =1; l<numberofruns; l++){
        if(percentageSTD[l][0] != 0) {
            miwi += percentageSTD[l][0]*percentageSTD[l][1];
            wi += percentageSTD[l][1];
        }
    }
    waverage = miwi/wi;
    ChannelList[k].sysTDCerror = waverage;*/

    sysErrors << "iEta iPhi DelayTimeErrors TDCErrors" << endl;
    sysErrorsWA << "DelayTimeErrors TDCErrors" << endl;
    float wiTDC = 0, wimiTDC = 0, waTDC, wiDT = 0, wimiDT = 0, waDT;

    for(int i1=0;i1<ChannelList.size();i1++) {
        tdcErr->Fill(ChannelList[i1].sysTDCerror);
        dtError->Fill(ChannelList[i1].sysDTerror);
        if(i1 == 55) sysErrors << RADDAM_CH[i1].ieta << " " <<RADDAM_CH[i1].iphi << " " << ChannelList[i1].sysDTerror << " " << ChannelList[i1].sysTDCerror ;
        else sysErrors << RADDAM_CH[i1].ieta << " " <<RADDAM_CH[i1].iphi <<" " << ChannelList[i1].sysDTerror << " " << ChannelList[i1].sysTDCerror << " " << endl;

        //cout << "tdc : " << wimiTDC << " / " << wiTDC << " error : " << ChannelList[i1].sysTDCerror << endl;
        if(!isnan(ChannelList[i1].sysTDCerror)){
            wimiTDC += 1/ChannelList[i1].sysTDCerror;
            wiTDC += 1/(ChannelList[i1].sysTDCerror*ChannelList[i1].sysTDCerror);
        }
        if(!isnan(ChannelList[i1].sysDTerror)){
            wimiDT += 1/ChannelList[i1].sysDTerror;
            wiDT += 1/(ChannelList[i1].sysDTerror*ChannelList[i1].sysDTerror);
        }

        /*if(i1%4 == 3){
            waTDC = wimiTDC/wiTDC;
            wimiTDC = 0; wiTDC = 0;

            waDT = wimiDT/wiDT;
            wimiDT = 0; wiDT = 0;

            if(i1 == 55) sysErrorsWA << RADDAM_CH[i1].ieta << " " << waDT << " " << waTDC ;
            else sysErrorsWA << RADDAM_CH[i1].ieta << " " << waDT << " " << waTDC << " " << endl;
        }*/
    }
    waTDC = wimiTDC/wiTDC; waDT = wimiDT/wiDT;
    sysErrorsWA << waDT << " " << waTDC << endl;

    sysErrors.close();
    sysErrorsWA.close();

    TCanvas* ca = new TCanvas("page","page",5000,6000);
    ca->Divide(1,2,0.001,0.001);

    ca->cd(1);
    tdcErr->Draw();

    ca->cd(2);
    dtError->Draw();

    ca->SaveAs("SysErrors.pdf");
    delete ca;
    //END OF SYSTEMATIC ERRORS FILE WRITE

    //RATIOS 3D ETA HISTOGRAM DRAW
    for(int k = 0; k<ChannelList.size(); k++){

        if(k%4 == 0){
            /* 3D Histogram of phi, ratio and days */
            char hname[100];
            sprintf(hname, "; Days; iPhi; Normalized Ratios");
            ChannelList[k].ratio2D = new TGraph2D();
            ChannelList[k].ratio2D->SetTitle(hname);
            ChannelList[k].ratio2D->GetXaxis()->CenterTitle(); ChannelList[k].ratio2D->GetYaxis()->CenterTitle(); ChannelList[k].ratio2D->GetZaxis()->CenterTitle();
        }

        float A = ChannelList[k].runForOneChannel[0].ratioValue;
        int numberofruns = ChannelList[k].runForOneChannel.size();

        for(int l =0; l<numberofruns; l++){

            float B = ChannelList[k].runForOneChannel[l].ratioValue;
            ChannelList[k-k%4].ratio2D->SetPoint(ChannelList[k-k%4].ratio2D->GetN(), ChannelList[k].runForOneChannel[l].datetime, RADDAM_CH[k].iphi,B/A);
        }

    }
    //END OF RATIOS ETA HISTOGRAMS DRAW

    //ALL RATIOS 3D HISTOGRAM DRAW
    TGraph2D* ratio2DHFP = new TGraph2D();
    ratio2DHFP->SetTitle("; Days; iEta; Weighted Averages");
    ratio2DHFP->GetXaxis()->CenterTitle(); ratio2DHFP->GetYaxis()->CenterTitle(); ratio2DHFP->GetZaxis()->CenterTitle();
    TGraph2D* ratio2DHFM = new TGraph2D();
    ratio2DHFM->SetTitle("; Days; iEta; Weighted Averages");
    ratio2DHFM->GetXaxis()->CenterTitle(); ratio2DHFM->GetYaxis()->CenterTitle(); ratio2DHFM->GetZaxis()->CenterTitle();

    int numberofruns = ChannelList[0].runForOneChannel.size();
    float weightedaverage[numberofruns][14] = {{0}};
    float WAnormalization[numberofruns][14] = {{0}};
    float A, DA, B, DB, er, erSquare;

    for(int k = 0; k<ChannelList.size(); k++){

        A = ChannelList[k].runForOneChannel[0].ratioValue;
        DA = ChannelList[k].runForOneChannel[0].ratioError;
        numberofruns = ChannelList[k].runForOneChannel.size();

        for(int l =0; l<numberofruns; l++){

            B = ChannelList[k].runForOneChannel[l].ratioValue;
            DB = ChannelList[k].runForOneChannel[l].ratioError;

            if(B != 0){
                er = (B/A)*sqrt((DA/A)*(DA/A) + (DB/B)*(DB/B));
                erSquare = 1/(er*er);
                WAnormalization[l][k/4] += erSquare;
                weightedaverage[l][k/4] += (B/A)*erSquare;
            }

            if(k%4 == 3){//The last ieta one

                int pointi2 = ChannelList[k].ratioEtas->GetN();
                ChannelList[k].ratioEtas->SetPoint(pointi2,ChannelList[k].runForOneChannel[l].datetime,(weightedaverage[l][k/4])/WAnormalization[l][k/4]);
                ChannelList[k].ratioEtas->SetPointError(pointi2,0,sqrt(1/WAnormalization[l][k/4]));

                if(RADDAM_CH[k].ieta > 0){
                    ratio2DHFP->SetPoint(ratio2DHFP->GetN(),ChannelList[k].runForOneChannel[l].datetime, RADDAM_CH[k].ieta,(weightedaverage[l][k/4])/WAnormalization[l][k/4]);
                } else {
                    ratio2DHFM->SetPoint(ratio2DHFM->GetN(),ChannelList[k].runForOneChannel[l].datetime, -1*RADDAM_CH[k].ieta,(weightedaverage[l][k/4])/WAnormalization[l][k/4]);
                }
            }
        }
    }

    //ETA HISTOGRAMS DRAW
    int dividexaxes3d = 1;
    int divideyaxes3d = 2;
    int pagenumber3d = 7;

    string histonames2d[]= {"Ratios2d.pdf","RatiosEta.pdf"};

    int histonumber2 =  sizeof(histonames2d)/sizeof(histonames2d[0]);

    //sprintf(histoname, "Ratios2d.pdf");
    for( int n=0 ; n<histonumber2 ; n++){

        int l = histonames2d[n].length();

        char histoname1[20]; char histoname[20] ;
        histonames2d[n].copy(histoname,l,0);

        for(int foo=l; foo<20;foo++){
            histoname[foo]= '\0';
        }

        TCanvas* canvas3d[pagenumber3d];

        for(int x=0;x<pagenumber3d;x++){
            char name[100];
            sprintf(name, "page%d",x);
            canvas3d[x] = new TCanvas(name,name,5000,6000);
            gStyle->SetOptStat(0);
            gStyle->SetTitleFontSize(0.1);
            canvas3d[x]->Divide(dividexaxes3d,divideyaxes3d,0.001,0.001);
        }
        int currentpage = 0; int currentcanv = 1;

        for(int o=0 ; o<14; o++){
            canvas3d[currentpage]->cd(currentcanv);
            currentcanv++;
            if(currentcanv> divideyaxes3d*dividexaxes3d ) {
                currentpage++;
                currentcanv = 1;
            }
            if ( n == 0 ){ //2D ETA HISROGRAM
                ChannelList[o*4].ratio2D->GetHistogram()->SetMaximum(1.15); ChannelList[o*4].ratio2D->GetHistogram()->SetMinimum(0.85);
                ChannelList[o*4].ratio2D->Draw("lego2z");
                ChannelList[o*4].ratio2D->Write();
            } else if ( n == 1 ) {
                ChannelList[o*4+3].ratioEtas->GetHistogram()->SetMaximum(1.15); ChannelList[o*4+3].ratioEtas->GetHistogram()->SetMinimum(0.85);
                ChannelList[o*4+3].ratioEtas->Draw("A*");
                ChannelList[o*4+3].ratioEtas->Write();
            }

        }
        for(int x=0;x<pagenumber3d; x++){

            char dummy[10];

            if( x==0 ) sprintf(dummy,"("); //First page of pdf file
            else if(x==pagenumber3d-1) sprintf(dummy,")"); //Last page of pdf file
            else sprintf(dummy,"");

            sprintf(histoname1,"%s%s",histoname,dummy);
            canvas3d[x]->SaveAs(histoname1);
            delete canvas3d[x];
        }
    }
    //END OF ETA HISTOGRAMS DRAW

    TCanvas* canvas3 = new TCanvas("page","page",5000,6000);
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.1);
    canvas3->Divide(1,2,0.001,0.001);

    canvas3->cd(1);
    ratio2DHFP->GetHistogram()->SetMaximum(1.15); ratio2DHFP->GetHistogram()->SetMinimum(0.85);
    ratio2DHFP->Draw("lego2Z");
    ratio2DHFP->Write();

    canvas3->cd(2);
    ratio2DHFM->GetHistogram()->SetMaximum(1.15); ratio2DHFM->GetHistogram()->SetMinimum(0.85);
    ratio2DHFM->Draw("lego2Z");
    ratio2DHFM->Write();

    canvas3->SaveAs("RatiosAll2D.pdf");
    delete canvas3;
    //END OF ALL RATIOS HISTOGRAMS DRAW

    float etasP[7] = {30,32,34,36,38,40,41};

    TGraph* dayZero1 = new TGraph(7,etasP,arrayHFP); dayZero1->SetTitle("");
    dayZero1->GetXaxis()->SetTitle("iEta"); dayZero1->GetYaxis()->SetTitle("Normalized Error");
    dayZero1->GetYaxis()->SetRangeUser(0,0.025); dayZero1->GetXaxis()->SetRangeUser(0,42);
    dayZero1->SetMarkerStyle(8); dayZero1->SetMarkerColor(46); dayZero1->SetMarkerSize(5);
    //dayZero1->Fit("expo");

    TGraph* dayZero2 = new TGraph(7,etasP,arrayHFM); //dayZero2->SetTitle("HF-");
    //dayZero2->GetXaxis()->SetTitle("Eta x(-1)"); dayZero2->GetYaxis()->SetTitle("Normalized Percentage Error");
    dayZero2->SetMarkerStyle(29); dayZero2->SetMarkerColor(8); dayZero2->SetMarkerSize(5);

    TGraph* dayZero3 = new TGraph(7,etasP,arrayHF); dayZero3->SetTitle("HF");
    dayZero3->GetYaxis()->SetRangeUser(0,0.025); dayZero3->GetXaxis()->SetRangeUser(0,42);
    dayZero3->SetMarkerStyle(3); dayZero3->SetMarkerColor(9); dayZero3->SetMarkerSize(5);

    TCanvas* canvas4 = new TCanvas("page","page",5000,6000);
    canvas4->Divide(1,2,0.001,0.001);
    //gStyle->SetOptFit(0111);
    canvas4->cd(1); dayZero3->Draw("P");//dayZero1->Write();
    canvas4->cd(2); dayZero1->Draw("A*"); dayZero2->Draw("P"); //dayZero2->Write();
    canvas4->SaveAs("day0Eta.pdf");

    system("mv Ratios.root ../Histos/Ratios.root");

    outfile->Close();
    return ;
}

int main(int argc, char *argv[]) {

    getmap();

    if(atoi(argv[1]) == 0) { //Draw the histograms of Q2/Q1 ratios of all runs
        plotratios();
    } else {
        RunNo=atoi(argv[1]);//Converts string to integer
        plotleds();
    }
}

