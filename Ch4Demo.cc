#include "prlite_genmat.hpp"
#include "emdw.hpp"
// #include "sqrtmvg.hpp"
#include "affinemvg.hpp"
#include "sqrtmvgtable.hpp"

#include "lbp_cg.hpp"
#include "lbu_cg.hpp"


#include <iostream>
#include <vector>
#include <fstream>

#include <chrono>
// #include <ctime>

using namespace std;
using namespace emdw;

typedef SqrtMVG SG;
typedef AffineMVG AG;
typedef unsigned T;
typedef DiscreteTable<T> DT;
typedef std::chrono::high_resolution_clock Clock;

vector<rcptr<Factor>> discreteGraph(uint discreteBase, uint timestep, uint objects, uint measurements, vector<rcptr<Factor>> messages);

int main(int argc, char *argv[]) {
  try {

    enum outputLevel {None, Error, Verbose, Debug, Periphrastic};
    outputLevel OL = Error;
    if(argc > 2){
      switch (argv[2][0])
      {
      case 'n':
      case 'N':
        OL = None;
        break;
      case 'e':
      case 'E':
        OL = Error;
        break;
      case 'd':
      case 'D':
        OL = Debug;
        break;
      case 'v':
      case 'V':
        OL = Verbose;
        break;
      case 'p':
      case 'P':
        OL = Periphrastic;
        break;
      default:
        break;
      }
    }

    //Multi  Object Discrete Setup:
    uint objects = atoi(argv[3]);
    uint measurements = atoi(argv[4]);

    vector<vector<double>> measXList(measurements);
    vector<vector<double>> measYList(measurements);

    char c;
    std::string str;

    ifstream inFile;

    inFile.open("Xpos.txt");
    if (!inFile) {
        if(OL >= Error) cout << "Unable to open file";
        exit(1); // terminate with error
    }

    int measurementReadCounter = -1;
    while (inFile >> c) {
      switch(c) {
        case '[':
          measurementReadCounter++;
          if(OL >= Verbose) cout << "start x parse" << measurementReadCounter << endl;
          break;
        case ']':
        case ',':
          measXList[measurementReadCounter].push_back(std::stod(str));
          str.clear();
          break;
        case ' ':
        case '\n':
          break;
        default:
          str.push_back(c);
      }
    }

    inFile.close();
    if(OL >= Verbose) cout << "end x parse" << endl;

    inFile.open("Ypos.txt");
    if (!inFile) {
        if(OL >= Error) cout << "Unable to open file";
        exit(1); // terminate with error
    }

    measurementReadCounter = -1;
    while (inFile >> c) {
      switch(c) {
        case '[':
          measurementReadCounter++;
          if(OL >= Verbose) cout << "start Y parse " << measurementReadCounter << endl;
          break;
        case ',':
        case ']':
          measYList[measurementReadCounter].push_back(std::stod(str));
          str.clear();
          break;
        case ' ':
        case '\n':
          break;
        default:
          str.push_back(c);
      }
    }

    inFile.close();
    if(OL >= Verbose) cout << "end y parse" << endl;

    if(OL >= Debug) cout << "measYList \n" << measYList << "\n measXList \n" << measXList << endl;

    uint timesteps = measXList[0].size();
    if (timesteps < 2){
      if(OL >= Error) cout << "At least 2 timesteps required";
        exit(1); // terminate with error
    }
    if (timesteps != measYList[0].size()){
      if(OL >= Error) cout << "measurement dimension mismatch";
        exit(1); // terminate with error
    }

    vector<rcptr<Factor>> posPdfs; //TODO convert this to a specific size for effeciencies
    vector<rcptr<Factor>> measHdfs;

    uint stateVectorSize = 6;

    // Discrete inits for hybrid Setup:
    uint DiscreteBase = stateVectorSize*timesteps*objects;
    rcptr< vector<T> > binDom( new vector<T>{0,1} ); // should we define these globally?
    double clutterProb = 0.05;
    prlite::ColVector<double> uninformMn(2);
    prlite::RowMatrix<double> uninforCov(2,2);
    uninformMn[0] = 621; uninformMn[1] = 187;
    uninforCov(0,0) = 2000.0; uninforCov(0,1) = 0.0;
    uninforCov(1,0) = 0.0; uninforCov(1,1) = 2000.0;

    //construct the cluster graph
    for(uint t = 0; t<timesteps-1; t++){
      uint tCID = t*(stateVectorSize*objects);
      uint FtCID = (t+1)*(stateVectorSize*objects);
      uint tDID = DiscreteBase+t*objects*measurements;
      for(uint o = 0; o<objects; o++){
        uint oCID = tCID+stateVectorSize*o;//fog
        uint FoCID = FtCID+stateVectorSize*o;//fog

        //constructs an affine gaussian for the point to point relations that includes an unknown velocity and acceleration
        prlite::ColMatrix<double> A(6,6);
        prlite::ColVector<double> C(6);
        prlite::ColMatrix<double> R(6,6);
        A(0,0)=1; A(0,1)=0; A(0,2)=1; A(0,3)=0; A(0,4)=0.5; A(0,5)=0;
        A(1,0)=0; A(1,1)=1; A(1,2)=0; A(1,3)=1; A(1,4)=0; A(1,5)=0.5;
        A(2,0)=0; A(2,1)=0; A(2,2)=1; A(2,3)=0; A(2,4)=1; A(2,5)=0;
        A(3,0)=0; A(3,1)=0; A(3,2)=0; A(3,3)=1; A(3,4)=0; A(3,5)=1;
        A(4,0)=0; A(4,1)=0; A(4,2)=0; A(4,3)=0; A(4,4)=1; A(4,5)=0;
        A(5,0)=0; A(5,1)=0; A(5,2)=0; A(5,3)=0; A(5,4)=0; A(5,5)=1;

        C[0]=0; C[1]=0; C[2]=0; C[3]=0; C[4]=0; C[5]=0;

        R(0,0)=1; R(0,1)=0; R(0,2)=0; R(0,3)=0; R(0,4)=0; R(0,5)=0;
        R(1,0)=0; R(1,1)=1; R(1,2)=0; R(1,3)=0; R(1,4)=0; R(1,5)=0;
        R(2,0)=0; R(2,1)=0; R(2,2)=0.01; R(2,3)=0; R(2,4)=0; R(2,5)=0;
        R(3,0)=0; R(3,1)=0; R(3,2)=0; R(3,3)=0.01; R(3,4)=0; R(3,5)=0;
        R(4,0)=0; R(4,1)=0; R(4,2)=0; R(4,3)=0; R(4,4)=0.01; R(4,5)=0;
        R(5,0)=0; R(5,1)=0; R(5,2)=0; R(5,3)=0; R(5,4)=0; R(5,5)=0.01;

        posPdfs.push_back(uniqptr<AG>(new AG({oCID,oCID+1,oCID+2,oCID+3,oCID+4,oCID+5},{FoCID,FoCID+1,FoCID+2,FoCID+3,FoCID+4,FoCID+5}, A, C, R)));

        for(uint m = 0; m<measurements;m++){
          uint DID = tDID+m*objects+o;//check that this aligns with the discrete structure we construct
          prlite::ColVector<double> mn(2);
          prlite::RowMatrix<double> cov(2,2);
          mn[0] = measXList[m][t]; mn[1] = measYList[m][t];
          cov(0,0) = 10.0; cov(0,1) = 0.0;
          cov(1,0) = 0.0; cov(1,1) = 10.0;
          map< vector<T>, tuple< double, prlite::ColVector<double>, prlite::RowMatrix<double> > > weightedComponents =
          {{{0}, { clutterProb, uninformMn, uninforCov}},
          {{1}, { 1.0-clutterProb, mn, cov}}};
          measHdfs.push_back(uniqptr<SqrtMVGTable<T>>( new SqrtMVGTable<T>({DID}, {{oCID,oCID+1}}, {binDom}, weightedComponents)));
        }
      }
    }

    vector<rcptr<Factor>> expect(objects);

    //constructs the final measurements mvg tables as there is one more than the number of clusters
    {//scope for final measurement clusters
      uint tCID = (timesteps-1)*objects*stateVectorSize;
      uint tDID = DiscreteBase+(timesteps-1)*objects*measurements;
      for(uint o = 0; o<objects; o++){
        uint oCID = tCID+stateVectorSize*o;
        for(uint m = 0; m<measurements;m++){
          uint DID = tDID+m*objects+o;
          prlite::ColVector<double> mn(2);
          prlite::RowMatrix<double> cov(2,2);
          mn[0] = measXList[m][timesteps-1]; mn[1] = measYList[m][timesteps-1];
          cov(0,0) = 20.0; cov(0,1) = 0.0;
          cov(1,0) = 0.0; cov(1,1) = 20.0;
          map< vector<T>, tuple< double, prlite::ColVector<double>, prlite::RowMatrix<double> > > weightedComponents =
          {{{0}, { clutterProb, uninformMn, uninforCov}},
          {{1}, { 1.0-clutterProb, mn, cov}}};
          measHdfs.push_back(uniqptr<SqrtMVGTable<T>>( new SqrtMVGTable<T>({DID}, {{oCID,oCID+1}}, {binDom}, weightedComponents)));
        }

        //construct general expectations of state variables
        oCID = stateVectorSize*o; //TODO Check that this is right
        prlite::ColVector<double> mnEx(6);
        prlite::RowMatrix<double> covEx(6,6);
        mnEx[0] = measXList[o][0]; mnEx[1] = measYList[o][0]; mnEx[2] = 0; mnEx[3] = 0; mnEx[4] = 0; mnEx[5] = 0;
        covEx(0,0) = 10.0; covEx(0,1) = 0.0; covEx(0,2) = 0.0; covEx(0,3) = 0.0; covEx(0,4) = 0.0; covEx(0,5) = 0.0;
        covEx(1,0) = 0.0; covEx(1,1) = 10.0; covEx(1,2) = 0.0; covEx(1,3) = 0.0; covEx(1,4) = 0.0; covEx(1,5) = 0.0;
        covEx(2,0) = 0.0; covEx(2,1) = 0.0; covEx(2,2) = 100.0; covEx(2,3) = 0.0; covEx(2,4) = 0.0; covEx(2,5) = 0.0;
        covEx(3,0) = 0.0; covEx(3,1) = 0.0; covEx(3,2) = 0.0; covEx(3,3) = 100.0; covEx(3,4) = 0.0; covEx(3,5) = 0.0;
        covEx(4,0) = 0.0; covEx(4,1) = 0.0; covEx(4,2) = 0.0; covEx(4,3) = 0.0; covEx(4,4) = 100.0; covEx(4,5) = 0.0;
        covEx(5,0) = 0.0; covEx(5,1) = 0.0; covEx(5,2) = 0.0; covEx(5,3) = 0.0; covEx(5,4) = 0.0; covEx(5,5) = 100.0;
        expect[o] = uniqptr<SG>(new SG({oCID,oCID+1,oCID+2,oCID+3,oCID+4,oCID+5}, mnEx, covEx));
      }
    }


    vector<rcptr<Factor>> clustBeliefs((timesteps-1)*objects+timesteps*measurements*objects); //TODO Really Really should be buildng these with multidimensional arrays.... as this is borderline stupid
    vector<rcptr<Factor>> sepsetBeliefs((timesteps-2)*objects+timesteps*measurements*objects); // (t-2)*o+t*m*o -> posSepsets+measSepsets
    vector<rcptr<Factor>> interfaceDownMessages(timesteps*measurements*objects);
    vector<rcptr<Factor>> interfaceUpMessages(timesteps*measurements*objects);
    vector<rcptr<Factor>> newUpMessages(measurements*objects);
    rcptr<Factor> oldSepBel;

    // construct cluster beliefs
    for(uint t = 0; t<timesteps-1;t++){
      uint tBelIndex = t*(1+measurements)*objects;
      uint tMeasIndex = t*measurements*objects;
      uint tPosIndex = t*objects;
      for(uint o = 0; o<objects; o++){
        uint oBelIndex = tBelIndex+(1+measurements)*o;
        uint oMeasIndex = tMeasIndex+measurements*o;
        clustBeliefs[oBelIndex] = posPdfs[tPosIndex+o];
        for(uint m = 0; m<measurements;m++){
          clustBeliefs[oBelIndex+1+m] = measHdfs[oMeasIndex+m];
        }
      }
    }

    //construct final measurment cluster beliefs
    uint tBelIndex = (timesteps-1)*(1+measurements)*objects;
    uint tMeasIndex = (timesteps-1)*measurements*objects;
    for(uint o = 0; o<objects; o++){
      uint oBelIndex = tBelIndex+measurements*o;
      uint oMeasIndex = tMeasIndex+measurements*o;
      for(uint m = 0; m<measurements;m++){
        clustBeliefs[oBelIndex+m] = measHdfs[oMeasIndex+m]; //Misty
      }
    }

    if(OL >= Verbose) cout << "completed construction" << endl;

    //init all objects
    for(uint o = 0; o<objects; o++){ //check this init indexing
      uint CID = stateVectorSize*o;
      uint FCID = stateVectorSize*objects+stateVectorSize*o;
      clustBeliefs[o*(1+measurements)] = clustBeliefs[o*(1+measurements)]->absorb(expect[o])->marginalize({CID,CID+1,CID+2,CID+3,CID+4,CID+5,FCID,FCID+1,FCID+2,FCID+3,FCID+4,FCID+5}); //fog
    }

    //send messages forwards
    for(uint t = 0; t<timesteps-2;t++){
      uint tBelIndex = t*(1+measurements)*objects;
      uint tMesIndex = t*measurements*objects;
      uint tCID = t*(stateVectorSize*objects);
      uint FtCID = (t+1)*(stateVectorSize*objects);
      uint FFtCID = (t+2)*(stateVectorSize*objects);
      uint tDID = DiscreteBase+t*objects*measurements;

      if(OL >= Verbose) cout << "\nforward loop " << t << endl;

      //update the hybrid measurement clusters from positional info and construct down messages
      for(uint o = 0; o<objects; o++){
        uint oBelIndex = tBelIndex+(1+measurements)*o;
        uint oMesIndex = tMesIndex+measurements*o;
        uint oCID = tCID+stateVectorSize*o;

        for(uint m = 0; m<measurements;m++){
          uint DID = tDID+m*objects+o;
          sepsetBeliefs[oBelIndex+m] = clustBeliefs[oBelIndex]->marginalize({oCID,oCID+1}); //TODO add special offset for the final clusts
          clustBeliefs[oBelIndex+1+m] = clustBeliefs[oBelIndex+1+m]->absorb(sepsetBeliefs[oBelIndex+m]);
          interfaceDownMessages[oMesIndex+m] = clustBeliefs[oBelIndex+1+m]->marginalize({DID}); //µ′a,b = Sum:Ψ /µb,a but µb,a uninformative at init
          // if(OL >= Debug) cout <<"DID " << DID << endl;
        }
      }
      if(OL >= Verbose) cout << "updated hybrid meas clusters and down interface messages" << endl;

      //run discrete belief update, extract up messages, update the hybrid measurement clusters from discrete inference info
      vector<rcptr<Factor>> DGBeliefs = discreteGraph(DiscreteBase, t, objects, measurements, interfaceDownMessages);
      for(uint o = 0; o < objects; o++){
        uint oMesIndex = o*measurements;
        uint oBelIndex = tBelIndex+(1+measurements)*o;
        if(OL >= Periphrastic) cout << "discrete extract " << o << endl;
        for(uint m = 0; m < measurements; m++){
          uint mMesIndex = oMesIndex+m;
          uint mBelIndex = oBelIndex+1+m; //+1 included for pos cluster offset
          newUpMessages[mMesIndex] = DGBeliefs[oMesIndex+m]->cancel(interfaceDownMessages[tMesIndex+oMesIndex+m]);//µ′a,b = Sum:Ψ /µb,a
          if(OL >= Periphrastic) cout << "discrete up message \n" << *newUpMessages[mMesIndex] << endl;
          clustBeliefs[mBelIndex] = clustBeliefs[mBelIndex]->absorb(newUpMessages[mMesIndex]); //Ψ′ = Ψ(µ′a,b/µa,b) but µa,b is uninformative
          interfaceUpMessages[tMesIndex+mMesIndex] = newUpMessages[mMesIndex]; //µ=µ′
        }
      }
      if(OL >= Verbose) cout << "updated discrete system and up interface messages" << endl;
      if(OL >= Verbose) cout << "updated hybrid meas clusters from discrete messages" << endl;

      //update current pos clusters
      for(uint o = 0; o < objects; o++){
        uint oBelIndex = tBelIndex+(1+measurements)*o;
        uint oCID = tCID+stateVectorSize*o;
        for(uint m = 0; m<measurements;m++){ //TODO check the values here
          uint mBelIndex = oBelIndex+m;
          oldSepBel = sepsetBeliefs[mBelIndex];
          sepsetBeliefs[mBelIndex] = clustBeliefs[mBelIndex+1]->marginalize({oCID,oCID+1});
          clustBeliefs[oBelIndex] = clustBeliefs[oBelIndex]->absorb(sepsetBeliefs[mBelIndex])->cancel(oldSepBel);
        }
      }
      if(OL >= Verbose) cout << "updated current pos cluster from measurements" << endl;

      //update next pos clusters
      for(uint o = 0; o < objects; o++){
        uint oBelIndex = tBelIndex+(1+measurements)*o;
        uint FoBelIndex = (t+1)*(1+measurements)*objects+(1+measurements)*o; //TODO check
        uint FoCID = FtCID+stateVectorSize*o;
        uint FFoCID = FFtCID+stateVectorSize*o;
        sepsetBeliefs[oBelIndex+measurements] = clustBeliefs[oBelIndex]->marginalize({FoCID,FoCID+1,FoCID+2,FoCID+3,FoCID+4,FoCID+5});
        clustBeliefs[FoBelIndex] = clustBeliefs[FoBelIndex]->absorb(sepsetBeliefs[oBelIndex+measurements])->marginalize({FoCID,FoCID+1,FoCID+2,FoCID+3,FoCID+4,FoCID+5,FFoCID,FFoCID+1,FFoCID+2,FFoCID+3,FFoCID+4,FFoCID+5});
      }
      if(OL >= Verbose) cout << "updated next pos cluster" << endl;
    }

    {//scope for final clusters
      uint tBelIndex = (timesteps-2)*(1+measurements)*objects;
      uint tMesIndex = (timesteps-2)*measurements*objects;
      uint tCID = (timesteps-2)*(stateVectorSize*objects);
      uint tDID = DiscreteBase+(timesteps-2)*objects*measurements;

      //update the second last hybrid measurement clusters from positional info
      for(uint o = 0; o<objects; o++){
        uint oClustBelIndex = tBelIndex+(measurements+1)*o;
        uint oSepBelIndex = tBelIndex+measurements*o;
        uint oMesIndex = tMesIndex+measurements*o;
        uint oCID = tCID+stateVectorSize*o;

        for(uint m = 0; m<measurements;m++){
          uint DID = tDID+m*objects+o;
          sepsetBeliefs[oSepBelIndex+m] = clustBeliefs[oClustBelIndex]->marginalize({oCID,oCID+1}); //TODO add special offset for the final clusts
          clustBeliefs[oClustBelIndex+1+m] = clustBeliefs[oClustBelIndex+1+m]->absorb(sepsetBeliefs[oSepBelIndex+m]);
          interfaceDownMessages[oMesIndex+m] = clustBeliefs[oClustBelIndex+1+m]->marginalize({DID}); //µ′a,b = Sum:Ψ /µb,a but µb,a uninformative at init
        }
      }
      if(OL >= Verbose) cout << "updated hybrid meas clusters and down interface messages" << endl;

      //run discrete belief update on second last hybrid measurement clusters
      vector<rcptr<Factor>> DGBeliefs = discreteGraph(DiscreteBase, (timesteps-2), objects, measurements, interfaceDownMessages);
      for(uint o = 0; o < objects; o++){
        uint oMesIndex = o*measurements;
        uint oBelIndex = tBelIndex+(1+measurements)*o;
        if(OL >= Periphrastic) cout << "discrete extract " << o << endl;
        for(uint m = 0; m < measurements; m++){
          uint mMesIndex = oMesIndex+m;
          uint mBelIndex = oBelIndex+1+m; //+1 included for pos cluster offset
          newUpMessages[mMesIndex] = DGBeliefs[mMesIndex]->cancel(interfaceDownMessages[tMesIndex+oMesIndex+m]);//µ′a,b = Sum:Ψ /µb,a
          if(OL >= Periphrastic) cout << "discrete up message \n" << *newUpMessages[mMesIndex] << endl;
          clustBeliefs[mBelIndex] = clustBeliefs[mBelIndex]->absorb(newUpMessages[mMesIndex]); //Ψ′ = Ψ(µ′a,b/µa,b) but µa,b is uninformative
          interfaceUpMessages[tMesIndex+mMesIndex] = newUpMessages[mMesIndex]; //µ=µ′
        }
      }
      if(OL >= Verbose) cout << "updated discrete system and up interface messages" << endl;
      if(OL >= Verbose) cout << "updated hybrid meas clusters from discrete messages" << endl;

      //update last pos clusters
      for(uint o = 0; o < objects; o++){
        uint oClustBelIndex = tBelIndex+(measurements+1)*o;
        uint oSepBelIndex = tBelIndex+measurements*o;
        uint oCID = tCID+stateVectorSize*o;
        for(uint m = 0; m<measurements;m++){
          uint mMesClustBelIndex = oClustBelIndex+1+m;
          uint mSepBelIndex = oSepBelIndex+m;
          oldSepBel = sepsetBeliefs[mSepBelIndex];
          sepsetBeliefs[mSepBelIndex] = clustBeliefs[mMesClustBelIndex]->marginalize({oCID,oCID+1});
          clustBeliefs[oClustBelIndex] = clustBeliefs[oClustBelIndex]->absorb(sepsetBeliefs[mSepBelIndex])->cancel(oldSepBel);
        }
      }
      if(OL >= Verbose) cout << "updated last pos cluster from second last measurements" << endl;

      tMesIndex = (timesteps-1)*measurements*objects;
      tCID = (timesteps-1)*(stateVectorSize*objects);
      tDID = DiscreteBase+(timesteps-1)*objects*measurements;

      //update the last hybrid measurement clusters from positional info
      for(uint o = 0; o<objects; o++){
        uint oClustBelIndex = tBelIndex+(measurements+1)*o;
        uint oSepBelIndex = tBelIndex+measurements*objects+measurements*o;
        uint oMesIndex = tMesIndex+measurements*o;//mist
        uint oCID = tCID+stateVectorSize*o;//mist

        for(uint m = 0; m<measurements;m++){
          uint mMesClustBelIndex = tBelIndex+(measurements+1)*objects+measurements*o+m;
          uint DID = tDID+m*objects+o;//mist
          sepsetBeliefs[oSepBelIndex+m] = clustBeliefs[oClustBelIndex]->marginalize({oCID,oCID+1}); //TODO add special offset for the final clusts
          clustBeliefs[mMesClustBelIndex] = clustBeliefs[mMesClustBelIndex]->absorb(sepsetBeliefs[oSepBelIndex+m]);
          interfaceDownMessages[oMesIndex+m] = clustBeliefs[mMesClustBelIndex]->marginalize({DID}); //µ′a,b = Sum:Ψ /µb,a but µb,a uninformative at init
        }
      }
      if(OL >= Verbose) cout << "updated hybrid meas clusters and down interface messages" << endl;

      //run discrete belief update on second last hybrid measurement clusters
      if(OL >= Verbose) cout << "updated discrete system and up interface messages" << endl;
      DGBeliefs = discreteGraph(DiscreteBase, (timesteps-1), objects, measurements, interfaceDownMessages);
      for(uint o = 0; o < objects; o++){
        uint oMesIndex = o*measurements;
        if(OL >= Periphrastic) cout << "discrete extract " << o << endl;
        for(uint m = 0; m < measurements; m++){
          uint mMesIndex = oMesIndex+m;
          uint mMesClustBelIndex = tBelIndex+(measurements+1)*objects+measurements*o+m;
          newUpMessages[mMesIndex] = DGBeliefs[mMesIndex]->cancel(interfaceDownMessages[tMesIndex+oMesIndex+m]);//µ′a,b = Sum:Ψ /µb,a
          if(OL >= Periphrastic) cout << "discrete up message \n" << *newUpMessages[mMesIndex] << endl;
          clustBeliefs[mMesClustBelIndex] = clustBeliefs[mMesClustBelIndex]->absorb(newUpMessages[mMesIndex]); //Ψ′ = Ψ(µ′a,b/µa,b) but µa,b is uninformative
          interfaceUpMessages[tMesIndex+mMesIndex] = newUpMessages[mMesIndex]; //µ=µ′
        }
      }
      if(OL >= Verbose) cout << "updated discrete system and up interface messages" << endl;
      if(OL >= Verbose) cout << "updated hybrid meas clusters from discrete messages" << endl;

      //update last pos clusters
      for(uint o = 0; o < objects; o++){ //TODO definitely needs to be checked
        uint oClustBelIndex = tBelIndex+(measurements+1)*o;
        uint oSepBelIndex = tBelIndex+measurements*objects+measurements*o;
        uint oCID = tCID+stateVectorSize*o;
        for(uint m = 0; m<measurements;m++){
          uint mMesClustBelIndex = tBelIndex+(measurements+1)*objects+measurements*o+m;
          uint mSepBelIndex = oSepBelIndex+m;
          oldSepBel = sepsetBeliefs[mSepBelIndex];
          sepsetBeliefs[mSepBelIndex] = clustBeliefs[mMesClustBelIndex]->marginalize({oCID,oCID+1});
          clustBeliefs[oClustBelIndex] = clustBeliefs[oClustBelIndex]->absorb(sepsetBeliefs[mSepBelIndex])->cancel(oldSepBel);
        }
      }
      if(OL >= Verbose) cout << "updated last pos cluster from last measurements" << endl;

      if(OL >= Verbose) cout << "completed forward messages" << endl;
    }//end scope for final clusters

    //send messages backwards
    for(int t = timesteps-2; t>0 ;t--){ //-1 from 0 count, -1 from clust amounts, -1 from sep amount !TODO this seems wrong?
      uint tBelIndex = t*(1+measurements)*objects;
      uint tMesIndex = t*measurements*objects;
      uint tCID = t*(stateVectorSize*objects);
      uint tDID = DiscreteBase+t*objects*measurements;


      if(OL >= Verbose) cout << "\nbackward loop " << t << endl;

      auto start_time = Clock::now();

      //update the hybrid measurement clusters from positional info and construct down messages
      for(uint o = 0; o<objects; o++){
        uint oBelIndex = tBelIndex+(1+measurements)*o;
        uint oMesIndex = tMesIndex+measurements*o;
        uint oCID = tCID+stateVectorSize*o;
        uint sepBelIndexOffset = 0;
        if(t == timesteps-2){sepBelIndexOffset = o;}

        for(uint m = 0; m<measurements;m++){
          uint DID = tDID+m*objects+o;
          oldSepBel = sepsetBeliefs[oBelIndex+m-sepBelIndexOffset];
          sepsetBeliefs[oBelIndex+m-sepBelIndexOffset] = clustBeliefs[oBelIndex]->marginalize({oCID,oCID+1});
          clustBeliefs[oBelIndex+1+m] = clustBeliefs[oBelIndex+1+m]->absorb(sepsetBeliefs[oBelIndex+m-sepBelIndexOffset])->cancel(oldSepBel);
          interfaceDownMessages[oMesIndex+m] = clustBeliefs[oBelIndex+1+m]->marginalize({DID})->cancel(interfaceUpMessages[oMesIndex+m]); //µ′a,b = Sum:Ψ /µb,a
        }
      }
      if(OL >= Verbose) cout << "updated hybrid meas clusters and down interface messages" << endl;

      //run discrete belief update, extract up messages, update the hybrid measurement clusters from discrete inference info
      vector<rcptr<Factor>> DGBeliefs = discreteGraph(DiscreteBase, t, objects, measurements, interfaceDownMessages);
      for(uint o = 0; o < objects; o++){
        uint oMesIndex = o*measurements;
        uint oBelIndex = tBelIndex+(1+measurements)*o;
        if(OL >= Periphrastic) cout << "discrete extract " << o << endl;
        for(uint m = 0; m < measurements; m++){
          uint mMesIndex = oMesIndex+m;
          uint mBelIndex = oBelIndex+1+m; //+1 included for pos cluster offset
          newUpMessages[mMesIndex] = DGBeliefs[oMesIndex+m]->cancel(interfaceDownMessages[tMesIndex+oMesIndex+m]);//µ′a,b = Sum:Ψ /µb,a
          if(OL >= Periphrastic) cout << "discrete up message \n" << *newUpMessages[mMesIndex] << endl;
          clustBeliefs[mBelIndex] = clustBeliefs[mBelIndex]->absorb(newUpMessages[mMesIndex])->cancel(interfaceUpMessages[tMesIndex+mMesIndex]); //Ψ′ = Ψ(µ′a,b/µa,b)
          interfaceUpMessages[tMesIndex+mMesIndex] = newUpMessages[mMesIndex]; //µ=µ′
        }
      }
      if(OL >= Verbose) cout << "updated discrete system and up interface messages" << endl;
      if(OL >= Verbose) cout << "updated hybrid meas clusters from discrete messages" << endl;

      //update current pos clusters
      for(uint o = 0; o < objects; o++){
        uint oBelIndex = tBelIndex+(1+measurements)*o;
        uint oCID = tCID+stateVectorSize*o;
        uint sepBelIndexOffset = 0;
        if(t == timesteps-2){sepBelIndexOffset = o;}
        for(uint m = 0; m<measurements;m++){ //TODO check the values here
          uint mBelIndex = oBelIndex+m;
          oldSepBel = sepsetBeliefs[mBelIndex-sepBelIndexOffset];
          sepsetBeliefs[mBelIndex-sepBelIndexOffset] = clustBeliefs[mBelIndex+1]->marginalize({oCID,oCID+1});
          clustBeliefs[oBelIndex] = clustBeliefs[oBelIndex]->absorb(sepsetBeliefs[mBelIndex-sepBelIndexOffset])->cancel(oldSepBel);
        }
      }
      if(OL >= Verbose) cout << "updated current pos cluster from measurements" << endl;

      //update previous pos clusters
      for(uint o = 0; o < objects; o++){
        uint oBelIndex = tBelIndex+(1+measurements)*o;
        uint PoBelIndex = (t-1)*(1+measurements)*objects+(1+measurements)*o; //TODO check
        uint oCID = tCID+stateVectorSize*o;
        oldSepBel = sepsetBeliefs[PoBelIndex+measurements];
        sepsetBeliefs[PoBelIndex+measurements] = clustBeliefs[oBelIndex]->marginalize({oCID,oCID+1,oCID+2,oCID+3,oCID+4,oCID+5});
        clustBeliefs[PoBelIndex] = clustBeliefs[PoBelIndex]->absorb(sepsetBeliefs[PoBelIndex+measurements])->cancel(oldSepBel);
      }
      if(OL >= Verbose) cout << "updated previous pos clusters" << endl;

      if(OL >= Verbose){
      auto end_time = Clock::now();
      std::cout << "Time difference:"
          << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << std::endl;
      }
    }

    {//scope for first clusters
      uint tBelIndex = 0;
      uint tMesIndex = 0;
      uint tCID = 0;
      uint tDID = DiscreteBase;

      //update the first hybrid measurement clusters from positional info
      for(uint o = 0; o<objects; o++){
        uint oBelIndex = tBelIndex+(1+measurements)*o;
        uint oMesIndex = tMesIndex+measurements*o;
        uint oCID = tCID+stateVectorSize*o;
        for(uint m = 0; m<measurements;m++){
          uint DID = tDID+m*objects+o;
          oldSepBel = sepsetBeliefs[oBelIndex+m];
          sepsetBeliefs[oBelIndex+m] = clustBeliefs[oBelIndex]->marginalize({oCID,oCID+1});
          clustBeliefs[oBelIndex+1+m] = clustBeliefs[oBelIndex+1+m]->absorb(sepsetBeliefs[oBelIndex+m])->cancel(oldSepBel);
          interfaceDownMessages[oMesIndex+m] = clustBeliefs[oBelIndex+1+m]->marginalize({DID})->cancel(interfaceUpMessages[oMesIndex+m]); //µ′a,b = Sum:Ψ /µb,a
        }
      }
      if(OL >= Verbose) cout << "updated hybrid meas clusters and down interface messages" << endl;

      //TODO this code is unchecked
      //run discrete belief update and update the first hybrid measurement clusters from discrete inference info
      vector<rcptr<Factor>> DGBeliefs = discreteGraph(DiscreteBase, 0, objects, measurements, interfaceDownMessages);
      for(uint o = 0; o < objects; o++){
        uint oMesIndex = o*measurements;
        uint oBelIndex = tBelIndex+(1+measurements)*o;
        if(OL >= Periphrastic) cout << "discrete extract " << o << endl;
        for(uint m = 0; m < measurements; m++){
          uint mMesIndex = oMesIndex+m;
          uint mBelIndex = oBelIndex+1+m; //+1 included for pos cluster offset
          newUpMessages[mMesIndex] = DGBeliefs[oMesIndex+m]->cancel(interfaceDownMessages[tMesIndex+oMesIndex+m]);//µ′a,b = Sum:Ψ /µb,a
          if(OL >= Periphrastic) cout << "discrete up message \n" << *newUpMessages[mMesIndex] << endl;
          clustBeliefs[mBelIndex] = clustBeliefs[mBelIndex]->absorb(newUpMessages[mMesIndex])->cancel(interfaceUpMessages[tMesIndex+mMesIndex]); //Ψ′ = Ψ(µ′a,b/µa,b)
          interfaceUpMessages[tMesIndex+mMesIndex] = newUpMessages[mMesIndex]; //µ=µ′
        }
      }
      if(OL >= Verbose) cout << "updated discrete system and up interface messages" << endl;
      if(OL >= Verbose) cout << "updated hybrid meas clusters from discrete messages" << endl;

      //update first pos clusters
      for(uint o = 0; o < objects; o++){
        uint oBelIndex = tBelIndex+(1+measurements)*o;
        uint oCID = tCID+stateVectorSize*o;
        for(uint m = 0; m<measurements;m++){
          uint mBelIndex = oBelIndex+m;
          oldSepBel = sepsetBeliefs[mBelIndex];
          sepsetBeliefs[mBelIndex] = clustBeliefs[mBelIndex+1]->marginalize({oCID,oCID+1});
          clustBeliefs[oBelIndex] = clustBeliefs[oBelIndex]->absorb(sepsetBeliefs[mBelIndex])->cancel(oldSepBel);
        }
      }
      if(OL >= Verbose) cout << "updated current pos cluster from measurements" << endl;

      if(OL >= Verbose) cout << "completed back messages" << endl;
    }//end scope for first clusters

    uint loopMax = 1;
    if(atoi(argv[1]) > 0) {loopMax = atoi(argv[1]);}
    else {if(OL >= Error) cout << "please enter a valid input argument" << endl;}

    // for(uint loop = 2; loop <= loopMax; loop++){

    //   if(OL >= Verbose) cout << "\nloop " << loop << '\n' << endl;

    //   //send messages forwards
    //   clustBeliefs[0] = clustBeliefs[0]->absorb(expect)->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5,bID+6,bID+7,bID+8,bID+9,bID+10,bID+11});
    //   for(uint t = 0; t<timesteps-2;t++){
    //     uint bID = stateVectorSize*(t);
    //     uint belIndex = t*(1+measurements);
    //     uint mesTIndex = t*measurements*objects;
    //     uint timeBID = DiscreteBase+t*objects*measurements;

    //     if(OL >= Verbose) cout << "\nforward loop" << t << endl;

    //     //update the hybrid measurement clusters from positional info
    //     for(uint m = 0; m<measurements;m++){
    //       uint measBID = timeBID+objects*m;
    //       oldSepBel = sepsetBeliefs[belIndex+m];
    //       sepsetBeliefs[belIndex+m] = clustBeliefs[belIndex]->marginalize({bID,bID+1});
    //       clustBeliefs[belIndex+1+m] = clustBeliefs[belIndex+1+m]->absorb(sepsetBeliefs[belIndex+m])->cancel(oldSepBel);
    //       interfaceDownMessages[mesTIndex+m*objects+0] = clustBeliefs[belIndex+1+m]->marginalize({measBID})->cancel(interfaceUpMessages[mesTIndex+m*objects+0]); //µ′a,b = Sum:Ψ /µb,a
    //     }
    //     if(OL >= Verbose) cout << "updated hybrid meas clusters and down interface messages" << endl;

    //     //run discrete belief update
    //     vector<rcptr<Factor>> DGBeliefs = discreteGraph(DiscreteBase, t, objects, measurements, interfaceDownMessages);
    //     for(uint m = 0; m < measurements; m++){
    //       uint mesMIndex = m*objects;
    //       for(uint o = 0; o < objects; o++){ //TODO check the indicies here
    //         newUpMessages[mesMIndex+o] = DGBeliefs[mesMIndex+o]->cancel(interfaceDownMessages[mesTIndex+mesMIndex+o]);//µ′a,b = Sum:Ψ /µb,a
    //       }
    //     }
    //     if(OL >= Verbose) cout << "updated discrete system and up interface messages" << endl;

    //     //update the hybrid measurement clusters from discrete inference info
    //     for(uint m = 0; m < measurements; m++){
    //       uint mesMIndex = m*objects;
    //       for(uint o = 0; o < objects; o++){
    //         clustBeliefs[belIndex+1+mesMIndex+o] = clustBeliefs[belIndex+1+mesMIndex+o]->absorb(newUpMessages[mesMIndex+o])->cancel(interfaceUpMessages[mesTIndex+mesMIndex+o]); //Ψ′ = Ψ(µ′a,b/µa,b)
    //         interfaceUpMessages[mesTIndex+mesMIndex+o] = newUpMessages[mesMIndex+o]; //µ=µ′
    //       }
    //     }
    //     if(OL >= Verbose) cout << "updated hybrid meas clusters from discrete messages" << endl;

    //     //update current pos clusters
    //     for(uint m = 0; m<measurements;m++){
    //       oldSepBel = sepsetBeliefs[belIndex+m];
    //       sepsetBeliefs[belIndex+m] = clustBeliefs[belIndex+1+m]->marginalize({bID,bID+1});
    //       clustBeliefs[belIndex] = clustBeliefs[belIndex]->absorb(sepsetBeliefs[belIndex+m])->cancel(oldSepBel);
    //     }
    //     if(OL >= Verbose) cout << "updated current pos cluster from measurements" << endl;

    //     //update next pos clusters
    //     oldSepBel = sepsetBeliefs[belIndex+measurements];
    //     sepsetBeliefs[belIndex+measurements] = clustBeliefs[belIndex]->marginalize({bID+6,bID+7,bID+8,bID+9,bID+10,bID+11});
    //     clustBeliefs[belIndex+measurements+1] = clustBeliefs[belIndex+measurements+1]->absorb(sepsetBeliefs[belIndex+measurements])->cancel(oldSepBel);
    //     if(OL >= Verbose) cout << "updated next pos cluster" << endl;
    //   }

    //   {//scope for final clusters
    //     bID = stateVectorSize*(timesteps-2);
    //     uint timeBID = DiscreteBase+(timesteps-2)*measurements*objects;
    //     uint mesTIndex = (timesteps-2)*measurements*objects;
    //     uint belIndex = (timesteps-2)*(1+measurements);

    //     //update the second last hybrid measurement clusters from positional info
    //     for(uint m = 0; m<measurements;m++){
    //       uint measBID = timeBID+objects*m+0;
    //       oldSepBel = sepsetBeliefs[belIndex+m];
    //       sepsetBeliefs[belIndex+m] = clustBeliefs[belIndex]->marginalize({bID,bID+1});
    //       clustBeliefs[belIndex+1+m] = clustBeliefs[belIndex+1+m]->absorb(sepsetBeliefs[belIndex+m])->cancel(oldSepBel);
    //       interfaceDownMessages[mesTIndex+m*objects+0] = clustBeliefs[belIndex+1+m]->marginalize({measBID})->cancel(interfaceUpMessages[mesTIndex+m*objects+0]); //µ′a,b = Sum:Ψ /µb,a
    //     }
    //     if(OL >= Verbose) cout << "updated hybrid meas clusters and down interface messages" << endl;

    //     //run discrete belief update on second last hybrid measurement clusters
    //     vector<rcptr<Factor>> DGBeliefs = discreteGraph(DiscreteBase, (timesteps-2), objects, measurements, interfaceDownMessages);
    //     for(uint m = 0; m < measurements; m++){
    //       uint mesMIndex = m*objects;
    //       for(uint o = 0; o < objects; o++){ //TODO check the indicies here
    //         newUpMessages[mesMIndex+o] = DGBeliefs[mesMIndex+o]->cancel(interfaceDownMessages[mesTIndex+mesMIndex+o]);//µ′a,b = Sum:Ψ /µb,a
    //       }
    //     }
    //     if(OL >= Verbose) cout << "updated discrete system and up interface messages" << endl;

    //     //update the hybrid measurement clusters from discrete inference info
    //     for(uint m = 0; m < measurements; m++){
    //       uint mesMIndex = m*objects;
    //       for(uint o = 0; o < objects; o++){
    //         clustBeliefs[belIndex+1+mesMIndex+o] = clustBeliefs[belIndex+1+mesMIndex+o]->absorb(newUpMessages[mesMIndex+o])->cancel(interfaceUpMessages[mesTIndex+mesMIndex+o]); //Ψ′ = Ψ(µ′a,b/µa,b)
    //         interfaceUpMessages[mesTIndex+mesMIndex+o] = newUpMessages[mesMIndex+o]; //µ=µ′
    //       }
    //     }
    //     if(OL >= Verbose) cout << "updated hybrid meas clusters from discrete messages" << endl;

    //     //update last pos clusters
    //     for(uint m = 0; m<measurements;m++){
    //       oldSepBel = sepsetBeliefs[belIndex+m];
    //       sepsetBeliefs[belIndex+m] = clustBeliefs[belIndex+1+m]->marginalize({bID,bID+1});
    //       clustBeliefs[belIndex] = clustBeliefs[belIndex]->absorb(sepsetBeliefs[belIndex+m])->cancel(oldSepBel);
    //     }
    //     if(OL >= Verbose) cout << "updated last pos cluster from second last measurements" << endl;

    //     bID = stateVectorSize*(timesteps-1);
    //     timeBID = DiscreteBase+(timesteps-1)*measurements*objects;
    //     mesTIndex = (timesteps-1)*measurements*objects;

    //     //update the last hybrid measurement clusters from positional info
    //     for(uint m = 0; m<measurements;m++){
    //       uint measBID = timeBID+objects*m+0;
    //       oldSepBel = sepsetBeliefs[belIndex+measurements+m];
    //       sepsetBeliefs[belIndex+measurements+m] = clustBeliefs[belIndex]->marginalize({bID,bID+1});
    //       clustBeliefs[belIndex+1+measurements+m] = clustBeliefs[belIndex+1+measurements+m]->absorb(sepsetBeliefs[belIndex+measurements+m])->cancel(oldSepBel);
    //       interfaceDownMessages[mesTIndex+m*objects+0] = clustBeliefs[belIndex+1+measurements+m]->marginalize({measBID})->cancel(interfaceUpMessages[mesTIndex+m*objects+0]); //µ′a,b = Sum:Ψ /µb,a
    //     }
    //     if(OL >= Verbose) cout << "updated hybrid meas clusters and down interface messages" << endl;

    //     //run discrete belief update on last hybrid measurement clusters
    //     DGBeliefs = discreteGraph(DiscreteBase, (timesteps-1), objects, measurements, interfaceDownMessages);
    //     for(uint m = 0; m < measurements; m++){
    //       uint mesMIndex = m*objects;
    //       for(uint o = 0; o < objects; o++){ //TODO check the indicies here
    //         newUpMessages[mesMIndex+o] = DGBeliefs[mesMIndex+o]->cancel(interfaceDownMessages[mesTIndex+mesMIndex+o]);//µ′a,b = Sum:Ψ /µb,a
    //       }
    //     }
    //     if(OL >= Verbose) cout << "updated discrete system and up interface messages" << endl;

    //     //update the last hybrid measurement clusters from discrete inference info
    //     for(uint m = 0; m < measurements; m++){
    //       uint mesMIndex = m*objects;
    //       for(uint o = 0; o < objects; o++){
    //         clustBeliefs[belIndex+1+measurements+mesMIndex+o] = clustBeliefs[belIndex+1+measurements+mesMIndex+o]->absorb(newUpMessages[mesMIndex+o])->cancel(interfaceUpMessages[mesTIndex+mesMIndex+o]); //Ψ′ = Ψ(µ′a,b/µa,b) TODO cancel index check, should be fine though
    //         interfaceUpMessages[mesTIndex+mesMIndex+o] = newUpMessages[mesMIndex+o]; //µ=µ′
    //       }
    //     }
    //     if(OL >= Verbose) cout << "updated hybrid meas clusters from discrete messages" << endl;

    //     //update last pos clusters
    //     for(uint m = 0; m<measurements;m++){
    //       oldSepBel = sepsetBeliefs[belIndex+measurements+m];
    //       sepsetBeliefs[belIndex+measurements+m] = clustBeliefs[belIndex+1+measurements+m]->marginalize({bID,bID+1});
    //       clustBeliefs[belIndex] = clustBeliefs[belIndex]->absorb(sepsetBeliefs[belIndex+measurements+m])->cancel(oldSepBel);
    //     }
    //     if(OL >= Verbose) cout << "updated last pos cluster from last measurements" << endl;

    //     if(OL >= Verbose) cout << "completed forward messages" << endl;
    //   }//end scope for final clusters

    //   //send messages backwards
    //   for(int t = timesteps-3; t>0 ;t--){ //-1 from 0 count, -1 from clust amounts, -1 from sep amount
    //     uint bID = stateVectorSize*(t);
    //     uint belIndex = t*(1+measurements);
    //     uint mesTIndex = t*measurements*objects;
    //     uint timeBID = DiscreteBase+t*objects*measurements;

    //     if(OL >= Verbose) cout << "\nbackward loop " << t << endl;

    //     //update the hybrid measurement clusters from positional info
    //     for(uint m = 0; m<measurements;m++){
    //       uint measBID = timeBID+objects*m;
    //       oldSepBel = sepsetBeliefs[belIndex+m];
    //       sepsetBeliefs[belIndex+m] = clustBeliefs[belIndex]->marginalize({bID,bID+1});
    //       clustBeliefs[belIndex+1+m] = clustBeliefs[belIndex+1+m]->absorb(sepsetBeliefs[belIndex+m])->cancel(oldSepBel);
    //       interfaceDownMessages[mesTIndex+m*objects+0] = clustBeliefs[belIndex+1+m]->marginalize({measBID})->cancel(interfaceUpMessages[mesTIndex+m*objects+0]); //µ′a,b = Sum:Ψ /µb,a
    //     }
    //     if(OL >= Verbose) cout << "updated hybrid meas clusters and down interface messages" << endl;

    //     //run discrete belief update
    //     vector<rcptr<Factor>> DGBeliefs = discreteGraph(DiscreteBase, t, objects, measurements, interfaceDownMessages);
    //     for(uint m = 0; m < measurements; m++){
    //       uint mesMIndex = m*objects;
    //       for(uint o = 0; o < objects; o++){ //TODO check the indicies here
    //         newUpMessages[mesMIndex+o] = DGBeliefs[mesMIndex+o]->cancel(interfaceDownMessages[mesTIndex+mesMIndex+o]);//µ′a,b = Sum:Ψ /µb,a
    //       }
    //     }
    //     if(OL >= Verbose) cout << "updated discrete system and up interface messages" << endl;

    //     //update the hybrid measurement clusters from discrete inference info
    //     for(uint m = 0; m < measurements; m++){
    //       uint mesMIndex = m*objects;
    //       for(uint o = 0; o < objects; o++){
    //         clustBeliefs[belIndex+1+mesMIndex+o] = clustBeliefs[belIndex+1+mesMIndex+o]->absorb(newUpMessages[mesMIndex+o])->cancel(interfaceUpMessages[mesTIndex+mesMIndex+o]); //Ψ′ = Ψ(µ′a,b/µa,b)
    //         interfaceUpMessages[mesTIndex+mesMIndex+o] = newUpMessages[mesMIndex+o]; //µ=µ′
    //       }
    //     }
    //     if(OL >= Verbose) cout << "updated hybrid meas clusters from discrete messages" << endl;

    //     //update current pos clusters
    //     for(uint m = 0; m<measurements;m++){
    //       oldSepBel = sepsetBeliefs[belIndex+m];
    //       sepsetBeliefs[belIndex+m] = clustBeliefs[belIndex+1+m]->marginalize({bID,bID+1});
    //       clustBeliefs[belIndex] = clustBeliefs[belIndex]->absorb(sepsetBeliefs[belIndex+m])->cancel(oldSepBel);
    //     }
    //     if(OL >= Verbose) cout << "updated current pos cluster from measurements" << endl;

    //     //update previous pos clusters
    //     oldSepBel = sepsetBeliefs[belIndex-1];
    //     sepsetBeliefs[belIndex-1] = clustBeliefs[belIndex]->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5});
    //     clustBeliefs[belIndex-(measurements+1)] = clustBeliefs[belIndex-(measurements+1)]->absorb(sepsetBeliefs[belIndex-1])->cancel(oldSepBel);
    //     if(OL >= Verbose) cout << "updated previous pos cluster" << endl;
    //   }

    //   {//scope for first cluster
    //     bID = 0;
    //     uint belIndex = 0;
    //     uint mesTIndex = 0;
    //     uint timeBID = DiscreteBase;

    //     //update the first hybrid measurement clusters from positional info
    //     for(uint m = 0; m<measurements;m++){
    //       uint measBID = timeBID+objects*m;
    //       oldSepBel = sepsetBeliefs[belIndex+m];
    //       sepsetBeliefs[belIndex+m] = clustBeliefs[belIndex]->marginalize({bID,bID+1});
    //       clustBeliefs[belIndex+1+m] = clustBeliefs[belIndex+1+m]->absorb(sepsetBeliefs[belIndex+m])->cancel(oldSepBel);
    //       interfaceDownMessages[mesTIndex+m*objects+0] = clustBeliefs[belIndex+1+m]->marginalize({measBID})->cancel(interfaceUpMessages[mesTIndex+m*objects+0]); //µ′a,b = Sum:Ψ /µb,a
    //     }
    //     if(OL >= Verbose) cout << "updated hybrid meas clusters and down interface messages" << endl;

    //     //run discrete belief update
    //     vector<rcptr<Factor>> DGBeliefs = discreteGraph(DiscreteBase, 0, objects, measurements, interfaceDownMessages);
    //     for(uint m = 0; m < measurements; m++){
    //       uint mesMIndex = m*objects;
    //       for(uint o = 0; o < objects; o++){ //TODO check the indicies here
    //         newUpMessages[mesMIndex+o] = DGBeliefs[mesMIndex+o]->cancel(interfaceDownMessages[mesTIndex+mesMIndex+o]);//µ′a,b = Sum:Ψ /µb,a
    //       }
    //     }
    //     if(OL >= Verbose) cout << "updated discrete system and up interface messages" << endl;

    //     //update the first hybrid measurement clusters from discrete inference info
    //     for(uint m = 0; m < measurements; m++){
    //     uint mesMIndex = m*objects;
    //     for(uint o = 0; o < objects; o++){
    //       clustBeliefs[belIndex+1+mesMIndex+o] = clustBeliefs[belIndex+1+mesMIndex+o]->absorb(newUpMessages[mesMIndex+o])->cancel(interfaceUpMessages[mesTIndex+mesMIndex+o]); //Ψ′ = Ψ(µ′a,b/µa,b)
    //       interfaceUpMessages[mesTIndex+mesMIndex+o] = newUpMessages[mesMIndex+o]; //µ=µ′
    //     }
    //   }
    //     if(OL >= Verbose) cout << "updated hybrid meas clusters from discrete messages" << endl;

    //     //update first pos clusters
    //     for(uint m = 0; m<measurements;m++){
    //       oldSepBel = sepsetBeliefs[belIndex+m];
    //       sepsetBeliefs[belIndex+m] = clustBeliefs[belIndex+1+m]->marginalize({bID,bID+1});
    //       clustBeliefs[belIndex] = clustBeliefs[belIndex]->absorb(sepsetBeliefs[belIndex+m])->cancel(oldSepBel);
    //     }
    //     if(OL >= Verbose) cout << "updated current pos cluster from measurements" << endl;

    //     if(OL >= Verbose) cout << "completed back messages" << endl;
    //   }//end scope for first cluster
    // }

    //NOTE: atm only utilising diagonal covariance elements as they aren't correlated, but this will likely not be the case with more complex model

    ofstream xPred("xPred.txt");
    ofstream yPred("yPred.txt");
    ofstream xCovs("xCovPred.txt");
    ofstream yCovs("yCovPred.txt");
    for(uint o = 0; o<objects;o++){
      for(uint t = 0; t<timesteps-1;t++){
          uint CID = t*stateVectorSize*objects+stateVectorSize*o;
          uint belIndex = t*(1+measurements)*objects+(1+measurements)*o;

          rcptr<Factor> qPtr = clustBeliefs[belIndex]->marginalize({CID})->normalize();
          xPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0] << ',';

          qPtr = clustBeliefs[belIndex]->marginalize({CID+1})->normalize();
          yPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0] << ',';

          qPtr = clustBeliefs[belIndex]->marginalize({CID,CID+1})->normalize();
          xCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[0][0] << ',';
          yCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[1][1] << ',';
      }

      uint CID = (timesteps-1)*stateVectorSize*objects+stateVectorSize*o;
      uint belIndex = (timesteps-2)*(1+measurements)*objects+(1+measurements)*o;

      rcptr<Factor> qPtr = clustBeliefs[belIndex]->marginalize({CID})->normalize();
      xPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0];

      qPtr = clustBeliefs[belIndex]->marginalize({CID+1})->normalize();
      yPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0];

      qPtr = clustBeliefs[belIndex]->marginalize({CID,CID+1})->normalize();
      xCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[0][0];
      yCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[1][1];

      xPred << "\n";
      yPred << "\n";
      xCovs << "\n";
      yCovs << "\n";

    }


    xPred.close();
    yPred.close();
    xCovs.close();
    yCovs.close();

    return 0;
  } // try

  catch (char msg[]) {
    cerr << msg << endl;
  } // catch

  // catch (char const* msg) {
  //   cerr << msg << endl;
  // } // catch

  catch (const string& msg) {
    cerr << msg << endl;
    throw;
  } // catch

  catch (const exception& e) {
    cerr << "Unhandled exception: " << e.what() << endl;
    throw e;
  } // catch

  catch(...) {
    cerr << "An unknown exception / error occurred\n";
    throw;
  } // catch

} // main

vector<rcptr<Factor>> discreteGraph(uint discreteBase, uint timestep, uint objects, uint measurements, vector<rcptr<Factor>> messages){
  // cout << "entered DG" << endl;
  double defProb = 0.0;
  rcptr<vector<T>> binDom(new vector<T>{0,1});
  vector<rcptr<Factor>> discreteFactorPtrs; //technically it is possible to calculate the size here already and allocate memory accordingly instead of using push methods as this is much more effecient
  vector<rcptr<Factor>> returnMessages(measurements*objects);

  uint timeBID = discreteBase+timestep*objects*measurements;
  // cout << "\time: " << t << endl;
  for(uint m = 0; m < measurements; m++){
      uint measBID = timeBID+objects*m;
      // cout << "\nmeas: " << m << endl;
      for(uint o = 0; o < objects; o++){
          // cout << "\nobj: " << m << endl;
          for(uint oo = o+1; oo < objects; oo++){
              // cout << "\nobj-obj: " << m << endl;
              // cout << "DID " << measBID+o << "-" << measBID+oo << endl;
              discreteFactorPtrs.push_back( uniqptr<DT>( new DT({measBID+o,measBID+oo}, {binDom,binDom}, defProb,
              {{{0,0}, 1.0},
              {{0,1}, 1.0},
              {{1,0}, 1.0},} )));
          }
          for(uint mm = m+1; mm < measurements; mm++){
              // cout << "\nmeas-meas: " << m << endl;
              discreteFactorPtrs.push_back( uniqptr<DT>( new DT({measBID+o,timeBID+objects*mm+o}, {binDom,binDom}, defProb,
              {{{0,0}, 1.0},
              {{0,1}, 1.0},
              {{1,0}, 1.0},} )));
          }
      }
  }

  uint tMesIndex = timestep*measurements*objects;
  for(uint o = 0; o < objects; o++){
    uint oMesIndex = tMesIndex+measurements*o;
    for(uint m = 0; m < measurements; m++){
      // cout << "object" << o << "\nmeas" << m << "\n\n\nindex" << oMesIndex+m << "\nmessage" << *messages[oMesIndex+m] << endl;
      discreteFactorPtrs.push_back(messages[oMesIndex+m]);
    }
  }

  // cout << "finished constructing DG" << endl;

  map<RVIdType, AnyType> obsv;

  // for a factor graph, use BETHE instead of LTRIP
  // for a junction tree, use JTREE instead of LTRIP
  ClusterGraph cg(ClusterGraph::LTRIP, discreteFactorPtrs, obsv);
  //cout << cg << endl;

  // cout << "finished constructing CG" << endl;

  // export the graph to graphviz .dot format
  // cg.exportToGraphViz("testGraph");

  map<Idx2, rcptr<Factor> > msgs;
  MessageQueue msgQ;

  unsigned nMsgs = loopyBU_CG(cg, msgs, msgQ);
  // cout << "Sent " << nMsgs << " messages before convergence\n";

  // cout << "finished LBU" << endl;

  for(uint o = 0; o < objects; o++){
    uint oMesIndex = measurements*o;
    for(uint m = 0; m < measurements; m++){
      uint DID = timeBID+objects*m+o;
      // cout << "object " << o << "\nmeas " << m << "\nindex " << oMesIndex+m << "\nDID " << DID << endl;
      rcptr<Factor> qPtr = queryLBU_CG(cg, msgs, {DID})->normalize();
      // cout << "object " << o << "\nmeas " << m << "\nindex " << oMesIndex+m << "\nRet message " << *qPtr << endl;
      returnMessages[oMesIndex+m] = qPtr;
      // cout << "object " << o << "\nmeas " << m << "\nindex " << oMesIndex+m << "\nRet message " << *qPtr << endl;
    }
  }

  return returnMessages;
}