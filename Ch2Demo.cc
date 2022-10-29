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

using namespace std;
using namespace emdw;

typedef SqrtMVG SG;
typedef AffineMVG AG;
typedef unsigned T;



int main(int argc, char *argv[]) {
  try {

    enum outputLevel {None, Error, Debug, Verbose, Periphrastic};
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

    vector<double> measXList;
    vector<double> measYList;

    char c;
    std::string str;

    ifstream inFile;

    inFile.open("Xpos.txt");
    if (!inFile) {
        if(OL >= Error) cout << "Unable to open file";
        exit(1); // terminate with error
    }

    while (inFile >> c) {
      switch(c) {
        case '[':
          if(OL >= Verbose) cout << "start x parse" << endl;
          break;
        case ',':
        case ']':
          measXList.push_back(std::stod(str));
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

    while (inFile >> c) {
      switch(c) {
        case '[':
          if(OL >= Verbose) cout << "start Y parse" << endl;
          break;
        case ',':
        case ']':
          measYList.push_back(std::stod(str));
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

    uint timesteps = measXList.size();
    if (timesteps < 2){
      if(OL >= Error) cout << "At least 2 timesteps required";
        exit(1); // terminate with error
    }
    if (timesteps != measYList.size()){
      if(OL >= Error) cout << "measurement dimension mismatch";
        exit(1); // terminate with error
    }

    vector<rcptr<Factor>> posPdfs;
    vector<rcptr<Factor>> measHdfs;

    uint stateVectorSize = 6;

    //Discrete setup:
    uint DiscreteBase = stateVectorSize*timesteps;
    rcptr< vector<T> > binDom( new vector<T>{0,1} );
    double clutterProb = 0.01;
    prlite::ColVector<double> uninformMn(2);
    prlite::RowMatrix<double> uninforCov(2,2);
    uninformMn[0] = 0; uninformMn[1] = 0;
    uninforCov(0,0) = 1000.0; uninforCov(0,1) = 0.0;
    uninforCov(1,0) = 0.0; uninforCov(1,1) = 1000.0;

    //construct the cluster graph
    for(uint i = 0; i<timesteps-1;i++){

        uint bID = stateVectorSize*i;

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

        posPdfs.push_back(uniqptr<AG>(new AG({bID,bID+1,bID+2,bID+3,bID+4,bID+5},{bID+6,bID+7,bID+8,bID+9,bID+10,bID+11}, A, C, R)));

        prlite::ColVector<double> mn(2);
        prlite::RowMatrix<double> cov(2,2);
        mn[0] = measXList[i]; mn[1] = measYList[i];
        cov(0,0) = 10.0; cov(0,1) = 0.0;
        cov(1,0) = 0.0; cov(1,1) = 10.0;
        map< vector<T>, tuple< double, prlite::ColVector<double>, prlite::RowMatrix<double> > > weightedComponents =
        {{{0}, { clutterProb, uninformMn, uninforCov}},
        {{1}, { 1.0-clutterProb, mn, cov}}};
        measHdfs.push_back(uniqptr<SqrtMVGTable<T>>( new SqrtMVGTable<T>({DiscreteBase+i}, {{bID,bID+1}}, {binDom}, weightedComponents)));
    }

    //constructs the final measurements mvg as there is one more than the number of clusters
    uint bID = stateVectorSize*(timesteps-1);
    prlite::ColVector<double> mn(2);
    prlite::RowMatrix<double> cov(2,2);
    mn[0] = measXList[timesteps-1]; mn[1] = measYList[timesteps-1];
    cov(0,0) = 10.0; cov(0,1) = 0.0;
    cov(1,0) = 0.0; cov(1,1) = 10.0;
    map< vector<T>, tuple< double, prlite::ColVector<double>, prlite::RowMatrix<double> > > weightedComponents =
    {{{0}, { clutterProb, uninformMn, uninforCov}},
    {{1}, { 1.0-clutterProb, mn, cov}}};
    measHdfs.push_back(uniqptr<SqrtMVGTable<T>>( new SqrtMVGTable<T>({DiscreteBase+timesteps-1}, {{bID,bID+1}}, {binDom}, weightedComponents)));

    //construct general expectations of state variables
    bID = 0;
    prlite::ColVector<double> mnEx(6);
    prlite::RowMatrix<double> covEx(6,6);
    mnEx[0] = 50; mnEx[1] = 50; mnEx[2] = 0; mnEx[3] = 0; mnEx[4] = 0; mnEx[5] = 0;
    covEx(0,0) = 100.0; covEx(0,1) = 0.0; covEx(0,2) = 0.0; covEx(0,3) = 0.0; covEx(0,4) = 0.0; covEx(0,5) = 0.0;
    covEx(1,0) = 0.0; covEx(1,1) = 100.0; covEx(1,2) = 0.0; covEx(1,3) = 0.0; covEx(1,4) = 0.0; covEx(1,5) = 0.0;
    covEx(2,0) = 0.0; covEx(2,1) = 0.0; covEx(2,2) = 100.0; covEx(2,3) = 0.0; covEx(2,4) = 0.0; covEx(2,5) = 0.0;
    covEx(3,0) = 0.0; covEx(3,1) = 0.0; covEx(3,2) = 0.0; covEx(3,3) = 100.0; covEx(3,4) = 0.0; covEx(3,5) = 0.0;
    covEx(4,0) = 0.0; covEx(4,1) = 0.0; covEx(4,2) = 0.0; covEx(4,3) = 0.0; covEx(4,4) = 100.0; covEx(4,5) = 0.0;
    covEx(5,0) = 0.0; covEx(5,1) = 0.0; covEx(5,2) = 0.0; covEx(5,3) = 0.0; covEx(5,4) = 0.0; covEx(5,5) = 100.0;
    rcptr<Factor> expect = uniqptr<SG>(new SG({bID,bID+1,bID+2,bID+3,bID+4,bID+5}, mnEx, covEx));


    vector<rcptr<Factor>> clustBeliefs(2*timesteps-1);
    vector<rcptr<Factor>> sepsetBeliefs(2*timesteps-2); // (t-2)+t posSepset+measSepset
    rcptr<Factor> oldSepBel;

    // construct cluster beliefs
    for(uint i = 0; i<timesteps-1;i++){
      clustBeliefs[i*2] = posPdfs[i];
      clustBeliefs[i*2+1] = measHdfs[i];
    }
    clustBeliefs[2*(timesteps-2)+2] = measHdfs[timesteps-1];

    if(OL >= Verbose) cout << "completed construction" << endl;

    //send messages forwards
    clustBeliefs[0] = clustBeliefs[0]->absorb(expect)->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5,bID+6,bID+7,bID+8,bID+9,bID+10,bID+11});
    for(uint i = 0; i<timesteps-2;i++){
      uint bID = stateVectorSize*(i);

      if(OL >= Verbose) cout << "\nforward loop" << i << endl;

      //update measurement clusters
      sepsetBeliefs[2*i] = clustBeliefs[2*i]->marginalize({bID,bID+1});
      clustBeliefs[2*i+1] = clustBeliefs[2*i+1]->absorb(sepsetBeliefs[2*i]);
      if(OL >= Verbose) cout << "updated meas cluster" << endl;

      //update current pos clusters
      oldSepBel = sepsetBeliefs[2*i];
      sepsetBeliefs[2*i] = clustBeliefs[2*i+1]->marginalize({bID,bID+1});
      clustBeliefs[2*i] = clustBeliefs[2*i]->absorb(sepsetBeliefs[2*i])->cancel(oldSepBel);
      if(OL >= Verbose) cout << "updated current pos cluster" << endl;

      //update next pos clusters
      sepsetBeliefs[2*i+1] = clustBeliefs[2*i]->marginalize({bID+6,bID+7,bID+8,bID+9,bID+10,bID+11});
      clustBeliefs[2*i+2] = clustBeliefs[2*i+2]->absorb(sepsetBeliefs[2*i+1])->marginalize({bID+6,bID+7,bID+8,bID+9,bID+10,bID+11,bID+12,bID+13,bID+14,bID+15,bID+16,bID+17});
      if(OL >= Verbose) cout << "updated next pos cluster" << endl;
    }

    bID = stateVectorSize*(timesteps-2);

    //update second last measurement cluster
    sepsetBeliefs[2*(timesteps-2)] = clustBeliefs[2*(timesteps-2)]->marginalize({bID,bID+1});
    clustBeliefs[2*(timesteps-2)+1] = clustBeliefs[2*(timesteps-2)+1]->absorb(sepsetBeliefs[2*(timesteps-2)]);
    if(OL >= Verbose) cout << "updated second last meas cluster" << endl;

    //update last pos clusters
    oldSepBel = sepsetBeliefs[2*(timesteps-2)];
    sepsetBeliefs[2*(timesteps-2)] = clustBeliefs[2*(timesteps-2)+1]->marginalize({bID,bID+1});
    clustBeliefs[2*(timesteps-2)] = clustBeliefs[2*(timesteps-2)]->absorb(sepsetBeliefs[2*(timesteps-2)])->cancel(oldSepBel);
    if(OL >= Verbose) cout << "updated last pos cluster" << endl;

    //update last measurement cluster
    sepsetBeliefs[2*(timesteps-2)+1] = clustBeliefs[2*(timesteps-2)]->marginalize({bID+6,bID+7});
    clustBeliefs[2*(timesteps-2)+2] = clustBeliefs[2*(timesteps-2)+2]->absorb(sepsetBeliefs[2*(timesteps-2)+1]);
    if(OL >= Verbose) cout << "updated last meas cluster" << endl;

    //update last pos clusters
    oldSepBel = sepsetBeliefs[2*(timesteps-2)+1];
    sepsetBeliefs[2*(timesteps-2)+1] = clustBeliefs[2*(timesteps-2)+2]->marginalize({bID+6,bID+7});
    clustBeliefs[2*(timesteps-2)] = clustBeliefs[2*(timesteps-2)]->absorb(sepsetBeliefs[2*(timesteps-2)+1])->cancel(oldSepBel);
    if(OL >= Verbose) cout << "updated last pos cluster " << endl;

    if(OL >= Verbose) cout << "completed forward messages" << endl;

    //send messages backwards
    for(int i = timesteps-2; i>0 ;i--){ //-1 from 0 count, -1 from clust amounts, -1 from sep amount
      uint bID = stateVectorSize*(i);

      if(OL >= Verbose) cout << "\nbackward loop " << i << endl;

      //update measurement clusters
      oldSepBel = sepsetBeliefs[2*i];
      sepsetBeliefs[2*i] = clustBeliefs[2*i]->marginalize({bID,bID+1});
      clustBeliefs[2*i+1] = clustBeliefs[2*i+1]->absorb(sepsetBeliefs[2*i])->cancel(oldSepBel);

      //update current pos clusters
      oldSepBel = sepsetBeliefs[2*i];
      sepsetBeliefs[2*i] = clustBeliefs[2*i+1]->marginalize({bID,bID+1});
      clustBeliefs[2*i] = clustBeliefs[2*i]->absorb(sepsetBeliefs[2*i])->cancel(oldSepBel);
      if(OL >= Verbose) cout << "updated current pos cluster" << endl;

      //update prev pos clusters
      oldSepBel = sepsetBeliefs[2*i-1];
      sepsetBeliefs[2*i-1] = clustBeliefs[2*i]->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5});
      clustBeliefs[2*i-2] = clustBeliefs[2*i-2]->absorb(sepsetBeliefs[2*i-1])->cancel(oldSepBel);
      if(OL >= Verbose) cout << "updated next pos cluster" << endl;
    }

    bID = 0;

    //update first measurement clusters
    oldSepBel = sepsetBeliefs[0];
    sepsetBeliefs[0] = clustBeliefs[0]->marginalize({bID,bID+1});
    clustBeliefs[1] = clustBeliefs[1]->absorb(sepsetBeliefs[0])->cancel(oldSepBel);
    if(OL >= Verbose) cout << "updated first meas cluster" << endl;

    //update first pos clusters
    oldSepBel = sepsetBeliefs[0];
    sepsetBeliefs[0] = clustBeliefs[1]->marginalize({bID,bID+1});
    clustBeliefs[1] = clustBeliefs[1]->absorb(sepsetBeliefs[0])->cancel(oldSepBel);
    if(OL >= Verbose) cout << "updated first pos cluster" << endl;

    if(OL >= Verbose) cout << "completed back messages" << endl;

    uint loopMax = 1;
    if(atoi(argv[1]) > 0) {loopMax = atoi(argv[1]);}
    else {if(OL >= Error) cout << "please enter a valid input argument" << endl;}

    for(uint loop = 2; loop <= loopMax; loop++){

      if(OL >= Verbose) cout << "\nloop " << loop << '\n' << endl;

      //send messages forwards
      for(uint i = 0; i<timesteps-2;i++){
        uint bID = stateVectorSize*(i);
        if(OL >= Verbose) cout << "\nforward loop " << i << endl;

        //update measurement clusters
        oldSepBel = sepsetBeliefs[2*i];
        sepsetBeliefs[2*i] = clustBeliefs[2*i]->marginalize({bID,bID+1});
        clustBeliefs[2*i+1] = clustBeliefs[2*i+1]->absorb(sepsetBeliefs[2*i])->cancel(oldSepBel);
        if(OL >= Verbose) cout << "meas updated" << endl;

        //update current pos clusters
        oldSepBel = sepsetBeliefs[2*i];
        sepsetBeliefs[2*i] = clustBeliefs[2*i+1]->marginalize({bID,bID+1});
        clustBeliefs[2*i] = clustBeliefs[2*i]->absorb(sepsetBeliefs[2*i])->cancel(oldSepBel);
        if(OL >= Verbose) cout << "current pos updated" << endl;

        //update next pos clusters
        oldSepBel = sepsetBeliefs[2*i+1];
        sepsetBeliefs[2*i+1] = clustBeliefs[2*i]->marginalize({bID+6,bID+7,bID+8,bID+9,bID+10,bID+11});
        clustBeliefs[2*i+2] = clustBeliefs[2*i+2]->absorb(sepsetBeliefs[2*i+1])->cancel(oldSepBel);
        if(OL >= Verbose) cout << "next pos updated" << endl;
      }

      bID = stateVectorSize*(timesteps-2);

      //update second last measurement cluster
      oldSepBel = sepsetBeliefs[2*(timesteps-2)];
      sepsetBeliefs[2*(timesteps-2)] = clustBeliefs[2*(timesteps-2)]->marginalize({bID,bID+1});
      clustBeliefs[2*(timesteps-2)+1] = clustBeliefs[2*(timesteps-2)+1]->absorb(sepsetBeliefs[2*(timesteps-2)])->cancel(oldSepBel);
      if(OL >= Verbose) cout << "updated second last meas cluster" << endl;

      //update last pos cluster from second last measurment
      oldSepBel = sepsetBeliefs[2*(timesteps-2)];
      sepsetBeliefs[2*(timesteps-2)] = clustBeliefs[2*(timesteps-2)+1]->marginalize({bID,bID+1});
      clustBeliefs[2*(timesteps-2)] = clustBeliefs[2*(timesteps-2)]->absorb(sepsetBeliefs[2*(timesteps-2)])->cancel(oldSepBel);
      if(OL >= Verbose) cout << "updated last pos cluster" << endl;

      //update last measurement cluster
      oldSepBel = sepsetBeliefs[2*(timesteps-2)+1];
      sepsetBeliefs[2*(timesteps-2)+1] = clustBeliefs[2*(timesteps-2)]->marginalize({bID+6,bID+7});
      clustBeliefs[2*(timesteps-2)+2] = clustBeliefs[2*(timesteps-2)+2]->absorb(sepsetBeliefs[2*(timesteps-2)+1])->cancel(oldSepBel);
      if(OL >= Verbose) cout << "updated last meas cluster" << endl;

      //update last pos cluster from last measurment
      oldSepBel = sepsetBeliefs[2*(timesteps-2)+1];
      sepsetBeliefs[2*(timesteps-2)+1] = clustBeliefs[2*(timesteps-2)+2]->marginalize({bID+6,bID+7});
      clustBeliefs[2*(timesteps-2)] = clustBeliefs[2*(timesteps-2)]->absorb(sepsetBeliefs[2*(timesteps-2)+1])->cancel(oldSepBel);
      if(OL >= Verbose) cout << "updated last pos cluster " << endl;

      if(OL >= Verbose) cout << "completed forward messages" << endl;

      //send messages backwards
      for(int i = timesteps-3; i>0 ;i--){ //-1 from 0 count, -1 from clust amounts, -1 from sep amount
        uint bID = stateVectorSize*(i);

        if(OL >= Verbose) cout << "\nbackward loop " << i << endl;

        //update measurement clusters
        oldSepBel = sepsetBeliefs[2*i];
        sepsetBeliefs[2*i] = clustBeliefs[2*i]->marginalize({bID,bID+1});
        clustBeliefs[2*i+1] = clustBeliefs[2*i+1]->absorb(sepsetBeliefs[2*i])->cancel(oldSepBel);
        if(OL >= Verbose) cout << "updated meas cluster" << endl;

        //update current pos clusters
        oldSepBel = sepsetBeliefs[2*i];
        sepsetBeliefs[2*i] = clustBeliefs[2*i+1]->marginalize({bID,bID+1});
        clustBeliefs[2*i] = clustBeliefs[2*i]->absorb(sepsetBeliefs[2*i])->cancel(oldSepBel);
        if(OL >= Verbose) cout << "updated current pos cluster" << endl;

        //update prev pos clusters
        oldSepBel = sepsetBeliefs[2*i-1];
        sepsetBeliefs[2*i-1] = clustBeliefs[2*i]->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5});
        clustBeliefs[2*i-2] = clustBeliefs[2*i-2]->absorb(sepsetBeliefs[2*i-1])->cancel(oldSepBel);
        if(OL >= Verbose) cout << "updated prev pos cluster" << endl;
      }

      bID = 0;

      //update first measurement clusters
      oldSepBel = sepsetBeliefs[0];
      sepsetBeliefs[0] = clustBeliefs[0]->marginalize({bID,bID+1});
      clustBeliefs[1] = clustBeliefs[1]->absorb(sepsetBeliefs[0])->cancel(oldSepBel);
      if(OL >= Verbose) cout << "updated first meas cluster" << endl;

      //update first pos clusters
      oldSepBel = sepsetBeliefs[0];
      sepsetBeliefs[0] = clustBeliefs[1]->marginalize({bID,bID+1});
      clustBeliefs[1] = clustBeliefs[1]->absorb(sepsetBeliefs[0])->cancel(oldSepBel);
      if(OL >= Verbose) cout << "updated first pos cluster" << endl;

      if(OL >= Verbose) cout << "completed back messages" << endl;
    }

    //NOTE: atm only utilising diagonal covariance elements as they aren't correlated, but this may not be the case with more complex model

    ofstream xPred("xPred.txt");
    ofstream yPred("yPred.txt");
    ofstream xCovs("xCovPred.txt");
    ofstream yCovs("yCovPred.txt");
    for(uint i = 0; i<timesteps-1;i++){

        uint bID = stateVectorSize*(i);

        rcptr<Factor> qPtr = clustBeliefs[2*i]->marginalize({bID})->normalize();
        xPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0] << ',';

        qPtr = clustBeliefs[2*i]->marginalize({bID+1})->normalize();
        yPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0] << ',';

        qPtr = clustBeliefs[2*i]->marginalize({bID,bID+1})->normalize();
        xCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[0][0] << ',';
        yCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[1][1] << ',';

    }

    bID = stateVectorSize*(timesteps-1);

    rcptr<Factor> qPtr = clustBeliefs[2*(timesteps-2)]->marginalize({bID})->normalize();
    xPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0];

    qPtr = clustBeliefs[2*(timesteps-2)]->marginalize({bID+1})->normalize();
    yPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0];

    qPtr = clustBeliefs[2*(timesteps-2)]->marginalize({bID,bID+1})->normalize();
    xCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[0][0];
    yCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[1][1];


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