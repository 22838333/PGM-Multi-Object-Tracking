#include "prlite_genmat.hpp"
#include "emdw.hpp"
// #include "sqrtmvg.hpp"
#include "affinemvg.hpp"

#include "lbp_cg.hpp"
#include "lbu_cg.hpp"


#include <iostream>
#include <vector>
#include <fstream>

using namespace std;
using namespace emdw;

typedef SqrtMVG SG;
typedef AffineMVG AG;



int main(int, char *argv[]) {
  try {

    vector<double> measXList;
    vector<double> measYList;

    char c;
    std::string str;

    ifstream inFile;

    inFile.open("Xpos.txt");
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }

    while (inFile >> c) {
      switch(c) {
        case '[':
          cout << "start x parse" << endl;
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
    cout << "end x parse" << endl;

    inFile.open("Ypos.txt");
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }

    while (inFile >> c) {
      switch(c) {
        case '[':
          cout << "start Y parse" << endl;
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
    cout << "end y parse" << endl;

    uint timesteps = measXList.size();
    if (timesteps < 2){
      cout << "At least 2 timesteps required";
        exit(1); // terminate with error
    }
    if (timesteps != measYList.size()){
      cout << "measurement dimension mismatch";
        exit(1); // terminate with error
    }

    vector<rcptr<Factor>> posPdfs;
    vector<rcptr<Factor>> measPdfs;
    map<RVIdType, AnyType> obsv;

    uint stateVectorSize = 4;

    //construct the cluster graph
    for(uint i = 0; i<timesteps-1;i++){

        uint bID = stateVectorSize*i;

        //constructs an affine gaussian for the point to point relations that includes a unknown velocity
        prlite::ColMatrix<double> A(4,4);
        prlite::ColVector<double> C(4);
        prlite::ColMatrix<double> R(4,4);
        A(0,0)=1; A(0,1)=0; A(0,2)=1; A(0,3)=0;
        A(1,0)=0; A(1,1)=1; A(1,2)=0; A(1,3)=1;
        A(2,0)=0; A(2,1)=0; A(2,2)=1; A(2,3)=0;
        A(3,0)=0; A(3,1)=0; A(3,2)=0; A(3,3)=1;
        C[0]=0; C[1]=0; C[2]=0; C[3]=0;
        R(0,0)=1; R(0,1)=0; R(0,2)=0; R(0,3)=0;
        R(1,0)=0; R(1,1)=1; R(1,2)=0; R(1,3)=0;
        R(2,0)=0; R(2,1)=0; R(2,2)=0.01; R(2,3)=0;
        R(3,0)=0; R(3,1)=0; R(3,2)=0; R(3,3)=0.01;
        posPdfs.push_back(uniqptr<AG>(new AG({bID,bID+1,bID+2,bID+3},{bID+4,bID+5,bID+6,bID+7}, A, C, R)));

        prlite::ColVector<double> mn(4);
        prlite::RowMatrix<double> cov(4,4);
        mn[0] = measXList[i]; mn[1] = measYList[i]; mn[2] = 0; mn[3] = 0;
        cov(0,0) = 1.0; cov(0,1) = 0.0; cov(0,2) = 0.0; cov(0,3) = 0.0;
        cov(1,0) = 0.0; cov(1,1) = 1.0; cov(1,2) = 0.0; cov(1,3) = 0.0;
        cov(2,0) = 0.0; cov(2,1) = 0.0; cov(2,2) = 10000.0; cov(2,3) = 0.0;
        cov(3,0) = 0.0; cov(3,1) = 0.0; cov(3,2) = 0.0; cov(3,3) = 10000.0;
        measPdfs.push_back(uniqptr<SG>(new SG({bID,bID+1,bID+2,bID+3}, mn, cov)));
    }

    //constructs the final measurements mvg as there is one more than the number of clusters
    uint bID = stateVectorSize*(timesteps-1);
    prlite::ColVector<double> mn(4);
    prlite::RowMatrix<double> cov(4,4);
    mn[0] = measXList[timesteps-1]; mn[1] = measYList[timesteps-1]; mn[2] = 0; mn[3] = 0;
    cov(0,0) = 1.0; cov(0,1) = 0.0; cov(0,2) = 0.0; cov(0,3) = 0.0;
    cov(1,0) = 0.0; cov(1,1) = 1.0; cov(1,2) = 0.0; cov(1,3) = 0.0;
    cov(2,0) = 0.0; cov(2,1) = 0.0; cov(2,2) = 10000.0; cov(2,3) = 0.0;
    cov(3,0) = 0.0; cov(3,1) = 0.0; cov(3,2) = 0.0; cov(3,3) = 10000.0;
    measPdfs.push_back(uniqptr<SG>(new SG({bID,bID+1,bID+2,bID+3}, mn, cov)));


    vector<rcptr<Factor>> clustBeliefs(timesteps-1);
    vector<rcptr<Factor>> forwardMessages(timesteps-2);
    vector<rcptr<Factor>> backMessages(timesteps-2);

    // Send messages forward
    for(uint i = 0; i<timesteps-2;i++){
      uint bID = stateVectorSize*i;
      if(i == 0){
        forwardMessages[i] = posPdfs[i]->absorb(measPdfs[i])->marginalize({bID+4,bID+5,bID+6,bID+7});
      }else{
        forwardMessages[i] = posPdfs[i]->absorb(measPdfs[i])->absorb(forwardMessages[i-1])->marginalize({bID+4,bID+5,bID+6,bID+7});
      }
    }

    //send messages backwards
    for(int i = timesteps-3; i>=0 ;i--){
      uint bID = stateVectorSize*(i+1);
      if((uint)i == timesteps-3){
        backMessages[i] = posPdfs[i+1]->absorb(measPdfs[i+1])->absorb(measPdfs[i+2])->marginalize({bID,bID+1,bID+2,bID+3});
      }else{

        backMessages[i] = posPdfs[i+1]->absorb(measPdfs[i+1])->absorb(backMessages[i+1])->marginalize({bID,bID+1,bID+2,bID+3});
      }
    }

    //construct cluster beliefs
    for(uint i = 0; i<timesteps-1;i++){
      if(i == 0){
        clustBeliefs[i] = posPdfs[i]->absorb(measPdfs[i])->absorb(backMessages[i]);
      }else if(i == timesteps-2){
        clustBeliefs[i] = posPdfs[i]->absorb(measPdfs[i])->absorb(measPdfs[i+1])->absorb(forwardMessages[i-1]);
      }else{
        clustBeliefs[i] = posPdfs[i]->absorb(measPdfs[i])->absorb(forwardMessages[i-1])->absorb(backMessages[i]);
      }
    }

    //NOTE: atm only utilising diagonal covariance elements as they aren't corrolated, but this will likely not be the case with more complex model

    ofstream xPred("xPred.txt");
    ofstream yPred("yPred.txt");
    ofstream xCovs("xCovPred.txt");
    ofstream yCovs("yCovPred.txt");
    for(uint i = 0; i<timesteps-1;i++){

        uint bID = stateVectorSize*(i);

        rcptr<Factor> qPtr = clustBeliefs[i]->marginalize({bID})->normalize();
        xPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0] << ',';

        qPtr = clustBeliefs[i]->marginalize({bID+1})->normalize();
        yPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0] << ',';

        qPtr = clustBeliefs[i]->marginalize({bID,bID+1})->normalize();
        xCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[0][0] << ',';
        yCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[1][1] << ',';

    }

    bID = stateVectorSize*(timesteps-1);

    rcptr<Factor> qPtr = clustBeliefs[timesteps-2]->marginalize({bID})->normalize();
    xPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0];

    qPtr = clustBeliefs[timesteps-2]->marginalize({bID+1})->normalize();
    yPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0];

    qPtr = clustBeliefs[timesteps-2]->marginalize({bID,bID+1})->normalize();
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