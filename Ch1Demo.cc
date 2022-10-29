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

    uint stateVectorSize = 6;

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

        if(argv[1][0] == 'v'){
            A(2,4)=0;A(3,5)=0;A(0,4)=0;A(1,5)=0;
            // cout << "v" << endl;
        }
        if(argv[1][0] == 'p'){
            A(2,4)=0;A(3,5)=0;A(0,4)=0;A(1,5)=0;
            A(0,2)=0;A(1,3)=0;
            // cout << "p\n" << A << endl;
        }


        C[0]=0; C[1]=0; C[2]=0; C[3]=0; C[4]=0; C[5]=0;

        R(0,0)=1; R(0,1)=0; R(0,2)=0; R(0,3)=0; R(0,4)=0; R(0,5)=0;
        R(1,0)=0; R(1,1)=1; R(1,2)=0; R(1,3)=0; R(1,4)=0; R(1,5)=0;
        R(2,0)=0; R(2,1)=0; R(2,2)=0.01; R(2,3)=0; R(2,4)=0; R(2,5)=0;
        R(3,0)=0; R(3,1)=0; R(3,2)=0; R(3,3)=0.01; R(3,4)=0; R(3,5)=0;
        R(4,0)=0; R(4,1)=0; R(4,2)=0; R(4,3)=0; R(4,4)=0.01; R(4,5)=0;
        R(5,0)=0; R(5,1)=0; R(5,2)=0; R(5,3)=0; R(5,4)=0; R(5,5)=0.01;

        posPdfs.push_back(uniqptr<AG>(new AG({bID,bID+1,bID+2,bID+3,bID+4,bID+5},{bID+6,bID+7,bID+8,bID+9,bID+10,bID+11}, A, C, R)));

        prlite::ColVector<double> mn(6);
        prlite::RowMatrix<double> cov(6,6);
        mn[0] = measXList[i]; mn[1] = measYList[i]; mn[2] = 0; mn[3] = 0; mn[4] = 0; mn[5] = 0;
        cov(0,0) = 3.0; cov(0,1) = 0.0; cov(0,2) = 0.0; cov(0,3) = 0.0; cov(0,4) = 0.0; cov(0,5) = 0.0;
        cov(1,0) = 0.0; cov(1,1) = 3.0; cov(1,2) = 0.0; cov(1,3) = 0.0; cov(1,4) = 0.0; cov(1,5) = 0.0;
        cov(2,0) = 0.0; cov(2,1) = 0.0; cov(2,2) = 10000.0; cov(2,3) = 0.0; cov(2,4) = 0.0; cov(2,5) = 0.0;
        cov(3,0) = 0.0; cov(3,1) = 0.0; cov(3,2) = 0.0; cov(3,3) = 10000.0; cov(3,4) = 0.0; cov(3,5) = 0.0;
        cov(4,0) = 0.0; cov(4,1) = 0.0; cov(4,2) = 0.0; cov(4,3) = 0.0; cov(4,4) = 10000.0; cov(4,5) = 0.0;
        cov(5,0) = 0.0; cov(5,1) = 0.0; cov(5,2) = 0.0; cov(5,3) = 0.0; cov(5,4) = 0.0; cov(5,5) = 10000.0;
        measPdfs.push_back(uniqptr<SG>(new SG({bID,bID+1,bID+2,bID+3,bID+4,bID+5}, mn, cov)));
    }

    //constructs the final measurements mvg as there is one more than the number of clusters
    uint bID = stateVectorSize*(timesteps-1);
    prlite::ColVector<double> mn(6);
    prlite::RowMatrix<double> cov(6,6);
    mn[0] = measXList[timesteps-1]; mn[1] = measYList[timesteps-1]; mn[2] = 0; mn[3] = 0; mn[4] = 0; mn[5] = 0;
    cov(0,0) = 3.0; cov(0,1) = 0.0; cov(0,2) = 0.0; cov(0,3) = 0.0; cov(0,4) = 0.0; cov(0,5) = 0.0;
    cov(1,0) = 0.0; cov(1,1) = 3.0; cov(1,2) = 0.0; cov(1,3) = 0.0; cov(1,4) = 0.0; cov(1,5) = 0.0;
    cov(2,0) = 0.0; cov(2,1) = 0.0; cov(2,2) = 10000.0; cov(2,3) = 0.0; cov(2,4) = 0.0; cov(2,5) = 0.0;
    cov(3,0) = 0.0; cov(3,1) = 0.0; cov(3,2) = 0.0; cov(3,3) = 10000.0; cov(3,4) = 0.0; cov(3,5) = 0.0;
    cov(4,0) = 0.0; cov(4,1) = 0.0; cov(4,2) = 0.0; cov(4,3) = 0.0; cov(4,4) = 1000.0; cov(4,5) = 0.0;
    cov(5,0) = 0.0; cov(5,1) = 0.0; cov(5,2) = 0.0; cov(5,3) = 0.0; cov(5,4) = 0.0; cov(5,5) = 1000.0;
    measPdfs.push_back(uniqptr<SG>(new SG({bID,bID+1,bID+2,bID+3,bID+4,bID+5}, mn, cov)));

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
      clustBeliefs[i*2+1] = measPdfs[i];
    }
    clustBeliefs[(timesteps-2)*2+2] = measPdfs[timesteps-1];

    cout << "completed construction " << endl;

    //send messages forwards
    clustBeliefs[0] = clustBeliefs[0]->absorb(expect)->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5,bID+6,bID+7,bID+8,bID+9,bID+10,bID+11});
    for(uint i = 0; i<timesteps-2;i++){
      uint bID = stateVectorSize*(i);

      // cout << "loop " << i << "\n" << endl;

      //update measurement clusters
      sepsetBeliefs[2*i] = clustBeliefs[2*i]->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5});
      clustBeliefs[2*i+1] = clustBeliefs[2*i+1]->absorb(sepsetBeliefs[2*i]);

      //update current pos clusters
      oldSepBel = sepsetBeliefs[2*i];
      sepsetBeliefs[2*i] = clustBeliefs[2*i+1]->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5});
      clustBeliefs[2*i] = clustBeliefs[2*i]->absorb(sepsetBeliefs[2*i])->cancel(oldSepBel);

      //update next pos clusters
      sepsetBeliefs[2*i+1] = clustBeliefs[2*i]->marginalize({bID+6,bID+7,bID+8,bID+9,bID+10,bID+11});
      clustBeliefs[2*i+2] = clustBeliefs[2*i+2]->absorb(sepsetBeliefs[2*i+1])->marginalize({bID+6,bID+7,bID+8,bID+9,bID+10,bID+11,bID+12,bID+13,bID+14,bID+15,bID+16,bID+17});
    }

    bID = stateVectorSize*(timesteps-2);

    //update second last measurement cluster
    sepsetBeliefs[2*(timesteps-2)] = clustBeliefs[2*(timesteps-2)]->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5});
    clustBeliefs[2*(timesteps-2)+1] = clustBeliefs[2*(timesteps-2)+1]->absorb(sepsetBeliefs[2*(timesteps-2)]);

    //update last pos clusters
    oldSepBel = sepsetBeliefs[2*(timesteps-2)];
    sepsetBeliefs[2*(timesteps-2)] = clustBeliefs[2*(timesteps-2)+1]->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5});
    clustBeliefs[2*(timesteps-2)] = clustBeliefs[2*(timesteps-2)]->absorb(sepsetBeliefs[2*(timesteps-2)])->cancel(oldSepBel);

    //update last measurement cluster
    sepsetBeliefs[2*(timesteps-2)+1] = clustBeliefs[2*(timesteps-2)]->marginalize({bID+6,bID+7,bID+8,bID+9,bID+10,bID+11});
    clustBeliefs[2*(timesteps-2)+2] = clustBeliefs[2*(timesteps-2)+2]->absorb(sepsetBeliefs[2*(timesteps-2)+1]);

    //update last pos clusters
    oldSepBel = sepsetBeliefs[2*(timesteps-2)+1];
    sepsetBeliefs[2*(timesteps-2)+1] = clustBeliefs[2*(timesteps-2)+2]->marginalize({bID+6,bID+7,bID+8,bID+9,bID+10,bID+11});
    clustBeliefs[2*(timesteps-2)] = clustBeliefs[2*(timesteps-2)]->absorb(sepsetBeliefs[2*(timesteps-2)+1])->cancel(oldSepBel);

    cout << "completed forward messages" << endl;

    //send messages backwards
    for(int i = timesteps-2; i>0 ;i--){ //-1 from 0 count, -1 from clust amounts, -1 from sep amount
      uint bID = stateVectorSize*(i);

      //update measurement clusters
      oldSepBel = sepsetBeliefs[2*i];
      sepsetBeliefs[2*i] = clustBeliefs[2*i]->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5});
      clustBeliefs[2*i+1] = clustBeliefs[2*i+1]->absorb(sepsetBeliefs[2*i])->cancel(oldSepBel);

      // //update current pos clusters (uneeded as current since there is no new info in the measuremnt cluster, however will be when iterative)
      // oldSepBel = sepsetBeliefs[2*i];
      // sepsetBeliefs[2*i] = clustBeliefs[2*i+1]->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5});
      // clustBeliefs[2*i+1] = clustBeliefs[2*i+1]->absorb(sepsetBeliefs[2*i])->cancel(oldSepBel);

      //update prev pos clusters
      oldSepBel = sepsetBeliefs[2*i-1];
      sepsetBeliefs[2*i-1] = clustBeliefs[2*i]->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5});
      clustBeliefs[2*i-2] = clustBeliefs[2*i-2]->absorb(sepsetBeliefs[2*i-1])->cancel(oldSepBel);

    }

    bID = 0;

    //update first measurement clusters
    oldSepBel = sepsetBeliefs[0];
    sepsetBeliefs[0] = clustBeliefs[0]->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5});
    clustBeliefs[1] = clustBeliefs[1]->absorb(sepsetBeliefs[0])->cancel(oldSepBel);

    // //update current pos clusters (uneeded as current since there is no new info in the measuremnt cluster, however will be when iterative)
    // oldSepBel = sepsetBeliefs[0];
    // sepsetBeliefs[0] = clustBeliefs[1]->marginalize({bID,bID+1,bID+2,bID+3,bID+4,bID+5});
    // clustBeliefs[1] = clustBeliefs[1]->absorb(sepsetBeliefs[0])->cancel(oldSepBel);

    cout << "completed back messages" << endl;

    //NOTE: atm only utilising diagonal covariance elements as they aren't corrolated, but this will likely not be the case with more complex model

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