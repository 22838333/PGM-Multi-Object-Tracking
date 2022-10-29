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

    double xVelo = 0;
    double yVelo = 0;

    inFile.open("dir.txt");
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }

    while (inFile >> c) {
      switch(c) {
        case '[':
          cout << "start dir parse" << endl;
          break;
        case ',':
          xVelo = std::stod(str);
          str.clear();
          break;
        case ']':
          yVelo = std::stod(str);
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
    cout << "end dir parse" << endl;

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

    //construct the cluster graph
    for(uint i = 0; i<timesteps-1;i++){

      uint bID = 2*i;

      //constructs an affine gaussian for the point to point relations
      prlite::ColMatrix<double> A(2,2);
      prlite::ColVector<double> C(2);
      prlite::ColMatrix<double> R(2,2);
      A(0,0)=1; A(0,1)=0; A(1,0)=0; A(1,1)=1;
      C[0]=xVelo; C[1]=yVelo;
      R(0,0)=1; R(0,1)=0; R(1,0)=0; R(1,1)=1;
      posPdfs.push_back(uniqptr<AG>(new AG({bID,bID+1},{bID+2,bID+3}, A, C, R)));

      //constructs a mvg based around the measurements
      prlite::ColVector<double> mn(2);
      prlite::RowMatrix<double> cov(2,2);
      mn[0] = measXList[i]; mn[1] = measYList[i];
      cov(0,0) = 1.0; cov(1,0) = 0.0; cov(0,1) = 0.0; cov(1,1) = 1.0;
      measPdfs.push_back(uniqptr<SG>(new SG({bID,bID+1}, mn, cov)));
    }

    //constructs the final measurements mvg as there is one more than the number of clusters
    uint bID = 2*(timesteps-1);
    prlite::ColVector<double> mn(2);
    prlite::RowMatrix<double> cov(2,2);
    mn[0] = measXList[timesteps-1]; mn[1] = measYList[timesteps-1];
    cov(0,0) = 1.0; cov(1,0) = 0.0; cov(0,1) = 0.0; cov(1,1) = 1.0;
    measPdfs.push_back(uniqptr<SG>(new SG({bID,bID+1}, mn, cov)));


    // vector<rcptr<Factor>> clustBeliefs(timesteps-1);
    // vector<rcptr<Factor>> forwardMessages(timesteps-2);
    // vector<rcptr<Factor>> backMessages(timesteps-2);

    // // Send messages forward
    // for(uint i = 0; i<timesteps-2;i++){
    //   uint bID = 2*i;
    //   if(i == 0){
    //     forwardMessages[i] = posPdfs[i]->absorb(measPdfs[i])->marginalize({bID+2,bID+3});
    //   }else{
    //     forwardMessages[i] = posPdfs[i]->absorb(measPdfs[i])->absorb(forwardMessages[i-1])->marginalize({bID+2,bID+3});
    //   }
    // }

    // //send messages backwards
    // for(int i = timesteps-3; i>=0 ;i--){
    //   uint bID = 2*i+2;
    //   if((uint)i == timesteps-3){
    //     backMessages[i] = posPdfs[i+1]->absorb(measPdfs[i+1])->absorb(measPdfs[i+2])->marginalize({bID,bID+1});
    //   }else{
    //     backMessages[i] = posPdfs[i+1]->absorb(measPdfs[i+1])->absorb(backMessages[i+1])->marginalize({bID,bID+1});
    //   }
    // }

    // //construct cluster beliefs
    // for(uint i = 0; i<timesteps-1;i++){
    //   if(i == 0){
    //     clustBeliefs[i] = posPdfs[i]->absorb(measPdfs[i])->absorb(backMessages[i]);
    //   }else if(i == timesteps-2){
    //     clustBeliefs[i] = posPdfs[i]->absorb(measPdfs[i])->absorb(measPdfs[i+1])->absorb(forwardMessages[i-1]);
    //   }else{
    //     clustBeliefs[i] = posPdfs[i]->absorb(measPdfs[i])->absorb(forwardMessages[i-1])->absorb(backMessages[i]);
    //   }
    // }

    // //NOTE: atm only utilising diagonal covariance elements as aren't corrolated, but this will not be the case with more complex model

    // ofstream xPred("xPred.txt");
    // ofstream yPred("yPred.txt");
    // ofstream xCovs("xCovPred.txt");
    // ofstream yCovs("yCovPred.txt");
    // for(uint i = 0; i<timesteps-1;i++){

    //   rcptr<Factor> qPtr = clustBeliefs[i]->marginalize({2*i})->normalize();
    //   xPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0] << ',';

    //   qPtr = clustBeliefs[i]->marginalize({2*i+1})->normalize();
    //   yPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0] << ',';

    //   qPtr = clustBeliefs[i]->marginalize({2*i,2*i+1})->normalize();
    //   xCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[0][0] << ',';
    //   yCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[1][1] << ',';

    // }

    // rcptr<Factor> qPtr = clustBeliefs[timesteps-2]->marginalize({2*(timesteps-1)})->normalize();
    // xPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0];

    // qPtr = clustBeliefs[timesteps-2]->marginalize({2*(timesteps-1)+1})->normalize();
    // yPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0];

    // qPtr = clustBeliefs[timesteps-2]->marginalize({2*(timesteps-1),2*(timesteps-1)+1})->normalize();
    // xCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[0][0];
    // yCovs << dynamic_pointer_cast<SG>(qPtr)->getCov()[1][1];


    // xPred.close();
    // yPred.close();
    // xCovs.close();
    // yCovs.close();

    // //instantiate base cluster belief
    // for(uint i = 0; i<timesteps-1;i++){
    //   clustBeliefs.push_back(posPdfs[i]);
    // }

    // //Include info from measurements
    // for(uint i = 0; i<timesteps-1;i++){
    //   clustBeliefs[i] = clustBeliefs[i]->absorb(measPdfs[i]);
    // }
    // clustBeliefs[timesteps-2] = clustBeliefs[timesteps-2]->absorb(measPdfs[timesteps-1]);

    // //Send messages forward
    // for(uint i = 0; i<timesteps-2;i++){
    //   uint bID = 2*i;
    //   forwardMessages.push_back(clustBeliefs[i]->marginalize({bID+2,bID+3}));
    //   clustBeliefs[i+1] = clustBeliefs[i+1]->absorb(forwardMessages[i]);
    // }

    // //send messages back to prev clust
    // for(uint i = 0; i<timesteps-2;i++){
    //   uint invNum = timesteps-2-i;
    //   uint bID = 2*(invNum);
    //   backMessages.push_back(clustBeliefs[invNum]->cancel(forwardMessages[invNum-1])->marginalize({bID,bID+1}));
    //   clustBeliefs[invNum-1] = clustBeliefs[invNum-1]->absorb(backMessages[i]);
    // }

    vector<rcptr<Factor>> clustBeliefs2;
    vector<rcptr<Factor>> forwardMessages2;
    vector<rcptr<Factor>> backMessages2;

    //instantiate base cluster belief
    for(uint i = 0; i<timesteps-1;i++){
      clustBeliefs2.push_back(posPdfs[i]);
    }

    //Include info from measurements
    for(uint i = 0; i<timesteps-1;i++){
      clustBeliefs2[i] = clustBeliefs2[i]->absorb(measPdfs[i]);
    }
    clustBeliefs2[timesteps-2] = clustBeliefs2[timesteps-2]->absorb(measPdfs[timesteps-1]);

    //Send messages forward
    for(uint i = 0; i<timesteps-2;i++){
      uint bID = 2*i;
      forwardMessages2.push_back(clustBeliefs2[i]->marginalize({bID+2,bID+3}));
      clustBeliefs2[i+1] = clustBeliefs2[i+1]->absorb(forwardMessages2[i]);
    }

    //send messages back to prev clust
    for(uint i = 0; i<timesteps-2;i++){
      uint invNum = timesteps-2-i;
      uint bID = 2*(invNum);
      // put this in to see how the absorb + cancel doesn't work on the affinemvgs
      if(i == 0){
        rcptr<Factor> ptr1 = (clustBeliefs2[invNum]->cancel(forwardMessages2[invNum-1]));
        rcptr<Factor> ptr2 = (posPdfs[invNum]->absorb(measPdfs[invNum])->absorb(measPdfs[invNum+1]));
        rcptr<Factor> ptr5 = (ptr2->absorb(forwardMessages2[invNum-1])->cancel(forwardMessages2[invNum-1]));
        rcptr<Factor> ptr3 = ptr1->marginalize({bID,bID+1,bID+2,bID+3});
        rcptr<Factor> ptr4 = ptr2->marginalize({bID,bID+1,bID+2,bID+3});
        rcptr<Factor> ptr6 = ptr5->marginalize({bID,bID+1,bID+2,bID+3});
        cout << "0 unmarg back Message: " << i << "\n" << *ptr1 << endl;
        cout << "0 abs+cans Message: " << i << "\n" << *ptr5 << endl;
        cout << "0 constructed back Message: " << i << "\n" << *ptr2 << endl;
        cout << "0 marg back Message: " << i << "\n" << *ptr3 << endl;
        cout << "0 marg abs+cans Message: " << i << "\n" << *ptr6 << endl;
        cout << "0 marg constructed back Message: " << i << "\n" << *ptr4 << endl;
      }
      backMessages2.push_back(clustBeliefs2[invNum]->cancel(forwardMessages2[invNum-1])->marginalize({bID,bID+1}));
      clustBeliefs2[invNum-1] = clustBeliefs2[invNum-1]->absorb(backMessages2[i]);
    }

    ofstream xPred2("xPred2.txt");
    ofstream yPred2("yPred2.txt");
    for(uint i = 0; i<timesteps-1;i++){

      rcptr<Factor> qPtr2 = clustBeliefs2[i]->marginalize({2*i})->normalize();
      xPred2 << dynamic_pointer_cast<SG>(qPtr2)->getMean()[0] << ',';

      qPtr2 = clustBeliefs2[i]->marginalize({2*i+1})->normalize();
      yPred2 << dynamic_pointer_cast<SG>(qPtr2)->getMean()[0] << ',';

    }

    rcptr<Factor> qPtr2 = clustBeliefs2[timesteps-2]->marginalize({2*(timesteps-1)})->normalize();
    xPred2 << dynamic_pointer_cast<SG>(qPtr2)->getMean()[0];

    qPtr2 = clustBeliefs2[timesteps-2]->marginalize({2*(timesteps-1)+1})->normalize();
    yPred2 << dynamic_pointer_cast<SG>(qPtr2)->getMean()[0];

    xPred2.close();
    yPred2.close();

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