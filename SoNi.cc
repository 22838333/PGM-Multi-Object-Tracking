#include "prlite_genmat.hpp"
#include "sqrtmvg.hpp"

#include "emdw.hpp"

#include "lbp_cg.hpp"
#include "lbu_cg.hpp"


#include <iostream>
#include <vector>

using namespace std;
using namespace emdw;

typedef SqrtMVG SG;

#include <fstream>



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
    cout << "end x parse" << endl;

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
          cout << "start Y parse" << endl;
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
    cout << "end x parse" << endl;

    uint timesteps = measXList.size();
    if (timesteps != measYList.size()){
      cout << "measurement dimension mismatch";
        exit(1); // terminate with error
    }

    vector<rcptr<Factor>> pdfs;
    map<RVIdType, AnyType> obsv;

    double covars = 1000.0;
    uint yIdOffset = 2*timesteps;
    for(uint i = 0; i<timesteps-1;i++){
      prlite::ColVector<double> mn1(2);
      mn1[0] = 50; mn1[1] = 50+xVelo;
      prlite::RowMatrix<double> cov1(2,2);
      cov1(0,0) = covars; cov1(0,1) = covars; cov1(1,0) = covars;  cov1(1,1) = covars+2.;
      pdfs.push_back(uniqptr<SG>(new SG({i,i+1}, mn1, cov1)));

      prlite::ColVector<double> mn2(2);
      mn2[0] = 50; mn2[1] = 50;
      prlite::RowMatrix<double> cov2(2,2);
      cov2(0,0) = covars; cov2(0,1) = covars; cov2(1,0) = covars;  cov2(1,1) = covars+1.;
      pdfs.push_back(uniqptr<SG>(new SG({i,timesteps+i}, mn2, cov2)));

      //I think this is pretty clever, but honestly I could easily be missing something stupid that's much simpler
      i += yIdOffset;

      prlite::ColVector<double> mn3(2);
      mn3[0] = 50; mn3[1] = 50+yVelo;
      prlite::RowMatrix<double> cov3(2,2);
      cov3(0,0) = covars; cov3(0,1) = covars; cov3(1,0) = covars;  cov3(1,1) = covars+2.;
      pdfs.push_back(uniqptr<SG>(new SG({i,i+1}, mn3, cov3)));

      prlite::ColVector<double> mn4(2);
      mn4[0] = 50; mn4[1] = 50;
      prlite::RowMatrix<double> cov4(2,2);
      cov4(0,0) = covars; cov4(0,1) = covars; cov4(1,0) = covars;  cov4(1,1) = covars+1.;
      pdfs.push_back(uniqptr<SG>(new SG({i,timesteps+i}, mn4, cov4)));

      i -= yIdOffset;

      obsv[timesteps+i] = measXList.at(i);
      obsv[yIdOffset+timesteps+i] = measYList.at(i);
    }

    prlite::ColVector<double> mn(2);
    mn[0] = 50; mn[1] = 50;
    prlite::RowMatrix<double> cov(2,2);
    cov(0,0) = covars; cov(0,1) = covars; cov(1,0) = covars;  cov(1,1) = covars+1.;
    pdfs.push_back(uniqptr<SG>(new SG({timesteps-1,timesteps+timesteps-1}, mn, cov)));

    prlite::ColVector<double> mn0(2);
    mn0[0] = 50; mn0[1] = 50;
    prlite::RowMatrix<double> cov0(2,2);
    cov0(0,0) = covars; cov0(0,1) = covars; cov0(1,0) = covars;  cov0(1,1) = covars+1.;
    pdfs.push_back(uniqptr<SG>(new SG({yIdOffset+timesteps-1,yIdOffset+timesteps+timesteps-1}, mn0, cov0)));

    obsv[timesteps+timesteps-1] = measXList.at(timesteps-1);
    obsv[yIdOffset+timesteps+timesteps-1] = measYList.at(timesteps-1);


    ClusterGraph cg(ClusterGraph::LTRIP, pdfs, obsv);

    // export the graph to graphviz .dot format
    cg.exportToGraphViz("SoNi2d");

    map<Idx2, rcptr<Factor> > msgs;
    MessageQueue msgQ;

    unsigned nMsgs = loopyBP_CG(cg, msgs, msgQ);
    cout << "Sent " << nMsgs << " messages before convergence\n";

    ofstream xPred("xPred.txt");
    ofstream yPred("yPred.txt");
    for(uint i = 0; i<timesteps-1;i++){

      // cout << "at pos: " << i << endl;

      // rcptr<Factor> qPtr = queryLBP_CG(cg, msgs, {i})->normalize();
      // cout << "mean: \n" << dynamic_pointer_cast<SG>(qPtr)->getMean() << endl;

      // qPtr = queryLBP_CG(cg, msgs, {i+yIdOffset})->normalize();
      // cout << "mean: \n" << dynamic_pointer_cast<SG>(qPtr)->getMean() << endl;

      rcptr<Factor> qPtr = queryLBP_CG(cg, msgs, {i})->normalize();
      xPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0] << ',';

      qPtr = queryLBP_CG(cg, msgs, {i+yIdOffset})->normalize();
      yPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0] << ',';

    }

    rcptr<Factor> qPtr = queryLBP_CG(cg, msgs, {timesteps-1})->normalize();
    xPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0];

    qPtr = queryLBP_CG(cg, msgs, {timesteps-1+yIdOffset})->normalize();
    yPred << dynamic_pointer_cast<SG>(qPtr)->getMean()[0];

    xPred.close();
    yPred.close();

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