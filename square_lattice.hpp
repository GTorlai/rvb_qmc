#ifndef SQUARELATTICE_HPP 
#define SQUARELATTICE_HPP 

#include <fstream>

class SquareLattice {

private:

  int L_;
  int Nsites_;
  int Nlinks_;

public:
        
  std::vector<std::vector<int> > Coordinates_;
  std::vector<std::vector<int> > Neighbours_;
  std::vector<std::vector<int> > LinksOnPlaquettes_;
  std::vector<std::vector<int> > SitesOnLinks_; 
  std::vector<std::vector<int> > LinksOnSites_;
  SquareLattice(int L):L_(L),
                       Nsites_(L*L),
                       Nlinks_(2*L*L){
    Init();
  
  }  
 
  inline int LinSize() const {return L_;}

  inline int Nsites() const {return Nsites_;}

  inline int Nlinks() const {return Nlinks_;}


  void Init(){
    int counter[2]={0,0};
    
    Coordinates_.resize(Nsites_,std::vector<int>(2));
    Neighbours_.resize(Nsites_,std::vector<int>(4));
    LinksOnPlaquettes_.resize(Nsites_,std::vector<int>(4));
    SitesOnLinks_.resize(Nlinks_,std::vector<int>(2));
    LinksOnSites_.resize(Nsites_,std::vector<int>(4));
  
    for (int i=1;i<Nsites_;i++){
      if ( (counter[0]+1) % L_ == 0){ //end of x-row
        //coordinate[2*i] = 0; //reset
        counter[0]=0;
          if ( (counter[1]+1) % L_ == 0)
            counter[1] = 0; //reset
          else{
            counter[1]++;
            //break;
          }
      }//if
      else {
        counter[0]++;
      }
      
      Coordinates_[i][0]=counter[0];
      Coordinates_[i][1]=counter[1];
    }
    
    //Neighbours
    for(int i=0;i<Nsites_;i++) {
      Neighbours_[i][0]=Index(Coordinates_[i][0]+1,Coordinates_[i][1]);
      Neighbours_[i][1]=Index(Coordinates_[i][0]  ,Coordinates_[i][1]+1);
      Neighbours_[i][2]=Index(Coordinates_[i][0]-1,Coordinates_[i][1]);
      Neighbours_[i][3]=Index(Coordinates_[i][0]  ,Coordinates_[i][1]-1);
    }
    
    //4 links on a plaquette
    for(int x=0;x<L_;x++) {
      for(int y=0;y<L_;y++) {
        LinksOnPlaquettes_[Index(x,y)][0]=2*Index(x,y);
        LinksOnPlaquettes_[Index(x,y)][3]=2*Index(x,y)+1;
        LinksOnPlaquettes_[Index(x,y)][1]=2*Index(x+1,y)+1;
        LinksOnPlaquettes_[Index(x,y)][2]=2*Index(x,y+1);
      }
    }
    //4 Links on a site
    for(int i=0;i<Nsites_;i++) {
        
        LinksOnSites_[i][0]=2*i;
        LinksOnSites_[i][1]=2*i+1;
        LinksOnSites_[i][2]=2*Neighbours_[i][2];
        LinksOnSites_[i][3]=2*Neighbours_[i][3]+1;
    }

    //2 sites on a link
    for(int x=0;x<L_;x++) {
      for(int y=0;y<L_;y++) {
        SitesOnLinks_[2*Index(x,y)][0] = Index(x,y);
        SitesOnLinks_[2*Index(x,y)][1] = Index(x+1,y);
        SitesOnLinks_[2*Index(x,y)+1][0] = Index(x,y);
        SitesOnLinks_[2*Index(x,y)+1][1] = Index(x,y+1);
      }
    }
  }

//Indexing of coordinates
  int Index(int x, int y) {
  
    if (x<0) x+= L_;
    if (x>=L_) x-= L_;
    if (y<0) y+= L_;
    if (y>=L_) y-= L_;
  
    return L_*y+x;
  
  }
  
//  vector<int> Decoder::generateError(MTRand & random, double p) {
//  
//      vector<int> E;
//      E.assign(Nqubits,0);
//  
//      for (int i=0; i<Nqubits; i++) {
//  
//          if (random.rand() < p) E[i] = 1;
//      }
//  
//      return E;
//  }
//  
//  
//  int Decoder::measureStar(vector<int> E, int star) {
//  
//      int stabilizer = 1;
//  
//      for (int i=0; i<4; i++) {
//  
//          stabilizer *= -(2*E[starQubits[star][i]]-1);
//      }
//  
//      return stabilizer;
//  }
//  
//  
//  vector<int> Decoder::getSyndrome(vector<int> E) {
//  
//      vector<int> S;
//      S.assign(N,1);
//  
//      int stabilizer;
//      
//      for (int i=0; i<N; i++) {
//          stabilizer = measureStar(E,i);
//          if (stabilizer == 1) S[i] = 0;
//      }
//  
//      return S;
//  }
//  
//  vector<int> Decoder::getCycle(vector<int> E0, vector<int> E) {
//  
//      vector<int> C;
//      C.assign(Nqubits,1);
//  
//      for (int i=0; i<Nqubits; i++) {
//  
//          if (E0[i] == E[i]) C[i] = 0;
//      }
//  
//      return C;
//  }
//  
//  
//  int Decoder::syndromeCheck(vector<int> E0, vector<int> E) {
//  
//      int status = 0;
//  
//      vector<int> S0 = getSyndrome(E0);
//      vector<int> S = getSyndrome(E);
//  
//      for (int i=0; i< N; i++) {
//  
//          if (S0[i] != S[i]) {
//              status = 1;
//              break;
//          }
//      }
//  
//      return status;
//  }
//  
//  int Decoder::getLogicalState(vector<int> C) {
//  
//      int status = 0;
//  
//      for (int x=0; x<L; x++) {
//          int temp = 0;
//          for (int y=0; y<L; y++) {
//  
//              temp += C[starQubits[index(x,y)][0]];
//          }
//  
//          if ((temp % 2) != 0) {
//              status = 1;
//              break;
//          }
//      }
//  
//      return status;
//  }
//  
//  int Decoder::getHomologyClass(vector<int> C) {
//  
//      int loop1=0;
//      int loop2=0;
//      int h = 0;
//  
//      for (int x=0; x<L; x++) {
//          int temp = 0;
//          for (int y=0; y<L; y++) {
//              temp += C[starQubits[index(x,y)][0]];
//          }
//          if ((temp % 2) != 0) {
//              loop1 = 1;
//          }
//      }
//      
//      for (int y=0; y<L; y++) {
//          int temp = 0;
//          for (int x=0; x<L; x++) {
//              temp += C[starQubits[index(x,y)][1]];
//          }
//          if ((temp % 2) != 0) {
//              loop2 = 1;
//          }
//      }
//  
//      if ((loop1==1) && (loop2==0)) {
//          h=1;
//      }
//      if ((loop1==0) && (loop2==1)) {
//          h=2;
//      }
//      if ((loop1==1) && (loop2==1)) {
//          h=3;
//      }
//  
//      return h;
//  }
//  
//  
//  void Decoder::testDecoder(MTRand & random) {
//  
//      vector<int> E;
//      vector<int> S;
//      vector<int> E0;
//      vector<int> S0;
//      vector<int> C;
//      
//  
//      int right_recovery = 1000;
//      int wrong_recovery = 888; 
//      
//      int right_counter=0;
//      int wrong_counter=0;
//  
//  
//      E0 = generateError(random,0.1);
//      S0 = getSyndrome(E0);
//      
//      E.assign(Nqubits,0);
//  
//      for (int i=0; i<Nqubits; i++) {
//          E[i] = E0[i];
//      }
//  
//      int plaq;
//      int loop;
//      int S_status;
//      int C_status;
//      for (int k=0; k<right_recovery; k++) {
//  
//          plaq = random.randInt(N-1);
//  
//          for (int i=0; i<4; i++) {
//              E[plaqQubits[plaq][i]] ^=1;
//          }
//          
//          C = getCycle(E0,E);
//          C_status = getLogicalState(C);
//          if (C_status ==0) right_counter++;
//      }
//      
//      for (int k=0; k<wrong_recovery; k++) {
//  
//          plaq = random.randInt(N-1);
//          loop = random.randInt(L-1);
//  
//          for (int i=0; i<L; i++) {
//              E[2*(L*loop + i)] ^= 1;
//          }
//          
//          C = getCycle(E0,E);
//          C_status = getLogicalState(C);
//          
//          if (C_status ==1) wrong_counter++;
//   
//          for (int i=0; i<L; i++) {
//              E[2*(L*loop + i)] ^= 1;
//          }
//  
//          for (int i=0; i<4; i++) {
//              E[plaqQubits[plaq][i]] ^=1;
//          }
//          
//      }
//  
//      cout << "Number of right states: " << right_counter << endl;
//      cout << "Number of wrong states: " << wrong_counter << endl;
//   
//      //S_status = syndromeCheck(E0,E);
//      //
//      //cout << "Syndrome ";
//      //if (S_status == 0) cout << " CORRECT" << endl << endl;
//      //else cout << " INCORRECT" << endl << endl;
//  
//      
//      
//      //cout << "Logical State: ";
//      //if (C_status == 0) cout << " PROTECTED" << endl << endl;
//      //else cout << "CORRUPTED" << endl << endl;
//      
//  
//  }
//  
//  int Decoder::qIndex(int j, int s) {
//  
//      int index;
//  
//      for (int i=0; i<4; i++) {
//  
//          if (starQubits[s][i] == j) {
//              index = i;
//          }
//      }
//      return index;
//  }
//  
  
  //Print Lattice Informations
  void PrintLattice() {
  
      std::cout << std::endl << std::endl << "Printing Indexing of Coordinates..." << std::endl << std::endl;
  
      for (int x=0; x<L_; x++) {
          for (int y=0; y<L_; y++) {
              std::cout << "Coordinate (" << x << "," << y << ")";
              std::cout << "  -->  Index: " << Index(x,y) << std::endl;
          }
      }
  
      std::cout << std::endl << std::endl << std::endl;
  
      std::cout << "Printing coordinates of indices..." << std::endl << std::endl;
      
      for (int i=0; i<Nsites_; i++) {
  
          std::cout << "Index "<<i;
          std::cout << "  --> Coordinates: (";
          std::cout << Coordinates_[i][0] <<","<<Coordinates_[i][1]<<")";
          std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;
  
      std::cout << "Printing neighbors..." << std::endl << std::endl;
      
      for (int i=0; i<Nsites_; i++) {
  
          std::cout << "Index "<<i;
          std::cout << "  --> Neighbors: ";
          
          for (int j=0; j<4;j++) {
              std::cout << Neighbours_[i][j] << " , ";
          }
          
          std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;
      
      std::cout << "Printing links on plaquettes..." << std::endl << std::endl;
   
      for (int i=0; i<Nsites_; i++) {
  
          std::cout << "Plaquette "<<i;
          std::cout << "  --> Links: ";
          
          for (int j=0; j<4;j++) {
              std::cout << LinksOnPlaquettes_[i][j] << " , ";
          }
          
          std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;
      
      std::cout << "Printing sites on links..." << std::endl << std::endl;
   
      for (int i=0; i<Nlinks_; i++) {
  
          std::cout << "Link "<<i;
          std::cout << "  --> Sites: ";
          
          for (int j=0; j<2;j++) {
              std::cout << SitesOnLinks_[i][j] << " , ";
          }
          
          std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;

      std::cout << "Printing links on sites..." << std::endl << std::endl;
   
      for (int i=0; i<Nsites_; i++) {
  
          std::cout << "Site "<<i;
          std::cout << "  --> Links: ";
          
          for (int j=0; j<4;j++) {
              std::cout << LinksOnSites_[i][j] << " , ";
          }
          
          std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;
      
      //cout << "Printing stars on links..." << endl << endl;
   
      //for (int i=0; i<Nqubits; i++) {
  
      //    cout << "Link "<<i;
      //    cout << "  --> Stars: ";
      //    
      //    for (int j=0; j<2;j++) {
      //        std::cout << qubitStars[i][j] << " , ";
      //    }
      //    
      //    std::cout << std::endl;
      //}
      //std::cout << std::endl << std::endl << std::endl;
  
  }
  
  
//  //Print Vector on terminal
//  void Decoder::printVector(vector<int> Vector){
//  
//      for (int i=0;i<Vector.size();i++){
//          cout<<Vector[i]<<" ";
//      }//i
//      cout<<endl;
//  
//  }//print
//  
//  //Write Vector on file
//  void Decoder::writeVector(vector<int> Vector, ofstream & file){
//  
//     for (int i=0;i<Vector.size();i++){
//          file<<Vector[i]<<" ";
//      }//i
//      file<<endl; 
//  }
//

};

#endif
