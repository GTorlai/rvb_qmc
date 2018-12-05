#ifndef CUBICLATTICE_HPP 
#define CUBICLATTICE_HPP 

#include <fstream>

class CubicLattice {

private:

  int L_;
  int Nsites_;
  int Nbonds_;
  int Nplaqs_;

public:
        
  std::vector<std::vector<int> > Coordinates_;        // (x,y) coordinates of each lattice site
  std::vector<std::vector<int> > neighbours_;         // 4 neighbours of each lattice site
  std::vector<std::vector<int> > SitesOnPlaquettes_;  // 4 sites on each plaquette
  
  CubicLattice(int L){
    
    L_ = L;
    Nsites_ = 2*L*L*L;  // The factor 2 is for the two replicas
    Nbonds_ = 3*Nsites_;
    Nplaqs_ = 2 * 2 * L*(L*L);
    std::cout<<Nplaqs_<<std::endl;
    Init();
  
  }  
 
  inline int LinSize() const {return L_;}

  inline int Nsites() const {return Nsites_;}

  inline int Nbonds() const {return Nbonds_;}

  inline int Nplaqs() const {return Nplaqs_;}

  void Init(){
    int counter[3]={0,0,0};
    
    Coordinates_.resize(Nsites_,std::vector<int>(3));
    neighbours_.resize(Nsites_,std::vector<int>(6));
    SitesOnPlaquettes_.resize(Nplaqs_,std::vector<int>(4));

    //Coordinates
    for (int i=1;i<Nsites_/2;i++){
      if ( (counter[0]+1) % L_ == 0){//end of x-row
        //coordinate[2*i] = 0; //reset
        counter[0]=0;
        if ( (counter[1]+1) % L_ == 0){
          counter[1] = 0; //reset
          if ( (counter[2]+1) % L_ == 0)
            counter[2] = 0;
          else{
            counter[2]++;
          }
        }
        else{
          counter[1]++;
        }
      }//if
      else {
        counter[0]++;
      }
      
      Coordinates_[i][0]=counter[0];
      Coordinates_[i][1]=counter[1];
      Coordinates_[i][2]=counter[2];
    }
    
    //Neighbours
    // Real
    for(int i=0;i<Nsites_/2;i++) {
      neighbours_[i][0]=Index(Coordinates_[i][0]+1,Coordinates_[i][1]  ,Coordinates_[i][2]);
      neighbours_[i][1]=Index(Coordinates_[i][0]  ,Coordinates_[i][1]+1,Coordinates_[i][2]);
      neighbours_[i][2]=Index(Coordinates_[i][0]  ,Coordinates_[i][1]  ,Coordinates_[i][2]+1);
      neighbours_[i][3]=Index(Coordinates_[i][0]-1,Coordinates_[i][1]  ,Coordinates_[i][2]);
      neighbours_[i][4]=Index(Coordinates_[i][0]  ,Coordinates_[i][1]-1,Coordinates_[i][2]);
      neighbours_[i][5]=Index(Coordinates_[i][0]  ,Coordinates_[i][1]  ,Coordinates_[i][2]-1);
      
    }
    // Replica
    for(int i=0;i<Nsites_/2;i++) {
      neighbours_[Nsites_/2+i][0]=Nsites_/2 + Index(Coordinates_[i][0]+1,Coordinates_[i][1]  ,Coordinates_[i][2]);
      neighbours_[Nsites_/2+i][1]=Nsites_/2 + Index(Coordinates_[i][0]  ,Coordinates_[i][1]+1,Coordinates_[i][2]);
      neighbours_[Nsites_/2+i][2]=Nsites_/2 + Index(Coordinates_[i][0]  ,Coordinates_[i][1]  ,Coordinates_[i][2]+1);
      neighbours_[Nsites_/2+i][3]=Nsites_/2 + Index(Coordinates_[i][0]-1,Coordinates_[i][1]  ,Coordinates_[i][2]);
      neighbours_[Nsites_/2+i][4]=Nsites_/2 + Index(Coordinates_[i][0]  ,Coordinates_[i][1]-1,Coordinates_[i][2]);
      neighbours_[Nsites_/2+i][5]=Nsites_/2 + Index(Coordinates_[i][0]  ,Coordinates_[i][1]  ,Coordinates_[i][2]-1);
      
    }

    //Sites on Plaquettes
    //Real
    int cnt=0;
    for(int z=0;z<L_;z++){
      //Plaquettes lying on the plane z (z is constant!)
      for(int y=0;y<L_;y++) {
        for(int x=0;x<L_;x++) {
          SitesOnPlaquettes_[cnt][0] = Index(x,y,z);
          SitesOnPlaquettes_[cnt][1] = Index(x+1,y,z);
          SitesOnPlaquettes_[cnt][2] = Index(x+1,y+1,z);
          SitesOnPlaquettes_[cnt][3] = Index(x,y+1,z);
          cnt++;
        }
      }
      //Plaquettes perpendicular to the plane z (y is constant!)
      for(int y=0;y<L_;y++) {
        for(int x=0;x<L_;x++) {
          SitesOnPlaquettes_[cnt][0] = Index(x,y,z+1);
          SitesOnPlaquettes_[cnt][1] = Index(x+1,y,z+1);
          SitesOnPlaquettes_[cnt][2] = Index(x+1,y,z);
          SitesOnPlaquettes_[cnt][3] = Index(x,y,z);
          cnt++;
        }
      }
    }
    //Replica
    for(int z=0;z<L_;z++){
      //Plaquettes lying on the plane z (z is constant!)
      for(int y=0;y<L_;y++) {
        for(int x=0;x<L_;x++) {
          SitesOnPlaquettes_[cnt][0] = Nplaqs_/2+Index(x,y,z);
          SitesOnPlaquettes_[cnt][1] = Nplaqs_/2+Index(x+1,y,z);
          SitesOnPlaquettes_[cnt][2] = Nplaqs_/2+Index(x+1,y+1,z);
          SitesOnPlaquettes_[cnt][3] = Nplaqs_/2+Index(x,y+1,z);
          cnt++;
        }
      }
      //Plaquettes perpendicular to the plane z (y is constant!)
      for(int y=0;y<L_;y++) {
        for(int x=0;x<L_;x++) {
          SitesOnPlaquettes_[cnt][0] = Nplaqs_/2+Index(x,y,z+1);
          SitesOnPlaquettes_[cnt][1] = Nplaqs_/2+Index(x+1,y,z+1);
          SitesOnPlaquettes_[cnt][2] = Nplaqs_/2+Index(x+1,y,z);
          SitesOnPlaquettes_[cnt][3] = Nplaqs_/2+Index(x,y,z);
          cnt++;
        }
      }
    }
  }

  //Indexing of coordinates
  int Index(int x, int y,int z) {
  
    if (x<0) x+= L_;
    if (x>=L_) x-= L_;
    if (y<0) y+= L_;
    if (y>=L_) y-= L_;
    if (z<0) z+= L_;
    if (z>=L_) z-= L_;
  
    return L_*L_*z + L_*y + x;
  
  }
  
  //Print Lattice Informations
  void Print() {
    std::cout << std::endl << std::endl << "Printing Indexing of Coordinates..." << std::endl << std::endl;

    for (int x=0; x<L_; x++) {
      for (int y=0; y<L_; y++) {
        for (int z=0; z<L_; z++) {
          std::cout << "Coordinate (" << x << "," << y << ","<<z<<")";
          std::cout << "  -->  Index: " << Index(x,y,z) << std::endl;
        }
      }
    }

    std::cout << std::endl << std::endl << std::endl;

    std::cout << "Printing coordinates of indices..." << std::endl << std::endl;
    
    for (int i=0; i<Nsites_; i++) {

        std::cout << "Index "<<i;
        std::cout << "  --> Coordinates: (";
        std::cout << Coordinates_[i][0] <<","<<Coordinates_[i][1]<<","<<Coordinates_[i][2]<<")";
        std::cout << std::endl;
    }
    std::cout << std::endl << std::endl << std::endl;
  
      std::cout << "Printing neighbors..." << std::endl << std::endl;
      
      for (int i=0; i<Nsites_; i++) {
  
          std::cout << "Index "<<i;
          std::cout << "  --> Neighbors: ";
          
          for (int j=0; j<6;j++) {
              std::cout << neighbours_[i][j] << " , ";
          }
          
          std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;
  
      std::cout << "Printing sites on plaquettes..." << std::endl << std::endl;
   
      for (int i=0; i<Nplaqs_; i++) {
  
          std::cout << "Plaquette "<<i;
          std::cout << "  --> Sites: ";
          
          for (int j=0; j<4;j++) {
              std::cout << SitesOnPlaquettes_[i][j] << " , ";
          }
          
          std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;
  }
  
};

#endif
