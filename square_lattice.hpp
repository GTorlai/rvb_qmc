#ifndef SQUARELATTICE_HPP 
#define SQUARELATTICE_HPP 

#include <fstream>

class SquareLattice {

private:

  int L_;
  int Lx_,Ly_;
  int Nsites_;
  int Nbonds_;

public:
        
  std::vector<std::vector<int> > Coordinates_;        // (x,y) coordinates of each lattice site
  std::vector<std::vector<int> > neighbours_;         // 4 neighbours of each lattice site
  std::vector<std::vector<int> > BondsOnPlaquettes_;  // 4 bonds in each plaquette
  std::vector<std::vector<int> > SitesOnBonds_;       // 2 sites on each bond
  std::vector<std::vector<int> > BondsOnSites_;       // 4 bonds on each site
  std::vector<std::vector<int> > SitesOnPlaquettes_;  // 4 sites on each plaquette
  
  SquareLattice(int L){
    
    L_ = L;
    Nsites_ = 2*L*L;//Lx_*Ly_;  // The factor 2 is for the two replicas
    Nbonds_ = 2*Nsites_;
    Init();
  
  }  
 
  inline int LinSize() const {return L_;}

  inline int Nsites() const {return Nsites_;}

  inline int Nbonds() const {return Nbonds_;}

  void Init(){
    int counter[2]={0,0};
    
    Coordinates_.resize(Nsites_,std::vector<int>(2));
    neighbours_.resize(Nsites_,std::vector<int>(4));
    BondsOnPlaquettes_.resize(Nsites_,std::vector<int>(4));
    SitesOnBonds_.resize(Nbonds_,std::vector<int>(2));
    BondsOnSites_.resize(Nsites_,std::vector<int>(4));
    SitesOnPlaquettes_.resize(Nsites_,std::vector<int>(4));

    //Coordinates
    for (int i=1;i<Nsites_/2;i++){
      if ( (counter[0]+1) % L_ == 0){//end of x-row
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
    // Real
    for(int i=0;i<Nsites_/2;i++) {
      neighbours_[i][0]=Index(Coordinates_[i][0]+1,Coordinates_[i][1]);
      neighbours_[i][1]=Index(Coordinates_[i][0]  ,Coordinates_[i][1]+1);
      neighbours_[i][2]=Index(Coordinates_[i][0]-1,Coordinates_[i][1]);
      neighbours_[i][3]=Index(Coordinates_[i][0]  ,Coordinates_[i][1]-1);
    }
    // Replica
    for(int i=0;i<Nsites_/2;i++) {
      neighbours_[Nsites_/2+i][0]=Nsites_/2+Index(Coordinates_[i][0]+1,Coordinates_[i][1]);
      neighbours_[Nsites_/2+i][1]=Nsites_/2+Index(Coordinates_[i][0]  ,Coordinates_[i][1]+1);
      neighbours_[Nsites_/2+i][2]=Nsites_/2+Index(Coordinates_[i][0]-1,Coordinates_[i][1]);
      neighbours_[Nsites_/2+i][3]=Nsites_/2+Index(Coordinates_[i][0]  ,Coordinates_[i][1]-1);
    }
    
    //4 bonds on a plaquette
    //Real
    for(int x=0;x<L_;x++) {
      for(int y=0;y<L_;y++) {
        BondsOnPlaquettes_[Index(x,y)][0]=2*Index(x,y);
        BondsOnPlaquettes_[Index(x,y)][3]=2*Index(x,y)+1;
        BondsOnPlaquettes_[Index(x,y)][1]=2*Index(x+1,y)+1;
        BondsOnPlaquettes_[Index(x,y)][2]=2*Index(x,y+1);
      }
    }
    //Replica
    for(int x=0;x<L_;x++) {
      for(int y=0;y<L_;y++) {
        BondsOnPlaquettes_[Nsites_/2+Index(x,y)][0]=Nbonds_/2+2*Index(x,y);
        BondsOnPlaquettes_[Nsites_/2+Index(x,y)][3]=Nbonds_/2+2*Index(x,y)+1;
        BondsOnPlaquettes_[Nsites_/2+Index(x,y)][1]=Nbonds_/2+2*Index(x+1,y)+1;
        BondsOnPlaquettes_[Nsites_/2+Index(x,y)][2]=Nbonds_/2+2*Index(x,y+1);
      }
    }

    //2 sites on a bond
    //Real
    for(int x=0;x<L_;x++) {
      for(int y=0;y<L_;y++) {
        if((x+y)%2==0){
          // Horizontal
          SitesOnBonds_[2*Index(x,y)][0] = Index(x,y);
          SitesOnBonds_[2*Index(x,y)][1] = Index(x+1,y);
          // Vertical
          SitesOnBonds_[2*Index(x,y)+1][0] = Index(x,y);
          SitesOnBonds_[2*Index(x,y)+1][1] = Index(x,y+1);
        }
        else{
          // Horizontal
          SitesOnBonds_[2*Index(x,y)][1] = Index(x,y);
          SitesOnBonds_[2*Index(x,y)][0] = Index(x+1,y);
          // Vertical
          SitesOnBonds_[2*Index(x,y)+1][1] = Index(x,y);
          SitesOnBonds_[2*Index(x,y)+1][0] = Index(x,y+1);
        }
      }
    }
    //Replica
    for(int x=0;x<L_;x++) {
      for(int y=0;y<L_;y++) {
        if((x+y)%2==0){
          SitesOnBonds_[Nbonds_/2+2*Index(x,y)][0]   = Nsites_/2+Index(x,y);
          SitesOnBonds_[Nbonds_/2+2*Index(x,y)][1]   = Nsites_/2+Index(x+1,y);
          SitesOnBonds_[Nbonds_/2+2*Index(x,y)+1][0] = Nsites_/2+Index(x,y);
          SitesOnBonds_[Nbonds_/2+2*Index(x,y)+1][1] = Nsites_/2+Index(x,y+1);
        }
        else{
          SitesOnBonds_[Nbonds_/2+2*Index(x,y)][1]   = Nsites_/2+Index(x,y);
          SitesOnBonds_[Nbonds_/2+2*Index(x,y)][0]   = Nsites_/2+Index(x+1,y);
          SitesOnBonds_[Nbonds_/2+2*Index(x,y)+1][1] = Nsites_/2+Index(x,y);
          SitesOnBonds_[Nbonds_/2+2*Index(x,y)+1][0] = Nsites_/2+Index(x,y+1);
        }
      }
    }

    //Sites on Plaquettes
    //Real
    for(int x=0;x<L_;x++) {
      for(int y=0;y<L_;y++) {
        SitesOnPlaquettes_[Index(x,y)][0] = Index(x,y); 
        SitesOnPlaquettes_[Index(x,y)][1] = Index(x+1,y);
        SitesOnPlaquettes_[Index(x,y)][2] = Index(x+1,y+1); 
        SitesOnPlaquettes_[Index(x,y)][3] = Index(x,y+1);
      }
    }
    //Replica
    for(int x=0;x<L_;x++) {
      for(int y=0;y<L_;y++) {
        SitesOnPlaquettes_[Nsites_/2+Index(x,y)][0] = Nsites_/2+Index(x,y); 
        SitesOnPlaquettes_[Nsites_/2+Index(x,y)][1] = Nsites_/2+Index(x+1,y);
        SitesOnPlaquettes_[Nsites_/2+Index(x,y)][2] = Nsites_/2+Index(x+1,y+1); 
        SitesOnPlaquettes_[Nsites_/2+Index(x,y)][3] = Nsites_/2+Index(x,y+1);
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
  
  //Print Lattice Informations
  void Print() {
  
      std::cout << std::endl << std::endl << "Printing Indexing of Coordinates..." << std::endl << std::endl;
      std::cout << "Printing neighbors..." << std::endl << std::endl;
      
      for (int i=0; i<Nsites_; i++) {
  
          std::cout << "Index "<<i;
          std::cout << "  --> Neighbors: ";
          
          for (int j=0; j<4;j++) {
              std::cout << neighbours_[i][j] << " , ";
          }
          
          std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;
      
      std::cout << "Printing Bonds on plaquettes..." << std::endl << std::endl;
   
      for (int i=0; i<Nsites_; i++) {
  
          std::cout << "Plaquette "<<i;
          std::cout << "  --> Bonds: ";
          
          for (int j=0; j<4;j++) {
              std::cout << BondsOnPlaquettes_[i][j] << " , ";
          }
          
          std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;

      std::cout << "Printing sites on plaquettes..." << std::endl << std::endl;
   
      for (int i=0; i<Nsites_; i++) {
  
          std::cout << "Plaquette "<<i;
          std::cout << "  --> Sites: ";
          
          for (int j=0; j<4;j++) {
              std::cout << SitesOnPlaquettes_[i][j] << " , ";
          }
          
          std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;
      
      std::cout << "Printing sites on Bonds..." << std::endl << std::endl;
   
      for (int i=0; i<Nbonds_; i++) {
  
          std::cout << "Link "<<i;
          std::cout << "  --> Sites: ";
          
          for (int j=0; j<2;j++) {
              std::cout << SitesOnBonds_[i][j] << " , ";
          }
          
          std::cout << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;
  
  }
  
};

#endif
