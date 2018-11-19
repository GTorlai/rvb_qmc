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
  
  }
  
};

#endif
