#include <iostream>
#include "TGeoTube_v.h"
#include "TGeoBBox_v.h"

#include "TBBStopWatch.h"

#include <boost/simd/sdk/simd/pack.hpp>
#include <boost/simd/include/functions/splat.hpp>
#include <boost/simd/include/functions/plus.hpp>
#include <boost/simd/include/functions/multiplies.hpp>
#include <boost/simd/include/functions/is_equal.hpp>
#include <boost/simd/include/functions/abs.hpp>
#include <boost/simd/include/functions/if_else.hpp>
#include <boost/simd/operator/functions/logical_and.hpp>
#include <boost/simd/include/functions/if_zero_else_one.hpp>
#include <boost/simd/include/functions/if_one_else_zero.hpp>
#include <boost/simd/include/functions/is_greater.hpp>
#include <boost/simd/include/functions/is_less.hpp>
#include  <boost/simd/memory/functions/load.hpp>
#include <boost/simd/memory/functions/aligned_load.hpp>
#include <boost/simd/memory/prefetch.hpp>

TGeoTube_v::TGeoTube_v(Double_t rmin, Double_t rmax, Double_t dz)
           :TGeoTube(rmin, rmax, dz)
{
// Default constructor specifying minimum and maximum radius
    SetShapeBit(TGeoShape::kGeoTube);
    SetTubeDimensions(rmin, rmax, dz);
    if ((fDz<0) || (fRmin<0) || (fRmax<0)) {
          SetShapeBit(kGeoRunTimeShape);
//    if (fRmax<=fRmin) SetShapeBit(kGeoInvalidShape);
//                //      printf("tube : dz=%f rmin=%f rmax=%f\n", dz, rmin, rmax);
              }
              ComputeBBox();
}

/*#ifndef VEC_EXTENSIONS
//_____________________________________________________________________________                                                                                                     
void TGeoTube_v::Contains_v(const StructOfCoord  & pointi, Bool_t * isin, Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      isin[i] = TGeoTube::Contains(point);
    }
}
#else
*/
// PUT VC CODE OR THE LIKE HERE
#ifdef BOOST_SIMD
void TGeoTube_v::Contains_v(const StructOfCoord  & pointi, Bool_t * isin, Int_t np) const
{
StopWatch tmr;
double accT=0.0;

 typedef boost::simd::pack<double> p_t; //Creating a simd vector of doubles
 const int lengthOfVector=sizeof(p_t)/sizeof(double); //Calculating length of Vector
 typedef boost::simd::pack<double> pd_m; //used as mask to receive the logical vector for doubles
 //typedef boost::simd::pack<int,lengthOfVector> p_m;

 //p_m zero_v(0); //kFALSE
 //p_m one_v(1);  //kTRUE

 //p_m one_i(1);
 //p_m zero_i(0);
 //std::cout<<"Length of double vector : "<<lengthOfVector<<std::endl;
 const int tailsize =np % lengthOfVector;
 p_t fDz_v(fDz);
 p_t fRmin_v(fRmin),fRmax_v(fRmax);

 //prefetching data outside
 boost::simd::prefetch_read<1>(pointi.x);
 boost::simd::prefetch_read<1>(pointi.y);
 boost::simd::prefetch_read<1>(pointi.z); 

 for(unsigned int i = 0; i < np-tailsize; i += lengthOfVector)

    { 

        pd_m in_mask(1);
	pd_m in_mask1(1);
	pd_m in_mask2(1);
        
	//p_t x_v(&pointi.x[i]); // will copy a certain number of x's into vc vector x;
        //p_t y_v(&pointi.y[i]);
        //p_t z_v(&pointi.z[i]);

	//Trying prefetching	
	//boost::simd::prefetch_read<1>(&pointi.x[i+1]);
	//boost::simd::prefetch_read<1>(&pointi.y[i+1]);
	//boost::simd::prefetch_read<1>(&pointi.z[i+1]);
        
	//Tracking time to load data
	tmr.Start();
	p_t x_v=boost::simd::aligned_load<p_t>(&pointi.x[i],0); // loading the simd vectors X component
	p_t y_v=boost::simd::aligned_load<p_t>(&pointi.y[i],0); // Y component
	p_t z_v=boost::simd::aligned_load<p_t>(&pointi.z[i],0); // Z component
        tmr.Stop();
	accT+=tmr.getDeltaSecs();

	//in_mask1=if_else((boost::simd::abs(z_v) > fDz_v),zero_v,one_v );
	in_mask1=if_zero_else_one(boost::simd::abs(z_v) > fDz_v);
	p_t r2_v(0);
	r2_v = x_v*x_v + y_v*y_v;
        //in_mask2=if_else(( (r2_v <  fRmin_v*fRmin_v) || (r2_v >  fRmax_v*fRmax_v) ),zero_v,one_v);
        in_mask2=if_zero_else_one(( (r2_v <  fRmin_v*fRmin_v) || (r2_v >  fRmax_v*fRmax_v) ));
	
	//in_mask=if_else((in_mask1 && in_mask2) ,one_i,zero_i);
	in_mask=if_one_else_zero((in_mask1 && in_mask2));
  	//std::cout<<"Mask : "<<in_mask<<std::endl;	
	for(unsigned int j = 0; j < lengthOfVector; ++j)
              {
                isin[i+j] = in_mask[j];//?kTRUE:kFALSE;
              }
  
  }
std::cout<<"Data loading time : "<<accT<<std::endl;
// tail part
for(unsigned int i = 0; i < tailsize; ++i)
      {
        Double_t point[3] = {pointi.x[np-tailsize+i], pointi.y[np-tailsize+i], pointi.z[np-tailsize+i]};
        isin[np-tailsize+i] = TGeoTube::Contains(point);
      }


//p_t v1(2.45,1.2,4.5,3.5);
//p_t v2(1.69,2.1,5.4,1.3);
//std::cout<< ( v1 > v2 )<<std::endl;
}
#endif
//#endif


//_____________________________________________________________________________

void TGeoTube_v::Safety_v(const StructOfCoord & pointi, Bool_t in, Double_t * safety, Int_t np ) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      safety[i] = TGeoTube::Safety(point, in);
    }
}

//_____________________________________________________________________________                                                                                                     
void TGeoTube_v::DistFromInside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/,const Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoTube::DistFromInside(point, dir, 3, step[i], 0);
    }
}

//_____________________________________________________________________________                                                                                                    

void TGeoTube_v::DistFromOutside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/, const  Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
    for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoTube::DistFromOutside(point, dir, 3, step[i], 0);
    }

}


//_____________________________________________________________________________                                                                                                     
void TGeoTubeSeg_v::Contains_v(const StructOfCoord  & pointi, Bool_t * isin, Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      isin[i] = TGeoTubeSeg::Contains(point);
    }
}

//_____________________________________________________________________________

void TGeoTubeSeg_v::Safety_v(const StructOfCoord & pointi, Bool_t in, Double_t * safety, Int_t np ) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      safety[i] = TGeoTubeSeg::Safety(point, in);
    }

}
//_____________________________________________________________________________                                                                                                     
void TGeoTubeSeg_v::DistFromInside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/,const Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoTubeSeg::DistFromInside(point, dir, 3, step[i], 0);
    }
}

//_____________________________________________________________________________                                                                                                    

void TGeoTubeSeg_v::DistFromOutside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/, const  Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
    for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoTubeSeg::DistFromOutside(point, dir, 3, step[i], 0);
    }

}
