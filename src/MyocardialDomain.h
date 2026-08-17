/**
 *  OpenCCO implementation 
 *  Copyright (C) 2023 B. Kerautret;  Phuc Ngo, N. Passat H. Talbot and C. Jaquet
 *  Modifications Copyright (C) 2024-2026 S. L. Vajire, J. S. Choy, G. S. Kassab and
 *  L.-C. Lee (Michigan State University; California Medical Innovations Institute),
 *  as part of HCT-CCO, a coronary-specific extension of OpenCCO. This file has been
 *  modified from the original; see README for the list of changes. Modifications are
 *  released under the same GNU General Public License v3.
 *
 *  This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 **/

#pragma once

#if defined(DOMAIN_CONTROLLER_RECURSE)
#error Recursive header files inclusion detected in MyocardialDomain.h
#else // defined(DOMAIN_CONTROLLER_RECURSE)
/** Prevents recursive inclusion of headers. */
#define DOMAIN_CONTROLLER_RECURSE

#if !defined DOMAIN_CONTROLLER_H
/** Prevents repeated inclusion of headers. */
#define DOMAIN_CONTROLLER_H


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/io/readers/GenericReader.h"
#include <array>
#include <limits>
#include <algorithm>
#include <cmath>


/**
 * Classes to control the construction domain of the coronary aretery tree.
 *  It is used during the reconstruction using the following points:
 *  - myCenter: used to orient the construction according the domain.
 *  - bool isInside(const TPoint &p):
 *    Test if a point is inside the domain or not (a point can be outside after the optmisation).
 *  - TPoint  randomPoint():
 *    To get a point inside the domaine used in the tree construction.
 *  - bool checkNoIntersectDomain(const TPoint &pt1, const TPoint &pt2):
 Usefull to ensure is a whole segment is inside the domain.
 *  - maxDistantPointFromBorder() const:
 *    usefull to determine a starting point and help to construct the tree.
 *  - TPoint firtCandidatePoint() const :
 *    use to initiate the reconstruction of the tree.
 **/


/**
 *  Domain controller based on a ball
 */
template<int TDim>
class CircularDomainCtrl {
public:
    typedef DGtal::PointVector<TDim, double> TPoint;
    typedef DGtal::PointVector<TDim, double> TPointD;

    typedef enum {NO_UPDATE, UPDATED} UPDATE_RAD_TYPE ;
    UPDATE_RAD_TYPE myUpdateType = UPDATED;
    double myRadius {1.0};
    TPoint myCenter;
    
    // Constructor for ImplicitCirc type
    CircularDomainCtrl(){};
    
    CircularDomainCtrl(double radius, const TPoint &center)
    {
        myRadius = radius;
        myCenter = center;
    };
    
    virtual bool isInside(const TPoint &p)
    {
        return (myCenter-p).norm() < myRadius;
    }
    
    TPoint
    randomPoint() {
        bool found = false;
        TPoint p;
        while(!found){
            double ss = 0.0;
            for(unsigned int i = 0;i<TDim; i++ ){
                p[i] = ((double)rand() / RAND_MAX)*2.0*myRadius - myRadius;
            }
            found = isInside(p);
        }
        return p + myCenter;
    }
    /**
     * Get the supported domain of the tree. By default it is defined from the circle center.
     * If the domain if defined from a mask image, the center if computed from the imate center.
     */
    TPoint getDomainCenter() const{
        return myCenter;
    }
    
    /**
     * Check if the segment defined by two points intersect the domain.
     *
     * @param pt1 first point of the segment
     * @param pt2  second point of the segment
     */
    bool
    checkNoIntersectDomain(const TPoint &pt1, const TPoint &pt2)
    {
        return isInside(pt1) && isInside(pt2);
    }
    
    TPoint
    maxDistantPointFromBorder() const {
        return myCenter;
    }
    
    TPoint
    firtCandidatePoint() const {
        TPoint res;
        if (TDim == 2){
            res[0] = 0.0;
            res[1] = myRadius;
        }else  if (TDim == 3){
            res[0] = 0.0;
            res[1] = myRadius;
            res[2] = 0.0;
            
        }
        
        return res;
        
    }
    std::vector<std::vector< TPoint > >
    contours()
    {
        std::vector<std::vector< TPoint > > res;
        return res;
    }
    TPoint
    lowerBound()
    {
        TPoint p = TPoint::diagonal(myRadius*0.01);
        return myCenter - p;
    }
    TPoint
    upperBound()
    {
        TPoint p = TPoint::diagonal(myRadius*0.01);
        return myCenter + p;
    }
    
};
template<int TDim>
class SquareDomainCtrl: public CircularDomainCtrl<TDim>{
    
public:
    // Constructor for ImplicitCirc type
    SquareDomainCtrl(){};
    
    SquareDomainCtrl(double radius,
                     const typename CircularDomainCtrl<TDim>::TPoint &center)
    {
        CircularDomainCtrl<TDim>::myRadius = radius;
        CircularDomainCtrl<TDim>::myCenter = center;
    };
    bool isInside(const typename SquareDomainCtrl<TDim>::TPoint &p)
    {
        bool res = true;
        for (unsigned int i=0; i<TDim; i++){
            res = res && (CircularDomainCtrl<TDim>::myCenter-p)[i] < CircularDomainCtrl<TDim>::myRadius;
        }
        return res;
    }
    
};

/**
 *  Domain controller based on a Image Mask
 */
template<int TDim>
class ImageMaskDomainCtrl {
public:
    
    typedef DGtal::PointVector<TDim, double> TPoint;
    typedef  DGtal::PointVector<TDim, int> TPointI;
    typedef DGtal::SpaceND< TDim, int >   SpaceCT;
    typedef DGtal::HyperRectDomain<SpaceCT> DomCT;
    typedef DGtal::ImageContainerBySTLVector< DomCT, unsigned char> Image;
    typedef DGtal::ImageContainerBySTLVector< DomCT, double> ImageD;
    typedef typename DGtal::DigitalSetSelector<DomCT,
    DGtal::BIG_DS+
    DGtal::HIGH_BEL_DS>::Type TDGset;
    typedef enum {NO_UPDATE, UPDATED} UPDATE_RAD_TYPE ;
    UPDATE_RAD_TYPE myUpdateType = NO_UPDATE;

    
    TPoint myDomPtLow, myDomPtUpper;
    TPointI myCenter;
    double myRadius {1.0};
    double minDistInitSegment {5.0};
    
    // Specific attributes to ImageMaskDomainCtr
    int myMaskThreshold {128};
    unsigned int myNbTry {100};
    double myMinDistanceToBorder {5.0};

    // Cache of valid domain points used by randomPoint() fallback.
    // Rebuilt lazily when empty or when myMinDistanceToBorder changes.
    std::vector<TPointI> myCachedValidPoints;
    double myCachedValidPointsForDist {-1.0};

public:
    Image myImage {Image(DomCT())} ;
    ImageD myDistanceImage {ImageD(DomCT())};
    ImageMaskDomainCtrl(const ImageMaskDomainCtrl&) {
        std::cout << "copy domain!!" << std::endl;
    }
    ImageMaskDomainCtrl(){};

    double sampleDistance(const TPoint &pt) const {
        return interpolateDistance(pt);
    }

    double sampleDistance(const TPointI &pt) const {
        TPoint cont;
        for (unsigned int i = 0; i < TDim; ++i) {
            cont[i] = static_cast<double>(pt[i]);
        }
        return interpolateDistance(cont);
    }

    
    // Constructor
    ImageMaskDomainCtrl(const std::string &fileImgDomain,
                        int maskThreshold, TPointI ptRoot,
                        unsigned int nbTry=100): myNbTry{nbTry}
    {
        myImage = DGtal::GenericReader<Image>::import(fileImgDomain,myMaskThreshold);
        myDistanceImage = CoronaryGeometry::getImageDistance<Image,ImageD>(myImage,
                                                                      myMaskThreshold );
        if ( !isInside(ptRoot) ){
            DGtal::trace.warning() << "ImageMaskDomainCtrl: Initial point given as input is not in domain." << std::endl;
            DGtal::trace.warning() << "ImageMaskDomainCtrl: Using default value from the maximal distant point." << std::endl;
            myCenter = maxDistantPointFromBorder();
        }else{
            myCenter = ptRoot;
        }

    }
    
    // Constructor
    ImageMaskDomainCtrl(const std::string &fileImgDomain,
                        int maskThreshold, unsigned int nbTry=100):
                                                    myNbTry{nbTry}
    {
        myImage = DGtal::GenericReader<Image>::import(fileImgDomain,myMaskThreshold);
        myDistanceImage = CoronaryGeometry::getImageDistance<Image,ImageD>(myImage,
                                                                      myMaskThreshold );
        myCenter = maxDistantPointFromBorder();
        checkImageDomain();
    };
    
    
    
    
    
    TPointI
    randomPoint()
    {
        bool found = false;
        TPointI pMin = myImage.domain().lowerBound();
        TPointI pMax = myImage.domain().upperBound();
        TPointI dp = pMax - pMin;
        TPointI pCand;
        unsigned int n = 0;
        while(!found && n < myNbTry){
            for (unsigned int i = 0; i< TDim; i++){
                pCand[i] = pMin[i]+(rand()%dp[i]);
            }
            double dist = sampleDistance(pCand);
            found = myImage(pCand)>=myMaskThreshold &&
            dist >= myMinDistanceToBorder;
            n++;
        }
        if (found){
            return pCand;
        }
        if (n >= myNbTry){
            // Fallback: previously this scanned the domain in point order
            // and returned the FIRST valid point, which is deterministic
            // and spatially biased (always the same corner). Instead,
            // collect all valid points once (cached) and pick uniformly.
            if (myCachedValidPoints.empty() ||
                myCachedValidPointsForDist != myMinDistanceToBorder) {
                myCachedValidPoints.clear();
                for (auto p : myImage.domain()) {
                    double dist = sampleDistance(p);
                    if (myImage(p) >= myMaskThreshold &&
                        dist >= myMinDistanceToBorder) {
                        myCachedValidPoints.push_back(p);
                    }
                }
                myCachedValidPointsForDist = myMinDistanceToBorder;
            }
            if (!myCachedValidPoints.empty()) {
                unsigned int idx = static_cast<unsigned int>(rand()) %
                                   static_cast<unsigned int>(myCachedValidPoints.size());
                return myCachedValidPoints[idx];
            }
        }
        return TPointI();
    }
    bool isInside(const TPointI &p){
      return myImage.domain().isInside(p) && myImage(p) > myMaskThreshold;
    }
    /**
     * Check if the segment defined by two points intersect the domain.
     *
     * @param pt1 first point of the segment
     * @param pt2  second point of the segment
     */
    bool
    checkNoIntersectDomain(const TPointI &pt1, const TPointI &pt2) const
    {
        if (!myImage.domain().isInside(pt1) ||
            !myImage.domain().isInside(pt2)){
            return false;
        }
        TPoint start, end;
        for(unsigned int i=0; i<TDim; i++ ){
            start[i] = static_cast<double>(pt1[i]);
            end[i] = static_cast<double>(pt2[i]);
        }
        TPoint dir = end - start;
        double length = dir.norm();
        if (length == 0.0){
            return sampleDistance(start) >= myMinDistanceToBorder;
        }
        dir /= length;
        for (double s = 0.0; s <= length; s += 1.0){
            TPoint p = start + dir * s;
            if (sampleDistance(p) < myMinDistanceToBorder){
                return false;
            }
        }
        return true;
    }
    
    TPointI
    maxDistantPointFromBorder() const {
        double m = 0.0;
        TPointI pM;
        for(auto p: myDistanceImage.domain()) {
            if (myDistanceImage(p) > m ){m = myDistanceImage(p);
                for (unsigned int i=0; i<TDim; i++){
                    pM[i] = p[i];
                }
            }
        }
        return pM;
    }
    TPoint
    firtCandidatePoint() const {
        TPointI res;
        bool find = searchRootFarthest(std::max(myDistanceImage(myCenter)/2.0, minDistInitSegment), res);
        assert(find);
        return res;
    }
    
    std::vector<std::vector< TPointI > >
    contours()
    {
        std::vector<std::vector< TPointI > > res;
        return res;
    }
    TPointI
    lowerBound()
    {
        
        return myImage.domain().lowerBound();
    }
    TPointI
    upperBound()
    {
        return myImage.domain().upperBound();
    }
    
private:
    double interpolateDistance(const TPoint &pt) const {
        const TPointI lower = myImage.domain().lowerBound();
        const TPointI upper = myImage.domain().upperBound();
        for (unsigned int i = 0; i < TDim; ++i) {
            if (pt[i] < lower[i] || pt[i] > upper[i]) {
                return -std::numeric_limits<double>::infinity();
            }
        }
        std::array<int, TDim> idx0;
        std::array<int, TDim> idx1;
        std::array<double, TDim> frac;
        for (unsigned int i = 0; i < TDim; ++i) {
            double lo = static_cast<double>(lower[i]);
            double hi = static_cast<double>(upper[i]);
            double coord = std::min(std::max(pt[i], lo), hi);
            double base = std::floor(coord);
            int i0 = static_cast<int>(base);
            if (i0 < lower[i]) {
                i0 = lower[i];
                base = static_cast<double>(i0);
            }
            if (i0 >= upper[i]) {
                if (upper[i] == lower[i]) {
                    i0 = upper[i];
                } else {
                    i0 = upper[i] - 1;
                }
                base = static_cast<double>(i0);
            }
            int i1 = std::min(i0 + 1, upper[i]);
            double t = coord - base;
            if (i1 == i0) {
                t = 0.0;
            }
            idx0[i] = i0;
            idx1[i] = i1;
            frac[i] = t;
        }
        double value = 0.0;
        const unsigned int cornerCount = 1u << TDim;
        for (unsigned int mask = 0; mask < cornerCount; ++mask) {
            double weight = 1.0;
            TPointI corner;
            for (unsigned int d = 0; d < TDim; ++d) {
                if (mask & (1u << d)) {
                    corner[d] = idx1[d];
                    weight *= frac[d];
                } else {
                    corner[d] = idx0[d];
                    weight *= (1.0 - frac[d]);
                }
            }
            value += weight * myDistanceImage(corner);
        }
        return value;
    }
    
private:
    // internal method
    
    bool
    searchRootFarthest(const double & d, TPointI &ptRoot ) const {
        typedef DGtal::SpaceND<TDim, int> Space;
        DomCT aDom;
        TDGset sPts = CoronaryGeometry::pointsOnSphere<TPointI, TDGset>(myCenter, d);
        for (const TPointI &p : sPts){
            if (checkNoIntersectDomain(p, myCenter)){
                for(unsigned int i = 0; i < TDim; i++){
                    ptRoot[i] = p[i];
                }
                return true;
            }
        }
        return false;
    }
    
    void checkImageDomain(){
        bool isOk = false;
        // Check if at least one pixel of with foreground value exist:
        for (auto p: myImage.domain()){
            if (myImage(p) >= myMaskThreshold){
                isOk = true;
                break;
            }
        }
        if (!isOk){
            std::cout << "ImageMaskDomainCtrl: domain non valid since no point are inside the mask image." <<  std::endl;
        }
    }
    
};

template<>
std::vector<std::vector< typename ImageMaskDomainCtrl<2>::TPointI > >
ImageMaskDomainCtrl<2>::contours()
{
    typedef  DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> TImage;
    typedef DGtal::functors::IntervalThresholder<typename TImage::Value> Binarizer;
    DGtal::Z2i::KSpace ks;
    if(! ks.init( myImage.domain().lowerBound(),
                 myImage.domain().upperBound(), true )){
        DGtal::trace.error() << "Problem in KSpace initialisation"<< std::endl;
    }
    
    Binarizer b(myMaskThreshold, 255);
    DGtal::functors::PointFunctorPredicate<TImage,Binarizer> predicate(myImage, b);
    DGtal::trace.info() << "DGtal contour extraction from thresholds ["<<  myMaskThreshold << "," << 255 << "]" ;
    DGtal::SurfelAdjacency<2> sAdj( true );
    std::vector<std::vector< typename ImageMaskDomainCtrl<2>::TPointI > > vectContoursBdryPointels;
    DGtal::Surfaces<DGtal::Z2i::KSpace>::extractAllPointContours4C( vectContoursBdryPointels, ks, predicate, sAdj );
    return vectContoursBdryPointels;
}





///////////////////////////////////////////////////////////////////////////////

#endif // !defined DOMAIN_CONTROLLER_H

#undef DOMAIN_CONTROLLER_RECURSE
#endif // else defined(DOMAIN_CONTROLLER_RECURSE)





