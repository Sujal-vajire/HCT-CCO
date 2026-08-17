/**
 *  HctCco: main program to generate 2D tree from OpenCCO implementation
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

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

#include <iostream>
#include <fstream>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/Mesh.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/images/IntervalForegroundPredicate.h"

#include "CLI11.hpp"

#include "HctccoModel.h"
#include "CoronaryGeometry.h"
#include "CoronaryGrowth.h"
#include "CoronaryTreeIO.h"

#include "MyocardialDomain.h"

#ifdef WITH_VISU3D_QGLVIEWER
#include "DGtal/io/viewers/Viewer3D.h"
#endif

// Forward declaration for wall analysis
template<typename TTree>
void analyzeWallDistribution(const TTree& tree, const std::string& xmlFile, const std::string& volFile);
/**
 * @brief main function call
 *
 */

/**
 * Function to construct the tree with the help ConstructionHelpers by using a domain of reconstruction defined fram an image (ImageMaskDomainCtrl)l
 */
template<typename TTree>
void
constructTreeMaskDomain(TTree &aTree,
                       bool verbose,
                       bool hasImportedTree = false)
{
    clock_t start, end;
    start = clock();
    // Only seed a fresh root when we have NOT imported one from XML.
    // Previously this path unconditionally re-seeded the root, which
    // overwrote the imported root's coordinate and produced growth that
    // spatially drifted away from the imported base tree's orientation.
    if (!hasImportedTree) {
        CoronaryGrowth::initFirtElemTree(aTree, verbose);
    }
    CoronaryGrowth::expandTree(aTree, verbose);
    end = clock();
    printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
}


template<typename TTree>
void
constructTreeImplicitDomain(TTree &aTree, std::string outputMeshName,
                            std::string exportXMLName, std::string importXMLName,
                            std::string exportDatName, bool verbose,
                            bool display3D, bool postAnalyzeWall)
{
    clock_t start, end;
    start = clock();
    if (importXMLName == "")
        CoronaryGrowth::initFirtElemTree(aTree, verbose);
    CoronaryGrowth::expandTree(aTree, verbose);
    end = clock();
    printf ("Execution time: %0.8f sec\n", ((double) end - start)/CLOCKS_PER_SEC);
    if (exportDatName != "") exportDat(aTree, exportDatName);
    if (exportXMLName != "") CoronaryTreeIO::writeTreeToXml(aTree,
                                                        exportXMLName.c_str());

    CoronaryTreeIO::writeTreeToXml(aTree, "tree_3D.xml");
    exportResultingMesh(aTree, outputMeshName);
    #ifdef WITH_VISU3D_QGLVIEWER
        if (display3D) display3DTree(aTree);
    #endif
    
    // Note: Wall analysis only available for image domains (requires .vol file)
    if (postAnalyzeWall) {
        std::cout << "Warning: Wall analysis (-p) only available with image domains (--organDomain option)" << std::endl;
    }
   }




template <typename TTree>
void exportResultingMesh(const TTree &tree, std::string outName)
{
    unsigned int i = 0;
    double thickness = 1;
    // Export 3D mesh of the tree
    DGtal::Mesh<DGtal::Z3i::RealPoint> aMesh(true);
    for (auto s : tree.myVectSegments) {
      // test if the segment is the root or its parent we do not display (already done).
      if (s.myIndex == 0 || s.myIndex == 1)
        continue;
      DGtal::Z3i::RealPoint distal = s.myCoordinate;
      DGtal::Z3i::RealPoint proxital = tree.myVectSegments[tree.myVectParent[s.myIndex]].myCoordinate;
      auto v = {distal, proxital};
      DGtal::Mesh<DGtal::Z3i::RealPoint>::createTubularMesh(aMesh, v, tree.myVectSegments[s.myIndex].myRadius*thickness, 0.05);
      i++;
    }
    aMesh >> outName;
}

template <typename TTree>
void
exportDat(const TTree &tree, std::string outName)
{
    unsigned int i =0;
    double thickness = 1;
    std::ofstream fout;
    fout.open(outName.c_str());
    for (auto s : tree.myVectSegments) {
      // test if the segment is the root or its parent we do not display (already done).
      if (s.myIndex == 0 || s.myIndex == 1)
        continue;
      DGtal::Z3i::RealPoint distal = s.myCoordinate;
      DGtal::Z3i::RealPoint proxital = tree.myVectSegments[tree.myVectParent[s.myIndex]].myCoordinate;
      fout << distal[0] << " " << distal[1] << " " << distal[2] << " ";
      fout << proxital[0] << " " << proxital[1] << " " << proxital[2] << " ";
      fout << tree.myVectSegments[s.myIndex].myRadius*thickness << std::endl;
      i++;
    }
    fout.close();
}

#ifdef WITH_VISU3D_QGLVIEWER
template<typename TTree>
int
display3DTree(const TTree &tree )
{
    int arg = 0;
      QApplication application(arg,NULL);
      typedef DGtal::Viewer3D<> MyViewer;
      MyViewer viewer;
      viewer.show();
      unsigned int i = 0;
      double thickness = 1;
      viewer << DGtal::CustomColors3D(DGtal::Color(0,0,250),DGtal::Color(0,0,250));
      DGtal::Z3i::RealPoint p1 = tree.myVectSegments[1].myCoordinate;
      DGtal::Z3i::RealPoint p2 = tree.myVectSegments[tree.myVectParent[tree.myVectSegments[1].myIndex]].myCoordinate;
      viewer.addBall(p2,tree.myVectSegments[1].myRadius);
      viewer << DGtal::CustomColors3D(DGtal::Color(0,250,0),DGtal::Color(0,250,0));
      viewer.addCylinder(p1,p2,tree.myVectSegments[1].myRadius*thickness);
      
      for (auto s : tree.myVectSegments) {
        // test if the segment is the root or its parent we do not display (already done).
        if (s.myIndex == 0 || s.myIndex == 1)
          continue;
        DGtal::Z3i::RealPoint distal = s.myCoordinate;
        DGtal::Z3i::RealPoint proxital = tree.myVectSegments[tree.myVectParent[s.myIndex]].myCoordinate;
        viewer << DGtal::CustomColors3D(DGtal::Color(250,0,0),DGtal::Color(250,0,0));
        viewer.addBall(distal,tree.myVectSegments[s.myIndex].myRadius);
        viewer << DGtal::CustomColors3D(DGtal::Color(0,250,0),DGtal::Color(0,250,0));
        //viewer.addBall(distal,1);
        viewer.addCylinder(distal,proxital,tree.myVectSegments[s.myIndex].myRadius*thickness);
        i++;
      }
      
      viewer<< MyViewer::updateDisplay;
      return application.exec();
    }
#endif


void runToolsOnXML(std::string exportXMLName, std::string domainVol) {
    // Consolidated post-processing in a single call: diameter-defined Strahler
    // ordering + order-wise diameter transform + hemodynamics + final .vtp / .obj
    // / plots + CSVs. See postprocessing/coronary_postprocess.py. POSTPROCESS_SCRIPT set by CMake.
    // The domain mask (if any) lets the hemodynamics scale flow to physiological perfusion.
#ifdef POSTPROCESS_SCRIPT
    std::string prefix = exportXMLName;
    size_t dot = prefix.rfind(".xml");
    if (dot != std::string::npos) prefix = prefix.substr(0, dot);
    std::stringstream cmd;
    cmd << "python3 \"" << POSTPROCESS_SCRIPT << "\" \"" << exportXMLName
        << "\" \"" << prefix << "\" \"" << domainVol << "\"";
    system(cmd.str().c_str());
#else
    std::cerr << "[postprocess] POSTPROCESS_SCRIPT not set at build time; skipping "
                 "post-processing of " << exportXMLName << std::endl;
#endif
}

/**
 * Wall Distribution Analysis Function
 * Analyzes vessel distribution across LV wall thickness (inner 50% vs outer 50%)
 */
template<typename TTree>
void analyzeWallDistribution(const TTree& tree, const std::string& xmlFile, const std::string& volFile) {
    std::cout << "\n=== Wall Distribution Analysis ===" << std::endl;
    std::cout << "Tree file: " << xmlFile << std::endl;
    std::cout << "Domain file: " << volFile << std::endl;
    
    try {
        // Load domain image
        typedef DGtal::ImageContainerBySTLVector<typename TTree::DomCT, unsigned char> ImageType;
        typedef DGtal::functors::IntervalForegroundPredicate<ImageType> Binarizer;
        
        ImageType domain = DGtal::VolReader<ImageType>::importVol(volFile);
        std::cout << "Domain loaded: " << domain.domain().size() << " voxels" << std::endl;
        
        // Calculate wall thickness using distance transform
        typedef DGtal::ExactPredicateLpSeparableMetric<typename TTree::SpaceCT, 2> L2Metric;
        typedef DGtal::DistanceTransformation<typename TTree::SpaceCT, Binarizer, L2Metric> DTL;
        
        L2Metric l2metric;
        Binarizer binarizer(domain, 128, 255);  // Tissue > 128
        DTL distanceTransform(&domain.domain(), &binarizer, &l2metric);
        
        // Find maximum distance (wall thickness)
        int maxThickness = 0;
        for (auto point : domain.domain()) {
            if (domain(point) > 128) {  // Tissue voxel
                int dist = distanceTransform(point);
                maxThickness = std::max(maxThickness, dist);
            }
        }
        
        std::cout << "Maximum wall thickness: " << maxThickness << " voxels" << std::endl;
        
        // Calculate erosion steps (50% of thickness)
        int erosionSteps = std::max(1, maxThickness / 2);
        std::cout << "Erosion steps for 50% split: " << erosionSteps << std::endl;
        
        // Create inner region by erosion
        ImageType innerRegion = domain;
        for (int step = 0; step < erosionSteps; step++) {
            ImageType temp(innerRegion.domain());
            
            // Morphological erosion: remove voxels that have any neighbor ≤ 128
            for (auto point : innerRegion.domain()) {
                if (innerRegion(point) > 128) {  // Current voxel is tissue
                    bool hasNonTissueNeighbor = false;
                    
                    // Check 6-connected neighbors
                    for (int dim = 0; dim < 3; dim++) {
                        for (int dir = -1; dir <= 1; dir += 2) {
                            auto neighbor = point;
                            neighbor[dim] += dir;
                            
                            if (innerRegion.domain().isInside(neighbor)) {
                                if (innerRegion(neighbor) <= 128) {
                                    hasNonTissueNeighbor = true;
                                    break;
                                }
                            } else {
                                hasNonTissueNeighbor = true;  // Outside domain
                                break;
                            }
                        }
                        if (hasNonTissueNeighbor) break;
                    }
                    
                    if (hasNonTissueNeighbor) {
                        temp.setValue(point, 0);  // Erode this voxel
                    } else {
                        temp.setValue(point, innerRegion(point));  // Keep this voxel
                    }
                } else {
                    temp.setValue(point, innerRegion(point));  // Copy non-tissue
                }
            }
            innerRegion = temp;
        }
        
        // Count inner region voxels
        int innerVoxels = 0;
        for (auto point : innerRegion.domain()) {
            if (innerRegion(point) > 128) {
                innerVoxels++;
            }
        }
        std::cout << "Inner region: " << innerVoxels << " voxels" << std::endl;
        
        // Classify vessels
        int innerCount = 0, outerCount = 0;
        int totalSegments = tree.myVectSegments.size();
        
        for (unsigned int i = 1; i < totalSegments; i++) {  // Skip root segment (index 0)
            auto vesselPos = tree.myVectSegments[i].myCoordinate;
            
            // Convert to voxel coordinates
            typename TTree::TPoint voxelPos;
            for (int dim = 0; dim < 3; dim++) {
                voxelPos[dim] = static_cast<int>(std::round(vesselPos[dim]));
            }
            
            // Check if vessel position is valid and in tissue
            if (domain.domain().isInside(voxelPos) && domain(voxelPos) > 128) {
                if (innerRegion(voxelPos) > 128) {
                    innerCount++;  // Vessel in inner 50% of wall
                } else {
                    outerCount++;  // Vessel in outer 50% of wall
                }
            }
        }
        
        // Calculate and display results
        int analyzedVessels = innerCount + outerCount;
        if (analyzedVessels > 0) {
            double innerPercent = (double)innerCount / analyzedVessels * 100.0;
            double outerPercent = (double)outerCount / analyzedVessels * 100.0;
            
            std::cout << "\nResults:" << std::endl;
            std::cout << "Total analyzed vessels: " << analyzedVessels << " (of " << (totalSegments-1) << " total)" << std::endl;
            std::cout << "Inner 50% of wall contains: " << innerPercent << "% (" << innerCount << " segments)" << std::endl;
            std::cout << "Outer 50% of wall contains: " << outerPercent << "% (" << outerCount << " segments)" << std::endl;
            
            if (innerPercent > outerPercent) {
                std::cout << "Analysis: Vessels are more concentrated in the inner half of the wall" << std::endl;
            } else if (outerPercent > innerPercent) {
                std::cout << "Analysis: Vessels are more concentrated in the outer half of the wall" << std::endl;
            } else {
                std::cout << "Analysis: Vessels are evenly distributed across wall thickness" << std::endl;
            }
        } else {
            std::cout << "No vessels found in tissue domain!" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "Error in wall analysis: " << e.what() << std::endl;
    }
    
    std::cout << "Analysis complete.\n" << std::endl;
}

int main(int argc, char **argv)
{
  srand ((int) time(NULL));
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  app.description("Generated a 3D tree using the CCO algorithm with coronary-specific constraints.\n"
                  "By default it generates a 3D mesh.\n\n"
                  "TWO SIMULATION MODES:\n"
                  "1. Normal Mode (default): 3-stage approach\n"
                  "   - Stage 1 (0-10%): Connect anywhere in tree\n"
                  "   - Stage 2 (10-35%): Connect to leaf + 3-6 generations above\n"
                  "   - Stage 3 (35-100%): Connect to leaf + immediate parent only\n\n"
                  "2. Extension Mode (--extensionModel): 2-stage approach\n"
                  "   - Stage 1 (0-25%): Connect to leaf + 1-3 generations above\n"
                  "   - Stage 2 (25-100%): Connect to leaf + immediate parent only\n\n"
                  "CONSTRAINTS:\n"
                  "- Maximum branch length: 25 units (prevents unrealistic long vessels)\n"
                  "- Minimum branch radius: 0.005mm\n"
                  "- No artificial angle constraints (natural optimization)\n"
                  "\nExamples:\n"
                  "  ./HctCco -n 1000 -x coronary_tree.xml\n"
                  "  ./HctCco -n 1000 --extensionModel -x extended_tree.xml\n"
                  "  ./HctCco -n 1000 -d myocardium.vol -x tree.xml -w");
  int nbTerm {500};
  double aPerf {20000};
  double gamma {3.0};
  double rootRadiusMax {1.588};
  double termFlowCV {0.0};
  double sideBranchBias {0.0};

  double minDistanceToBorder {5.0};
  bool verbose {false};
  bool tools {true};   // post-processing runs by default (see --no-postprocess)
  bool display3D {false};
  bool extensionMode {false};
  bool postAnalyzeWall {false};
  bool anywhereMode {false};
    
  std::string nameImgDom {""};
  std::string outputMeshName {"result.off"};
  std::string exportDatName {""};
  std::string exportXMLName {""};
  std::string importXMLName {""};
  std::string toolsDirName {""};
  std::vector<int> postInitV {-1,-1,-1};
  bool squaredImplDomain {false};

  app.add_option("-n,--nbTerm,1", nbTerm, "Set the number of terminal segments.", true);
  app.add_option("-a,--aPerf,2", aPerf, "The value of perfusion volume.", true);
  app.add_option("-g,--gamma", gamma, "Bifurcation (Murray) exponent. Coronary literature ~2.3-3.0; lower values give a faster taper.", true);
  app.add_option("--rootRadiusMax", rootRadiusMax, "Cap on root radius (mm) anchoring absolute scale. Default 1.588 = literature order-11 radius. Set <=0 to disable.", true);
  app.add_option("--termFlowCV", termFlowCV, "Coefficient of variation of per-terminal perfusion flow (0=uniform classic CCO). >0 injects territory asymmetry to raise the Strahler diameter ratio toward coronary literature.", true);
  app.add_option("--sideBranchBias", sideBranchBias, "Topological asymmetry (0=classic min-volume). >0 prefers attaching new terminals as side branches on large high-flow vessels, raising the Strahler diameter/bifurcation ratios toward coronary morphometry. Try 0.3-1.0.", true);
  app.add_option("--organDomain,-d", nameImgDom, "Define the organ domain using a mask image (organ=255).");
  app.add_option("-m,--minDistanceToBorder", minDistanceToBorder, "Set the minimal distance to border. Works only  with option --organDomain else it has not effect", true);
  app.add_option("-o,--outputName", outputMeshName, "Output the 3D mesh into OFF format", true);
  app.add_option("-e,--export", exportDatName, "Output the 3D mesh into text file", true);
  app.add_option("-x,--exportXML", exportXMLName, "Output the resulting graph as xml file", true);
  app.add_option("-i,--importXML", importXMLName, "Segmented epicardial base-tree XML. The tree is initialized from this base tree and grown (extended) beyond it into the myocardial domain.", true);
  app.add_flag("-s,--squaredDom",squaredImplDomain , "Use a squared implicit domain instead a sphere (is used only without --organDomain)");
  auto pInit = app.add_option("-p,--posInit", postInitV, "Initial position of root, if not given the position of point is determined from the image center")
    ->expected(3);

#ifdef WITH_VISU3D_QGLVIEWER
  app.add_flag("--view", display3D, "display 3D view using QGLViewer");
#endif
  app.add_flag("-v,--verbose", verbose);
  app.add_flag("--postprocess,!--no-postprocess", tools, "Run post-processing (DDS ordering + diameter transform + final .vtp/.obj/plots) automatically after generation. Default on.");
  app.add_flag("--extensionModel", extensionMode, "Skip stage 1, use modified 2-stage approach (leaf + 1-3 generations, then leaf + immediate parent)");
  app.add_flag("-w,--wallAnalysis", postAnalyzeWall, "Analyze vessel distribution across wall thickness (inner 50% vs outer 50%)");
  app.add_flag("--anywhereModel", anywhereMode, "Unrestricted model: allow connecting anywhere (disables coronary restrictions)");
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  DGtal::Z3i::Point ptRoot(postInitV[0], postInitV[1], postInitV[2]);

  if(nameImgDom != "" ){
    typedef ImageMaskDomainCtrl<3> TImgContrl;
    typedef  HctccoModel<TImgContrl, 3> TTree;
    TImgContrl aDomCtr;
    TImgContrl::TPointI pM;
      if (!pInit->empty())
      {
          pM[0] = postInitV[0];
          pM[1] = postInitV[1];
          pM[2] = postInitV[2];
          aDomCtr = TImgContrl(nameImgDom, 128, pM, 100);
      }
      else
      {
          aDomCtr = TImgContrl(nameImgDom, 128, 100);
      }
    
    aDomCtr.myMinDistanceToBorder = minDistanceToBorder;
    TTree tree  (aPerf, nbTerm, aDomCtr, 1.0, extensionMode, anywhereMode);
    if (importXMLName != "") { 
      CoronaryTreeIO::readTreeFromXml(tree, importXMLName.c_str());
    }
    tree.my_gamma = gamma;
    tree.my_rootRadiusMax = rootRadiusMax;
    tree.my_termFlowCV = termFlowCV;
    tree.my_sideBranchBias = sideBranchBias;

    constructTreeMaskDomain(tree, verbose, importXMLName != "");
    
    CoronaryTreeIO::writeTreeToXml(tree, "tree_3D.xml");
    exportResultingMesh(tree, outputMeshName);
    #ifdef WITH_VISU3D_QGLVIEWER
    if (display3D) display3DTree(tree);
    #endif
    if (exportDatName != "") exportDat(tree, exportDatName);
    if (exportXMLName != "") CoronaryTreeIO::writeTreeToXml(tree, exportXMLName.c_str());
    if (postAnalyzeWall && exportXMLName != "") analyzeWallDistribution(tree, exportXMLName, nameImgDom);
  }
  else if (squaredImplDomain)
  {
      typedef SquareDomainCtrl<3> SqDomCtrl;
      typedef  HctccoModel<SqDomCtrl, 3> TTree;
      SqDomCtrl::TPoint pCenter (0,0,0);
      SqDomCtrl aCtr(1.0 ,pCenter);
      TTree tree  (aPerf, nbTerm, aCtr, 1.0, extensionMode, anywhereMode);
      if (importXMLName != "") { 
          CoronaryTreeIO::readTreeFromXml(tree, importXMLName.c_str());
      }
      tree.my_gamma = gamma;
    tree.my_rootRadiusMax = rootRadiusMax;
    tree.my_termFlowCV = termFlowCV;
    tree.my_sideBranchBias = sideBranchBias;
      
      constructTreeImplicitDomain(tree, outputMeshName,
                                    exportXMLName, importXMLName,
                                    exportDatName, verbose, display3D, postAnalyzeWall);
  }
  else
  {
    typedef CircularDomainCtrl<3> SphereDomCtrl;
    typedef  HctccoModel<SphereDomCtrl, 3> TTree;
    SphereDomCtrl::TPoint pCenter (0,0,0);
    SphereDomCtrl aCtr(1.0 ,pCenter);
    TTree tree  (aPerf, nbTerm, aCtr, 1.0, extensionMode, anywhereMode);
    if (importXMLName != "") { 
        CoronaryTreeIO::readTreeFromXml(tree, importXMLName.c_str());
    }
    tree.my_gamma = gamma;
    tree.my_rootRadiusMax = rootRadiusMax;
    tree.my_termFlowCV = termFlowCV;
    tree.my_sideBranchBias = sideBranchBias;
    
    constructTreeImplicitDomain(tree, outputMeshName,
                                  exportXMLName, importXMLName,
                                  exportDatName, verbose, display3D, postAnalyzeWall);
  }
  if (tools && exportXMLName != "") {
      runToolsOnXML(exportXMLName, nameImgDom);
  }
  return EXIT_SUCCESS;
}
