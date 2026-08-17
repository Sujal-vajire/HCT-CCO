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

#if defined(XMLHELPERS_RECURSES)
#error Recursive header files inclusion detected in Xmlhelpers.h
#else // defined(XMLHELPERS_RECURSES)
/** Prevents recursive inclusion of headers. */
#define XMLHELPERS_RECURSES

#if !defined XMLHELPERS_h
/** Prevents repeated inclusion of headers. */
#define XMLHELPERS_h
#include <string.h>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include "HctccoModel.h"
using namespace std;

class NodeTable {
public:
    enum Possibilite {TERM, ROOT, BIF};
    /*
     //constants used for the type field of a node
     static int TERM;
     static int ROOT;
     static int BIF;
     static int FIELDS;
     
     NodeTable();
     */
};

namespace CoronaryTreeIO {
/**
 * prints a node into XML/GXL format from the NodeTable
 */
template <typename TPoint>
void subPrint_node(int nodeType, TPoint pos,
                   int idNode, ofstream &os){
    os<<"  <node id=\"n"<<idNode<<"\">"<<endl;
    os<<"    <attr name=\" nodeType\">"<<endl;
    if(nodeType == NodeTable::ROOT){
        os<<"      <string> root node </string>"<<endl;
    } else if(nodeType == NodeTable::TERM){
        os<<"      <string> terminal node </string>"<<endl;
    } else if(nodeType == NodeTable::BIF){
        os<<"      <string> bifurication </string>"<<endl;
    } else {
        os<<"      <string> unknown type </string>"<<endl;
    }
    os<<"    </attr>"<<endl;
    
    os<<"    <attr name=\" position\">"<<endl;
    os<<"      <tup>"<<endl;
    
    for (auto i=0; i <pos.dimension; i++){
        os<<"        <float>"<<pos[i]<<"</float>"<<endl;
    }
    
    os<<"      </tup>"<<endl;
    os<<"    </attr>"<<endl;
    os<<"  </node>"<<endl;
}

/**
 * prints an edge into XML/GXL format from a node table
 */

void subPrint_edge(int idSeg, int idSegPar,
                   double flow, double radius,double resist,
                   ofstream &os){
    
    if(idSeg != 0){
        // Sanitize numeric values to avoid NaN/Inf in XML
        auto safeVal = [](double v) {
            if (!std::isfinite(v)) return 0.0;
            return v;
        };
        flow = safeVal(flow);
        radius = safeVal(radius);
        resist = safeVal(resist);
        os<<"  <edge id=\"e"<<idSeg<<"\" to=\"n"<<idSeg<<"\" from=\"n"<<idSegPar<<"\">"<<endl;
        os<<"    <attr name=\" flow\">"<<endl;
        os<<"      <float>"<<flow<<"</float>"<<endl;
        os<<"    </attr>"<<endl;
        os<<"    <attr name=\" resistance\">"<<endl;
        os<<"      <float>"<<resist<<"</float>"<<endl;
        os<<"    </attr>"<<endl;
        
        os<<"    <attr name=\" radius\">"<<endl;
        os<<"      <float>"<<radius<<"</float>"<<endl;
        os<<"    </attr>"<<endl;
        
        os<<"  </edge>"<<endl;
    }
}


template<typename TTree>
inline
void writeTreeToXml(const TTree &tree, const char * filePath) {
    ofstream output;
    //writing the tree structure as GXL to the filePath specified
    output.open(filePath);
    output<<"<gxl><graph id=\""<<filePath<<"\" edgeids=\" true\" edgemode=\" directed\" hypergraph=\" false\">"<<endl;
    output<<"<info_graph>"<< endl;
    output<<"    <attr name=\" pPerf\">"<<endl;
    output<<"      <float>"<<tree.my_pPerf<<"</float>"<<endl;
    output<<"    </attr>"<<endl;
    output<<"    <attr name=\" pTerm\">"<<endl;
    output<<"      <float>"<<tree.my_pTerm<<"</float>"<<endl;
    output<<"    </attr>"<<endl;
    output<<"</info_graph>"<< endl;
    
    //writing tree's nodes
    for(auto s : tree.myVectSegments) {
        // test if the segment is the root or its parent
        if (tree.myVectParent[s.myIndex]==0) //root node
            subPrint_node(NodeTable::ROOT, s.myCoordinate, s.myIndex, output);
        else {
            if(std::find(begin(tree.myVectTerminals), end(tree.myVectTerminals), s.myIndex) != end(tree.myVectTerminals)) { //terminal node
                subPrint_node(NodeTable::TERM, s.myCoordinate, s.myIndex, output);
            }
            else { // bif node
                subPrint_node(NodeTable::BIF, s.myCoordinate, s.myIndex, output);
            }
        }
    }
    //writing tree's edges
    for(auto s : tree.myVectSegments) {
        subPrint_edge(s.myIndex,tree.myVectParent[s.myIndex], s.myFlow, s.myRadius, s.myResistance, output);
    }
    
    output<<"</graph></gxl>"<<endl;
    output.close();
}

template<typename TTree>
inline
void readTreeFromXml(TTree &tree, const char * filePath) {
    ifstream file(filePath);
    if(!file.is_open()) {
        throw std::runtime_error("Could not open XML file");
    }
    // helper variables
    string line;
    vector<string> lines;
    typename TTree::Segment segment;
    vector<typename TTree::Segment> segments;
    unsigned int currentNodeID = 0;
    NodeTable::Possibilite type; 
    float x, y, z;
    int from, to;
    tree.myVectSegments.clear();

    cout << "Reading XML File..." << endl;
    for (; getline(file, line); )
        lines.push_back(line);
    file.close();

    for (unsigned int i = 0; i < lines.size(); i++) {
        DGtal::trace.progressBar(i, lines.size());
        line = lines[i];
        if (line.find("<node") != std::string::npos) {
            if (sscanf(line.c_str(), " <node id=\"n%d\">", &currentNodeID) == 1) {
                // init empty node
                segment = {};
                segment.myIndex = currentNodeID;
            }
        } else if (line.find("<attr name=\" nodeType\">") != std::string::npos) {
            i++;
            line = lines[i];
            if (line.find("<string> terminal node </string>") != std::string::npos) {
                // if term type then add index to vector
                type = NodeTable::Possibilite::TERM;
                tree.myVectTerminals.push_back(segment.myIndex);
            }
            // All other attributes besides nodeType
        } else if (line.find("<attr name=\"") != std::string::npos) {
            char val[64];
            if (sscanf(line.c_str(), " <attr name=\" %63[^\"]\">", &val) == 1) {
                string value(val);
                if (value == "pPerf") {
                    i++;
                    line = lines[i];
                    float pPerf;
                    if (sscanf(line.c_str(), " <float>%f</float>", &pPerf) == 1) {
                        tree.my_pPerf = pPerf;
                    }
                } else if (value == "pTerm") {
                    i++;
                    line = lines[i];
                    float pTerm;
                    if (sscanf(line.c_str(), " <float>%f</float>", &pTerm) == 1) {
                        tree.my_pTerm = pTerm;
                    }
                } else if (value == "position") {
                    i += 2;
                    line = lines[i];
                    if (sscanf(line.c_str(), " <float>%f</float>", &x) == 1) {
                        segment.myCoordinate[0] = x;
                    }
                    i++;
                    line = lines[i];
                    if (sscanf(line.c_str(), " <float>%f</float>", &y) == 1) {
                        segment.myCoordinate[1] = y;
                    }
                    i++;
                    line = lines[i];
                    if (sscanf(line.c_str(), " <float>%f</float>", &z) == 1) {
                        segment.myCoordinate[2] = z;
                    }
                } else if (value == "flow") {
                    i++;
                    line = lines[i];
                    float flow;
                    if (sscanf(line.c_str(), " <float>%f</float>", &flow) == 1) {
                        for(typename TTree::Segment& s : segments) { 
                            if (s.myIndex == to) {
                                s.myFlow = flow;
                                break;
                            }
                        }
                    }
                } else if (value == "resistance") {
                    i++;
                    line = lines[i];
                    float resistance;
                    if (sscanf(line.c_str(), " <float>%f</float>", &resistance) == 1) {
                        for(typename TTree::Segment& s : segments) { 
                            if (s.myIndex == to) {
                                s.myResistance = resistance;
                                break;
                            }
                        }
                    }
                } else if (value == "radius") {
                    i++;
                    line = lines[i];
                    float radius;
                    if (sscanf(line.c_str(), " <float>%f</float>", &radius) == 1) {
                        for(typename TTree::Segment& s : segments) { 
                            if (s.myIndex == to) {
                                s.myRadius = radius;
                                break;
                            }
                        }
                    }
                }
            }
        } else if (line.find("</node>") != std::string::npos) {
            if (segment.myIndex > 1) {
                // extend the vectors
                tree.myVectParent.push_back(0);
                tree.myVectChildren.push_back(std::pair<unsigned int, unsigned int>(0, 0));
            }
            segments.push_back(segment);
        } else if (line.find("<edge") != std::string::npos) {
            // detect from and to in edge, and update parents accordingly
            if (sscanf(line.c_str(), "  <edge id=\"e%*d\" to=\"n%d\" from=\"n%d\">", &to, &from) == 2) {
                tree.myVectParent[to] = from;
                // if there is no first child, add there
                // otherwise add to the second
                if (tree.myVectChildren[from].first == 0)
                    tree.myVectChildren[from].first = to;
                else
                    tree.myVectChildren[from].second = to;
            }
        }
    }
    for(typename TTree::Segment s : segments) {
        // extend the myVectSegments to the proper size
        tree.myVectSegments.push_back(s);
    }
    // set the first parent to 1, always the case
    tree.myVectParent[0] = 1;

    // ------------------------------------------------------------------
    // Reconcile internal tree state with the imported topology.
    //
    // Why: the constructor built the tree's internal accounting (myKTerm,
    // per-segment kTerm/flow/resistance/beta, length factor) for a fresh
    // 1-segment tree. readTreeFromXml then cleared myVectSegments and
    // repopulated them from the XML, but did NOT refresh that internal
    // state, so subsequent growth ran against stale global accounting
    // (progress=0-or-wrong, length factor for kTerm=1, mixed resistances).
    //
    // Steps:
    //   1. Reset per-segment kTerm/flow, then propagate each terminal up
    //      to populate ancestor kTerm/flow consistently.
    //   2. Set global myKTerm to the number of imported terminals.
    //   3. Recompute the length factor from (myKTerm, myRsupp, my_rPerf).
    //   4. Full root-down resistance recompute and radii from the root.
    // ------------------------------------------------------------------
    for (auto &seg : tree.myVectSegments) {
        seg.myKTerm = 0;
        seg.myFlow  = 0.0;
    }
    // Imported terminals use uniform (unit-weight) perfusion flow; their real
    // territories are unknown, so we anchor them to q_term and let subsequent
    // grown terminals carry any non-uniform weighting.
    for (unsigned int termIdx : tree.myVectTerminals) {
        tree.updateFlow(termIdx, tree.my_qTerm);
    }
    tree.myKTerm = static_cast<unsigned int>(tree.myVectTerminals.size());
    tree.updateLengthFactor();
    tree.updateResistanceFromRoot(1);
    tree.updateRootRadius();
    tree.myLastUpdatedSegment = 1;

    cout << endl << "XML reading finished" << endl;
    writeTreeToXml(tree, "test.xml");
}

}; // end namespace CoronaryTreeIO

#endif // !defined XMLHELPERS_h

#undef XMLHELPERS_RECURSES
#endif // else defined(XMLHELPERS_RECURSES)
