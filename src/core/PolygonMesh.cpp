//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 16:34:29 taubin>
//------------------------------------------------------------------------
//
// PolygonMesh.cpp
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//     * Redistributions of source code must retain the above
//       copyright notice, this list of conditions and the following
//       disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials
//       provided with the distribution.
//     * Neither the name of the Brown University nor the names of its
//       contributors may be used to endorse or promote products
//       derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GABRIEL
// TAUBIN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include <iostream>
#include "PolygonMesh.hpp"
#include "Partition.hpp"

PolygonMesh::PolygonMesh(const int nVertices, const vector<int>& coordIndex):
  HalfEdges(nVertices,coordIndex),
  _nPartsVertex(),
  _isBoundaryVertex()
{
  int nV = getNumberOfVertices();
  int nE = getNumberOfEdges(); // Edges method
  // int nF = getNumberOfFaces();
  int nC = getNumberOfCorners();

  // 1) classify the vertices as boundary or internal
  int iV;
  for(iV=0;iV<nV;iV++)
    _isBoundaryVertex.push_back(false);
  // - for edge boundary iE label its two end vertices as boundary

  for(int e=0;e<nE;e++){
      int k = getNumberOfEdgeHalfEdges(e);
      if(k==1){
        int v0 = getVertex0(e);
        int v1 = getVertex1(e);
        _isBoundaryVertex[v0] = true;
        _isBoundaryVertex[v1] = true;
      }
  }
  
  // 2) create a partition of the corners in the stack
  Partition partition(nC);
  // 3) for each regular edge
  //    - get the two half edges incident to the edge
  //    - join the two pairs of corresponding corners accross the edge
  //    - you need to take into account the relative orientation of
  //      the two incident half edges

  for(int e = 0; e < nE; ++e){
      int k = getNumberOfEdgeHalfEdges(e);
      if(k == 2){
          int a = getEdgeHalfEdge(e, 0);
          int b = getEdgeHalfEdge(e, 1);

          int a_next = getNext(a);
          int b_next = getNext(b);

          partition.join(a_next, b);
          partition.join(b_next, a);
      }
  }

  // consistently oriented
  /* \                  / */
  /*  \ iC01 <-- iC00  /  */
  /*   X ---- iE ---- X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // oposite orientation
  /* \                  / */
  /*  \ iC01 --> iC00  /  */
  /*   X ---- iE ---- X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // a decision has to be made about inconsistently oriented faces
  // incident to the same edge, as well as how to deal with singular
  // edges; for the moment let's assume that the mesh does not have
  // singular edges, and that pairs of corners corresponding to the
  // same vertex across inconsistently oriented faces will be joined

  // note that the partition will end up with the corner separators as
  // singletons, but it doesn't matter for the last step, and
  // the partition will be deleteted upon return
  
  // 4) count number of parts per vertex
  //    - initialize _nPartsVertex array to 0's
  //    - for each corner iC which is a representative of its subset, 
  //    - get the corresponding vertex index iV and increment _nPartsVertex[iV]
  //    - note that all the corners in each subset share a common
  //      vertex index, but multiple subsets may correspond to the
  //      same vertex index, indicating that the vertex is singular

  _nPartsVertex.assign(nV, 0);
  vector<bool> seen(nC, false);

  for(int ic = 0; ic < nC; ++ic){
      if(_coordIndex[ic] < 0) continue;
      int rep = partition.find(ic);
      if(seen[rep]) continue;
      seen[rep] = true;
      int v = getSrc(ic);
      _nPartsVertex[v] += 1;
  }
}

int PolygonMesh::getNumberOfFaces() const {
    int count = 0;
    for(int i = 0; i < _coordIndex.size(); ++i) if(_coordIndex[i] == -1) ++count;
    return count;
}

int PolygonMesh::getNumberOfEdgeFaces(const int iE) const {
  return getNumberOfEdgeHalfEdges(iE);
}

int PolygonMesh::getEdgeFace(const int iE, const int j) const {
    if(iE < 0 || iE >= getNumberOfEdges()) return -1;
    int n = getNumberOfEdgeHalfEdges(iE);
    if(j < 0 || j >= n) return -1;
    int c = getEdgeHalfEdge(iE, j);
    if(c < 0) return -1;
    return getFace(c);
}

bool PolygonMesh::isEdgeFace(const int iE, const int iF) const {
    if(iE < 0 || iE >= getNumberOfEdges()) return false;
    int n = getNumberOfEdgeHalfEdges(iE);
    for(int j=0;j<n;j++){
        int f = getEdgeFace(iE,j);
        if(f == iF) return true;
    }
    return false;
}

// classification of edges

bool PolygonMesh::isBoundaryEdge(const int iE) const {
    if(iE < 0 || iE >= getNumberOfEdges()) return false;
    return getNumberOfEdgeHalfEdges(iE) == 1;
}

bool PolygonMesh::isRegularEdge(const int iE) const {
    if(iE < 0 || iE >= getNumberOfEdges()) return false;
    return getNumberOfEdgeHalfEdges(iE) == 2;
}

bool PolygonMesh::isSingularEdge(const int iE) const {
    if(iE < 0 || iE >= getNumberOfEdges()) return false;
    return getNumberOfEdgeHalfEdges(iE) >= 3;
}

// classification of vertices

bool PolygonMesh::isBoundaryVertex(const int iV) const {
  int nV = getNumberOfVertices();
  return (0<=iV && iV<nV)?_isBoundaryVertex[iV]:false;
}

bool PolygonMesh::isSingularVertex(const int iV) const {
  int nV = getNumberOfVertices();
  return (0<=iV && iV<nV && _nPartsVertex[iV]>1);
}

// properties of the whole mesh

bool PolygonMesh::isRegular() const {
    int nE = getNumberOfEdges();
    for(int e=0;e<nE;e++) if(isSingularEdge(e)) return false;

    int nV = getNumberOfVertices();
    for(int v=0; v<nV; ++v) if(isSingularVertex(v)) return false;

    return true;
}

bool PolygonMesh::hasBoundary() const {
    int nE = getNumberOfEdges();
    for(int e=0;e<nE;e++) if(isBoundaryEdge(e)) return true;
    return false;
}
