//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 16:34:28 taubin>
//------------------------------------------------------------------------
//
// HalfEdges.cpp (Assignment 2)
//
// Written by: <Your Name>
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

#include "io/StrException.hpp"
#include <math.h>
#include "HalfEdges.hpp"
#include "Graph.hpp"

// 1) all half edges corresponding to regular mesh edges are made twins
// 2) all the other edges are made boundary half edges (twin==-1)

HalfEdges::HalfEdges(const int nVertices, const vector<int>&  coordIndex):
  Edges(nVertices), // a graph with no edges is created here
  _coordIndex(coordIndex),
  _twin(),
  _face(),
  _firstCornerEdge(),
  _cornerEdge()
{

  // - both the _twin array and the _face array should end up being of
  //   the same size as the _coordIndex array
  // - for each corner index iC contained in face iF
  // - the half edge src is iC, and the half edge dst is iC+1 if iC is
  //   not the last corner of the face; otherwise is iC0
  //
  //   if _coordIndex[iC]>=0 then
  //     _face[iC] should be equal to iF
  //     if 
  //     _twin[iC] should be equal to the corner index of the twin half edge 
  //   if _coordIndex[iC]<0 then
  //   _face[
  int nV = nVertices;
  int nC = static_cast<int>(_coordIndex.size()); // number of corners

  // 0) just to be safe, verify that for each corner iC that
  //    -1<=iV && iV<nV, where iV=coordIndex[iC]
  //
  //    if you find a violation, you can increment and the variable
  //    nV, and then use the method Edges::_reset() to adjust the
  //    number of vertices of the graph, if necessary; or you can
  //    abort throwing an StrException
  
  // 1) create an empty vector<int> to count the number of incident
  //    faces per edge; size is not known at this point because the
  //    edges have not been created yet
  vector<int> nFacesEdge;

  // 2) insert all the edges in the graph; at the same time initialize
  //    the _twin array so that all the half edges are boundary, count
  //    the number of incident faces per edge, fill the _face
  //    array, and count the number of faces incident to each edge
  int iV0,iV1,iF=0,iE,iC,iC0=0,iC1;
  for(iC=0;iC<nC;iC++) {

      if(_coordIndex[iC]>=nV || _coordIndex[iC]<-1 ){
          throw new StrException("Invalid coordIndex");
      }

    // face iF comprises corners iC0<=iC<iC1
    // - each corner in this range corresponds to one half edge
    // - find the next corner within the face
    // - get the two vertex indices and insert an edge in the graph if
    //   not already there
    if(_coordIndex[iC] == -1){
        iC0 = iC+1;
        _face.push_back(-1);
        iF++;
        continue;
    }


    if(_coordIndex[iC+1]==-1)
        iC1 = iC0;
    else
        iC1 = iC+1;


    iV0 = _coordIndex[iC];
    iV1 = _coordIndex[iC1];
    iE = _insertEdge(iV0,iV1);
    if(iE >= nFacesEdge.size())
        nFacesEdge.push_back(1);
    else
        nFacesEdge[iE]++;

    // increment variables to continue processing next face
    _face.push_back(iF);

  }
  _nF = iF;

  int nE = getNumberOfEdges();
  
  // 3) create an array to hold the first twin corner for each edge
  vector<int> twinCorner;
  twinCorner.assign(nE, -1);
  // - the size of this array should be equal to the number of edges
  // - initialize it with -1's

  // 4) fill the _twin array
  // - visit all the half edges using a loop similar to the one used in step 2)
  // - for each half edge iC, get the src and dst vertex indices, and
  //   from them the index iE of the corresponding edge
  // - if twinCorner[iE]<1 save iC in twinCorner[iE]
  //   otherwise save the value stored in twinCorner[iE] in _twin[iC]
  //   and iC in _twin[_twin[iC]]

  _twin.assign(nC, -1);
  int i0 = 0;
  int i1;
  int currentFaceSize = 0;

  for (int i = 0; i < nC; ++i) {

      if(_coordIndex[i] == -1){
          i0 = i+1;
          _twin[i] = currentFaceSize;
          currentFaceSize = 0;
          continue;
      }
      currentFaceSize++;

      int src = _coordIndex[i];
      if(_coordIndex[i+1]==-1)
          i1 = i0;
      else
          i1 = i+1;
      int dst = _coordIndex[i1];

      int edge = getEdge(src,dst);
      if(twinCorner[edge]<0){
          twinCorner[edge] = i;
      }
      else{
          int other = twinCorner[edge];
          _twin[i] = other;
          _twin[other] = i;
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

  // a decision has to be made about inconsistently oriented half
  // edges incident to the same edge, as well as how to deal with
  // singular edges; for the moment let's assume that the mesh does
  // not have singular edges, but inconsistently oriented half edges
  // incident to the same edge are made twins (i.e. we do not have to
  // check for orientation here); later on we may want to modify this
  // class to have the option to do one thing or the other, and
  // methods to indicate which case we have.

  // get everything up to here implemented, debugged, and commited
  // before continuing

  // 5) initialize the array of arrays representing the half-edge to
  //    edge incident relationships
  //    _firstCornerEdge, and _cornerEdge
  //    - the size of _firstCornerEdge should be equal to nE+1
  //    - the size of _cornerEdge should be equal to the number of valid corners
  //      nC-nF
  //    - set boundaries
  //      _firstCornerEdge[0]=0
  //      _firstCornerEdge[iE+1] = _firstCornerEdge[iE]+nFacesEdge[iE] (1<=iE<nE)
  _firstCornerEdge.assign(nE+1, 0);
  for(int e = 0; e < nE; ++e) {
      _firstCornerEdge[e+1] = _firstCornerEdge[e] + nFacesEdge[e];
  }
  _cornerEdge.assign(nC-_nF,-1);
  int ic0=0,ic1;
  for (int ic = 0; ic < _coordIndex.size(); ++ic) {
      if (_coordIndex[ic] == -1) {
          ic0 = ic + 1;
          continue;
      }

      int src = _coordIndex[ic];
      if(_coordIndex[ic+1]==-1)
          ic1 = ic0;
      else
          ic1 = ic+1;
      int dst = _coordIndex[ic1];

      int edge = getEdge(src,dst);

      int halfEdgeToEdgeStart = _firstCornerEdge[edge];

      while(_cornerEdge[halfEdgeToEdgeStart]!=-1)halfEdgeToEdgeStart++;
      _cornerEdge[halfEdgeToEdgeStart] = ic;
  }

  // 6) fill the array of arrays - the indices of corners incident to
  //    edge iE (1 if boundary, 2 if regular, >2 if singular) should
  //    be stored consecutively in _cornerEdge starting at the
  //    location _firstCornerEdge[iE]
}

int HalfEdges::getNumberOfCorners() const {
  return static_cast<int>(_coordIndex.size());
}

// in all subsequent methods check that the arguments are valid, and
// return -1 if any argument is out of range

// half-edge method srcVertex()
bool HalfEdges::_invalidCorner(const int iC) const{
    return iC<0 || iC >= getNumberOfCorners() || _coordIndex[iC] == -1;
}

int HalfEdges::getFace(const int iC) const {
    if(_invalidCorner(iC))
        return -1;

  return _face[iC];
}

// half-edge method srcVertex()
int HalfEdges::getSrc(const int iC) const {
    if(_invalidCorner(iC))
        return -1;
  return _coordIndex[iC];
}

// half-edge method dstVertex()
int HalfEdges::getDst(const int iC) const {
    if(_invalidCorner(iC))
        return -1;
    if(_coordIndex[iC+1]==-1){
        int faceSize = _twin[iC+1];
        return _coordIndex[iC-faceSize+1];
    }
    return _coordIndex[iC+1];
}

// half-edge method next()
int HalfEdges::getNext(const int iC) const {
  // if iC is the last corner of its face, use the face size
  // stored in _twin[iC+1] to locate the first corner of the face
  if(_invalidCorner(iC))
      return -1;
  if(_coordIndex[iC+1]==-1){
      int faceSize = _twin[iC+1];
      return iC-faceSize+1;
  }
  return iC+1;
}

// half-edge method prev()
int HalfEdges::getPrev(const int iC) const {
  
  // if iC is the first corner of its face, since the face size is
  // stored at the end of the face in the _twin array, you will have
  // to search forward for the last corner of the face; you can use
  // the fact that all the faces have at least 3 corners to start the
  // search for the face separator at iC+3
  
  if(_invalidCorner(iC))
      return -1;

  if(iC==0 || _coordIndex[iC-1]==-1){
      int i = iC+3;
      while(_coordIndex[i]!=-1)i++;
      return i-1;
  }
  return iC-1;
}

int HalfEdges::getTwin(const int iC) const {
    if(_invalidCorner(iC))
        return -1;
    return _twin[iC];
}

// represent the half edge as an array of lists, with one list
// associated with each edge

bool HalfEdges::_invalidEdge(const int iE) const{
    return iE<0 || iE >= getNumberOfEdges();
}

int HalfEdges::getNumberOfEdgeHalfEdges(const int iE) const {
    if(_invalidEdge(iE)) return 0;
    return _firstCornerEdge[iE+1] - _firstCornerEdge[iE];
}

int HalfEdges::getEdgeHalfEdge(const int iE, const int j) const {
    if(_invalidEdge(iE)) return -1;
    int n = getNumberOfEdgeHalfEdges(iE);
    if(j < 0 || j >= n) return -1;
    return _cornerEdge[_firstCornerEdge[iE] + j];
}
