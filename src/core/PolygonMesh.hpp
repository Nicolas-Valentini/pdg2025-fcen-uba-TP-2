//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-05 16:34:26 taubin>
//------------------------------------------------------------------------
//
// PolygonMesh.hpp
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

#ifndef _POLYGONMESH_HPP_
#define _POLYGONMESH_HPP_

#include <vector>
#include "HalfEdges.hpp"

using namespace std;

class PolygonMesh : public HalfEdges {

public:

  // inherits from Edges
  //
  // void    reset(const int nV);
  // int     getNumberOfVertices()                     const;
  // int     getNumberOfEdges()                        const;
  // int     getEdge(const int iV0, const int iV1)     const;
  // int     getVertex0(const int iE)                  const;
  // int     getVertex1(const int iE)                  const;

  // inherits from HalfEdges
  //
  // int     getNumberOfCorners();
  // int     getFace(const int iC) const;
  // int     getSrc(const int iC) const;
  // int     getDst(const int iC) const;
  // int     getNext(const int iC) const;
  // int     getPrev(const int iC) const;
  // int     getTwin(const int iC) const;
  // int     getNumberOfEdgeHalfEdges(const int iE);
  // int     getEdgeHalfEdge(const int iE, const int j);

             PolygonMesh(const int nV, const vector<int>& coordIndex);

  // number of -1's in the coordIndex argument

     int     getNumberOfFaces()                        const;

  // number of faces incident to each edge; note that this is equal to
  // the number of half edges incident to each edge

     int     getNumberOfEdgeFaces(const int iE)        const;

  // if the arguments fall within their respective ranges, this method
  // returns the j-th face in the list of faces incident to the edge
  // iE; and it returns -1 if either argument is out of range

     int     getEdgeFace(const int iE, const int j)    const;

  // if the arguments fall within their respective ranges, this method
  // returns returns true if iF is found in the list of faces incident
  // to the edge iE; otherwise it returns false

     bool    isEdgeFace(const int iE, const int iF)    const;

  // edges are classified as boundary, regular, or singular depending
  // on the number of incident faces: 1=boundary, 2=regular, 3 or
  // more=singular

     bool    isBoundaryEdge(const int iE)               const;
     bool    isRegularEdge(const int iE)               const;
     bool    isSingularEdge(const int iE)              const;

  // a vertex is boundary if and only if it is the end of a boundary
  // edge

     bool    isBoundaryVertex(const int iV)            const;

  // a vertex is internal if and only if it is not a boundary edge

     bool    isInternalVertex(const int iV)            const;

  // a vertex is singular if the number of connected components in the
  // subgraph of the dual graph defined by the subset of faces
  // incident to the vertex is larger than 1; otherwise it is regular

     bool    isSingularVertex(const int iV)            const;

  // a way to determine which vertices are singular and which are
  // regular is to construct a partition of the corners of the mesh,
  // i.e. the elements of the coordIndex array, including the -1
  // separators. For each regular edge there are exactly two
  // half-edges incident to the edge.

  // If the two incident faces are consistently oriented, the
  // neighborhood of the edge will look like this

  //
  /* \                  / */
  /*  \ iC01 <-- iC00  /  */
  /*   X ------------ X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // If the two incident faces are not consistently oriented, the
  // neighborhood of the edge will look like this

  //
  /* \                  / */
  /*  \ iC00 --> iC01  /  */
  /*   X ------------ X   */
  /*  / iC10 --> iC11  \  */
  /* /                  \ */

  // In either case the two pairs of corners which are opposite to
  // each other accross the edge jshould be joined. Once all of these
  // join operations have been applied, the corners are partitioned
  // into parts. Some of these parts are singletons corresponding to
  // the face separators. Any other part is composed of corners which
  // point to the same vertex. Counting how many of these parts point
  // to the same vertex determines which vertices are regular and
  // which are singular. If only one of these parts point to a given
  // vertex, then the vertex is regular. If two or more of these parts
  // point to a given vertex, then the vertex is singular. 

  // the polygon mesh is regular if and only if it does not have any
  // singular edges and it does not have any singular vertices

     bool    isRegular()                               const;

  // the polygon mesh has boundary if and only if it has at least 1
  // boundary edge

     bool    hasBoundary()                             const;
  
private:

  // consider these private variables a suggestion
  // feel free to decide how to implement this class

  vector<int>      _nPartsVertex;
  vector<bool> _isBoundaryVertex;
  bool _invalidFace(const int iF) const;
  bool _hasAllRegularEdges() const;
  bool _hasAllRegularVertex() const;


  
};

#endif /* _POLYGONMESH_HPP_ */
