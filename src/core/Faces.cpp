//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-04 22:09:56 gtaubin>
//------------------------------------------------------------------------
//
// Faces.cpp
//
// Written by: <Your Name>
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Brown University nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL GABRIEL TAUBIN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <math.h>
#include "Faces.hpp"

void Faces::createInvalidFaces_(){
    //If the coordIndex in the constructor is invalid, Faces shouldnt work
    m_CoordIndex = vector<int>();
    m_FacesIndex = vector<int>();
    m_nVerts = 0;
}
  
Faces::Faces(const int nV, const vector<int>& coordIndex) {
    int currentFaceSize = 0;
    for (int i = 0; i < coordIndex.size(); ++i) {
        int currentCorner = coordIndex[i];
        if(currentCorner >= nV){
            createInvalidFaces_();
            break;
        }
        else{
            if(currentCorner == -1){
                //Its impossible to have a face with less than 2 vertex
                if( currentFaceSize < 2){
                    createInvalidFaces_();
                    break;
                }
                else{
                    m_CoordIndex.push_back(-1);
                    //m_FacesIndex has the index of the end of the face i
                    m_FacesIndex.push_back(i);
                    currentFaceSize = 0;
                }
            }
            else{
                m_CoordIndex.push_back(currentCorner);
                currentFaceSize++;
            }

        }
    }
    m_nVerts = nV;


}

int Faces::getNumberOfVertices() const {
  return m_nVerts;
}

int Faces::getNumberOfFaces() const {
  return m_FacesIndex.size();
}

int Faces::getNumberOfCorners() const {
  return m_CoordIndex.size();
}

bool Faces::outOfRangeFace_(const int iF) const{
    return iF < 0 || iF >= getNumberOfFaces();
}

int Faces::getFaceSize(const int iF) const {
    if (outOfRangeFace_(iF))
        return 0;
    else{
        if(iF == 0)
            return m_FacesIndex[0];
        return m_FacesIndex[iF] - m_FacesIndex[iF-1]-1;
    }
}

int Faces::getFaceFirstCorner(const int iF) const {
    if (outOfRangeFace_(iF))
        return -1;
    else{
        int faceSize = getFaceSize(iF);
        int faceIndex = m_FacesIndex[iF];
        return m_CoordIndex[faceIndex-faceSize];
    }
}

int Faces::getFaceVertex(const int iF, const int j) const {
    int faceSize = getFaceSize(iF);
    if (outOfRangeFace_(iF) || j<0 || j>=faceSize)
        return -1;
    else{
        int faceIndex = m_FacesIndex[iF];
        return m_CoordIndex[faceIndex-faceSize+j];
    }
}
bool Faces::invalidCorner_(const int iC) const{
    return iC < 0 || iC >= getNumberOfCorners() || m_CoordIndex[iC] == -1;
}

int Faces::getCornerFace(const int iC) const {
    if(invalidCorner_(iC))
        return -1;
    else{
        for (int i = 0; i < m_FacesIndex.size(); ++i) {
            if(i == 0 && m_FacesIndex[0] > iC)
                return 0;
            else{
                if(m_FacesIndex[i-1] < iC && m_FacesIndex[i] > iC)
                    return i;
            }

        }
    }
    return -1;

}

int Faces::getNextCorner(const int iC) const {
    if(invalidCorner_(iC))
        return -1;
    else{
        if(m_CoordIndex[iC+1]==-1){
            int face = getCornerFace(iC);
            return getFaceFirstCorner(face);
        }
        else{
            return m_CoordIndex[iC];
        }
    }
}

