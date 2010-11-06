
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr, Greg Humphreys and Bohumir Zamecnik.

    This file is part of pbrt and BokehLab.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_THINLENSTILTSHIFT_H
#define PBRT_CAMERAS_THINLENSTILTSHIFT_H

// cameras/thinlenstiltshift.h*
#include "pbrt.h"
#include "camera.h"
#include "film.h"

// ThinLensTiltShiftCamera Declarations
class ThinLensTiltShiftCamera : public ProjectiveCamera {
public:
    // ThinLensTiltShiftCamera Public Methods
    /**
     * // TODO: in future use a 2D vector for horizonal/vertical tilt and shift
     * cam2world
     * screenWindow
     * shutterOpen
     * shutterClose
     * lensRadius
     * focalLength
     * imageDistance
     * filmtiltx - in degrees
     * filmshifty
     * fov
     * film
     */
    ThinLensTiltShiftCamera(
        const AnimatedTransform &cam2world,
        const float screenWindow[4],
        float shutterOpen, float shutterClose,
        float lensRadius,
        float focalLength,
        float imageDistance,
        float filmtiltx,
        const float filmshift[2],
        float fov,
        Film *film
    );

    float GenerateRay(
        const CameraSample &sample,
        Ray *
    ) const;

    //float GenerateRayDifferential(
    //    const CameraSample &sample,
    //    RayDifferential *ray
    //) const;
private:
    // ThinLensTiltShiftCamera Private Data
    Vector dxCamera, dyCamera;
    float filmTiltX;
    Vector filmShift;
    // focal length - the internal property of the lens
    float focalLength;
    // Note: in contrast to it the "focal distance" is the distance of the
    // focused plane (with no tilt) to the center of the lens.
    // Iit depends on the current image distance.

    // internal data for tilt-shift computation
    Ray scheimpflugLine;
    Ray hingeLine;
    boolean tiltEnabled;
    Transform filmTiltRotate;

    static float computeFocalDistance(float focalLength, float imageDistance);
    Point getCameraSample(const CameraSample &sample) const;
    Point getFilmSample(const Point &cameraSample) const;
    void modifyRayForDof(Ray &ray, const CameraSample &sample) const;
};


ThinLensTiltShiftCamera *CreateThinLensTiltShiftCamera(
    const ParamSet &params,
    const AnimatedTransform &cam2world,
    Film *film
);

#endif // PBRT_CAMERAS_THINLENSTILTSHIFT_H
