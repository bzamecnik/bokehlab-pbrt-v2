
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


// cameras/thinlenstiltshift.cpp*
#include "stdafx.h"
#include "cameras/thinlenstiltshift.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"

// ThinLensTiltShiftCamera Method Definitions
ThinLensTiltShiftCamera:: ThinLensTiltShiftCamera(
        const AnimatedTransform &cam2world,
        const float screenWindow[4],
        float shutterOpen, float shutterClose,
        float lensRadius,
        float focalLength_,
        float imageDistance,
        float filmtiltx, float filmshifty,
        float fov,
        Film *film
) : ProjectiveCamera(cam2world, Perspective(fov, 1e-2f, 1000.f),
        screenWindow, shutterOpen, shutterClose, lensRadius,
        computeFocalDistance(focalLength_, imageDistance),
        film),
    filmTiltX(filmtiltx), filmShift(0., filmshifty, -imageDistance),
    focalLength(focalLength_),
    scheimpflugLineX(0.), hingeLineX(0.),
    filmTiltRotate(
        Translate(Vector(0, 0, -imageDistance)) *
        RotateY(filmtiltx) *
        Translate(Vector(0, 0, imageDistance)))
{
    // Compute differential changes in origin for perspective camera rays
    dxCamera = RasterToCamera(Point(1,0,0)) - RasterToCamera(Point(0,0,0));
    dyCamera = RasterToCamera(Point(0,1,0)) - RasterToCamera(Point(0,0,0));

    tiltEnabled = (filmtiltx > 0.) || (filmtiltx < 0.);

    if (tiltEnabled) {
        float sine_tilt_angle = sinf(Radians(filmtiltx));
        scheimpflugLineX = imageDistance / sine_tilt_angle;
        hingeLineX = focalLength / sine_tilt_angle;
    }
}


float ThinLensTiltShiftCamera::GenerateRay(
    const CameraSample &sample,
    Ray *ray
) const {
    // Get the ray going through the lens center
    Point filmSample = getFilmSample(getCameraSample(sample));
    *ray = Ray(Point(0, 0, 0), Vector(filmSample), 0.f, INFINITY);
    
    // Modify ray for depth of field
    modifyRayForDof(*ray, sample);
    
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    CameraToWorld(*ray, ray);
    return 1.f;
}


float ThinLensTiltShiftCamera::GenerateRayDifferential(
    const CameraSample &sample,
    RayDifferential *ray
) const {
    // Get the ray going through the lens center
    Point cameraSample = getCameraSample(sample);
    Point filmSample = getFilmSample(cameraSample);
    *ray = RayDifferential(Point(0,0,0), Normalize(Vector(filmSample)), 0.f, INFINITY);
    
    // Modify ray for depth of field
    modifyRayForDof(*ray, sample);

    // Compute offset rays for _ThinLensTiltShiftCamera_ ray differentials
    ray->rxOrigin = ray->ryOrigin = ray->o;
    ray->rxDirection = Normalize(Vector(cameraSample) + dxCamera);
    ray->ryDirection = Normalize(Vector(cameraSample) + dyCamera);
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    CameraToWorld(*ray, ray);
    ray->hasDifferentials = true;
    return 1.f;
}

float ThinLensTiltShiftCamera::computeFocalDistance(
    float focalLength,
    float imageDistance
) {
    return 1 / ((1 / focalLength) - (1 / imageDistance));
}

/**
 * Convert sample from Raster space to Camera space.
 * cameraSample->z should be mapped from 0 to the near plane.
 */
inline Point ThinLensTiltShiftCamera::getCameraSample(
    const CameraSample &sample
) const {
    return RasterToCamera(Point(sample.imageX, sample.imageY, 0));
}

/**
 * Get film sample in Camera space.
 */
inline Point ThinLensTiltShiftCamera::getFilmSample(
    const Point &cameraSample
) const {
    // project the point from the near plane to the untilted film
    float filmT = filmShift.z / cameraSample.z;
    Point filmSample = Point(cameraSample.x * filmT, cameraSample.y * filmT, filmShift.z);
    // tilt the film
    return filmTiltRotate(filmSample);
}

/**
 * Modify the ray for depth of field.
 */
inline void ThinLensTiltShiftCamera::modifyRayForDof(
    Ray &ray,
    const CameraSample &sample
) const {
    // TODO: use an epsilon for float comparison
    if (lensRadius <= 0.) {
        return;
    }

    // Sample point on lens
    float lensU, lensV;
    ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
    lensU *= lensRadius;
    lensV *= lensRadius;

    // Compute point on plane of focus
    float ft; // ray parameter
    if (tiltEnabled) {
        ft = ((scheimpflugLineX - lensU) * focalLength) /
            ((ray.d.x * focalLength) - (hingeLineX - scheimpflugLineX) * ray.d.z);
    } else {
        ft = focalDistance / ray.d.z;
    }
    Point focusPoint = ray(ft);

    // Update the ray for effect of lens
    ray.o = Point(lensU, lensV, 0.f);
    Vector outputDirection = focusPoint - ray.o;
    if (focusPoint.z < 0.) {
        outputDirection = -outputDirection;
    }
    ray.d = Normalize(outputDirection);
}

ThinLensTiltShiftCamera *CreateThinLensTiltShiftCamera(
    const ParamSet &params,
    const AnimatedTransform &cam2world,
    Film *film
) {
    // Extract common camera parameters from _ParamSet_
    float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    if (shutterclose < shutteropen) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                shutterclose, shutteropen);
        swap(shutterclose, shutteropen);
    }
    float lensradius = params.FindOneFloat("lensradius", 0.f);
    float focallength = params.FindOneFloat("focallength", 1.);
    float imagedistance = params.FindOneFloat("imagedistance", 1.);
    float frame = params.FindOneFloat("frameaspectratio",
        float(film->xResolution)/float(film->yResolution));
    float screen[4];
    if (frame > 1.f) {
        screen[0] = -frame;
        screen[1] =  frame;
        screen[2] = -1.f;
        screen[3] =  1.f;
    }
    else {
        screen[0] = -1.f;
        screen[1] =  1.f;
        screen[2] = -1.f / frame;
        screen[3] =  1.f / frame;
    }
    int swi;
    const float *sw = params.FindFloat("screenwindow", &swi);
    if (sw && swi == 4) {
        memcpy(screen, sw, 4*sizeof(float));
    }
    float fov = params.FindOneFloat("fov", 90.);
    float halffov = params.FindOneFloat("halffov", -1.f);
    if (halffov > 0.f) {
        // hack for structure synth, which exports half of the full fov
        fov = 2.f * halffov;
    }
    // TODO: in future use a 2D vector for horizonal/vertical tilt and shift
    // - [*] horizontal tilt = tilt around the X direction of the film center
    // - vertical tilt = swing = tilt around the Y direction of the film center
    // - horizontal shift = shift in the X direction
    // - [*] vertical shift = shift in the Y direction
    // Note:
    // - the film tilt is given and stored in degrees and
    // - it rotates the film around its center in the X direction
    float filmtiltx = params.FindOneFloat("filmtiltx", 0.);
    // Note:
    // - film is shifted in the Y direction
    float filmshifty = params.FindOneFloat("filmshifty", 0.);
    return new ThinLensTiltShiftCamera(cam2world, screen, shutteropen,
        shutterclose, lensradius, focallength, imagedistance,
        filmtiltx, filmshifty, fov, film);
}
