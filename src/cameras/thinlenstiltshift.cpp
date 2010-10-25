
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
    // Generate raster and camera samples
    // sample in Raster space
    Point pointRaster(sample.imageX, sample.imageY, 0);
    Point pointCamera;
    // convert it to a sample in Camera space
    // pointCamera->z should be mapped from 0 to the near plane
    RasterToCamera(pointRaster, &pointCamera);
    // project the point from the near plane to the untilted film
    float filmT = filmShift.z / pointCamera.z;
    Point pointFilm = Point(pointCamera.x * filmT, pointCamera.y *filmT, filmShift.z);
    // tilt the film
    Point pointFilmTilted;
    filmTiltRotate(pointFilm, &pointFilmTilted);
    *ray = Ray(Point(0,0,0), Vector(pointFilm), 0.f, INFINITY);
    // Modify ray for depth of field
    if (lensRadius > 0.) {
        // Sample point on lens
        float lensU, lensV;
        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
        lensU *= lensRadius;
        lensV *= lensRadius;

        // Compute point on plane of focus
        float ft;
        if (tiltEnabled) {
            ft = ((scheimpflugLineX - lensU) * focalLength) /
                ((ray->d.x * focalLength) - (hingeLineX - scheimpflugLineX) * ray->d.z);
        } else {
            ft = focalDistance / ray->d.z;
        }
        Point Pfocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point(lensU, lensV, 0.f);
        ray->d = Normalize(Pfocus - ray->o);
    }
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    CameraToWorld(*ray, ray);
    return 1.f;
}


float ThinLensTiltShiftCamera::GenerateRayDifferential(
    const CameraSample &sample,
    RayDifferential *ray
) const {
    // Generate raster and camera samples
    Point pointRaster(sample.imageX, sample.imageY, 0);
    Point pointCamera;
    RasterToCamera(pointRaster, &pointCamera);
    // project the point from the near plane to the untilted film
    float filmT = filmShift.z / pointCamera.z;
    Point pointFilm = Point(pointCamera.x * filmT, pointCamera.y *filmT, filmShift.z);
    // tilt the film
    Point pointFilmTilted;
    filmTiltRotate(pointFilm, &pointFilmTilted);
    Vector dir = Normalize(Vector(pointFilmTilted.x, pointFilmTilted.y, pointFilmTilted.z));
    *ray = RayDifferential(Point(0,0,0), dir, 0.f, INFINITY);
    // Modify ray for depth of field
    if (lensRadius > 0.) {
        // Sample point on lens
        float lensU, lensV;
        ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
        lensU *= lensRadius;
        lensV *= lensRadius;

        // Compute point on plane of focus
        float ft;
        if (tiltEnabled) {
            ft = ((scheimpflugLineX - lensU) * focalLength) /
                ((ray->d.x * focalLength) - (hingeLineX - scheimpflugLineX) * ray->d.z);
        } else {
            ft = focalDistance / ray->d.z;
        }
        Point Pfocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point(lensU, lensV, 0.f);
        Vector outputDirection = Pfocus - ray->o;
        if (Pfocus.z < 0.) {
            outputDirection = -outputDirection;
        }
        ray->d = Normalize(outputDirection);
    }

    // Compute offset rays for _ThinLensTiltShiftCamera_ ray differentials
    ray->rxOrigin = ray->ryOrigin = ray->o;
    ray->rxDirection = Normalize(Vector(pointCamera) + dxCamera);
    ray->ryDirection = Normalize(Vector(pointCamera) + dyCamera);
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    CameraToWorld(*ray, ray);
    ray->hasDifferentials = true;
    return 1.f;
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

float ThinLensTiltShiftCamera::computeFocalDistance(float focalLength, float imageDistance) {
    return 1 / ((1 / focalLength) - (1 / imageDistance));
}