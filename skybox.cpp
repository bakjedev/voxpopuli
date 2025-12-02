#include "precomp.h"
#include "skybox.h"
// jacco's skydome implementation but with surfacefloat

Skybox::Skybox( const char* filepath )
{
	skydome = new SurfaceFloat( filepath );

	/*for ( int i = 0; i < skydome->width * skydome->height * skydome->channels; i++ )
		skydome->pixels[i] = sqrtf( skydome->pixels[i] ); */// Gamma Adjustment for Reduced HDR Range
}

Skybox::~Skybox()
{
	delete skydome;
	skydome = nullptr;
}

float3 Skybox::render( const Ray& ray ) const
{
	const float3 dir = normalize( ray.direction );

	// Sample Sky
	const uint u = static_cast<uint>( skydome->width * atan2f( dir.z, dir.x ) * INV2PI - 0.5f );
	const uint v = static_cast<uint>( skydome->height * acosf( dir.y ) * INVPI - 0.5f );

	const uint sky_idx = ( u + v * skydome->width ) % ( skydome->width * skydome->height );
	return float3( skydome->pixels[sky_idx * 3], skydome->pixels[sky_idx * 3 + 1], skydome->pixels[sky_idx * 3 + 2] );
}