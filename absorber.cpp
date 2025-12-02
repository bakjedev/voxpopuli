#include "precomp.h"
#include "absorber.h"

Absorber::Absorber( float3 position, uint primitiveIdx, float rate, uint color ) : GameObject( position ), primitiveIdx( primitiveIdx ), rate( rate ), color( color )
{
}

void Absorber::Init()
{
}

void Absorber::Update()
{
}

void Absorber::Shutdown()
{
}

void Absorber::SetPrimitiveIdx( uint idx )
{
	primitiveIdx = idx;
}

uint Absorber::GetPrimitiveIdx() const
{
	return primitiveIdx;
}

void Absorber::SetRate( const float newRate )
{
	rate = newRate;
}

float Absorber::GetRate() const
{
	return rate;
}

void Absorber::SetColor( const uint newColor )
{
	color = newColor;
}

uint Absorber::GetColor() const
{
	return color;
}
