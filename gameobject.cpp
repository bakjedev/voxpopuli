#include "precomp.h"
#include "gameobject.h"

GameObject::GameObject( float3 position ) : position( position )
{
}

void GameObject::SetPosition( const float3 newPosition )
{
	position = newPosition;
}

float3 GameObject::GetPosition() const
{
	return position;
}
