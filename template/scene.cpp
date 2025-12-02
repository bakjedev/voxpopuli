#include "precomp.h"

float3 Ray::GetNormal() const
{
	// return the voxel normal at the nearest intersection
	const float3 I1 = ( origin + t * direction ) * WORLDSIZE; // our scene size is (1,1,1), so this scales each voxel to (1,1,1)
	const float3 fG = fracf( I1 );
	const float3 d = min3( fG, 1.0f - fG );
	const float mind = min( min( d.x, d.y ), d.z );
	const float3 sign = Dsign * 2 - 1;
	return float3( mind == d.x ? sign.x : 0, mind == d.y ? sign.y : 0, mind == d.z ? sign.z : 0 );
	// TODO:
	// - *only* in case the profiler flags this as a bottleneck:
	// - This function might benefit from SIMD.
}
//sven
float2 Ray::GetUV() const
{
	// Calculate the normal at the intersection point
	float3 normal = GetNormal();
	float3 hitPos = origin + t * direction;
	hitPos = hitPos * WORLDSIZE;

	// Calculate UV coordinates based on the face normal
	float u, v;
	if ( fabs( normal.x ) > fabs( normal.y ) && fabs( normal.x ) > fabs( normal.z ) )
	{
		// The hit is on a face perpendicular to the X-axis
		u = fmod( hitPos.z, 1.0f );
		v = fmod( hitPos.y, 1.0f );
	}
	else if ( fabs( normal.y ) > fabs( normal.x ) && fabs( normal.y ) > fabs( normal.z ) )
	{
		// The hit is on a face perpendicular to the Y-axis
		u = fmod( hitPos.x, 1.0f );
		v = fmod( hitPos.z, 1.0f );
	}
	else
	{
		// The hit is on a face perpendicular to the Z-axis
		u = fmod( hitPos.x, 1.0f );
		v = fmod( hitPos.y, 1.0f );
	}

	// Adjust UV coordinates based on the face normal direction
	if ( normal.x < 0 || normal.y < 0 || normal.z < 0 )
	{
		u = 1.0f - u;
		v = 1.0f - v;
	}

	// flip UV coordinates to match the texture (idk why this is necessary)
	v = 1.0f - v;

	// Clamp UV coordinates to avoid artifacts
	u = clamp( u, 0.0f, 1.0f );
	v = clamp( v, 0.0f, 1.0f );

	return float2( u, v );
}

Cube::Cube( const float3 pos, const float3 size, float3 color, uint material )
{
	// set cube bounds
	b[0] = pos;
	b[1] = pos + size;
	this->color = color;
	this->material = material;
}

float Cube::Intersect( const Ray& ray ) const
{
	// test if the ray intersects the cube
	const int signx = ray.direction.x < 0, signy = ray.direction.y < 0, signz = ray.direction.z < 0;
	float tmin = ( b[signx].x - ray.origin.x ) * ray.rD.x;
	float tmax = ( b[1 - signx].x - ray.origin.x ) * ray.rD.x;
	const float tymin = ( b[signy].y - ray.origin.y ) * ray.rD.y;
	const float tymax = ( b[1 - signy].y - ray.origin.y ) * ray.rD.y;
	if ( tmin > tymax || tymin > tmax ) goto miss;
	tmin = max( tmin, tymin ), tmax = min( tmax, tymax );
	const float tzmin = ( b[signz].z - ray.origin.z ) * ray.rD.z;
	const float tzmax = ( b[1 - signz].z - ray.origin.z ) * ray.rD.z;
	if ( tmin > tzmax || tzmin > tmax ) goto miss; // yeah c has 'goto' ;)
	if ( ( tmin = max( tmin, tzmin ) ) > 0 ) return tmin;
miss:
	return 1e34f;
}

bool Cube::Contains( const float3& pos ) const
{
	// test if pos is inside the cube
	return pos.x >= b[0].x && pos.y >= b[0].y && pos.z >= b[0].z &&
		pos.x <= b[1].x && pos.y <= b[1].y && pos.z <= b[1].z;
}

Sphere::Sphere( const float3 pos, const float radius, float3 color, uint material )
	: pos( pos ), radius( radius ), color( color ), material( material )
{
}


bool Sphere::Intersect( const Ray& ray, float& t, float3& n ) const
{
	float3 AC = ray.origin - pos;
	float a = dot( ray.direction, ray.direction );
	float b = 2.0f * dot( ray.direction, AC );
	float c = dot( AC, AC ) - radius * radius;

	float discriminant = b * b - 4 * a * c;

	if ( discriminant > 0 )
	{
		float nrm_sign = c < 0.0f ? -1.f : 1.f;

		float tempSolve = ( -b - nrm_sign * sqrtf( discriminant ) ) / ( 2.f * a );

		if ( tempSolve > 0.f && tempSolve < 100000.f )
		{
			t = tempSolve;
			float3 intersection = ray.origin + ray.direction * t;
			n = ( ( intersection - pos ) / radius );
			return true;
		}
	}

	return false;
}


BVH::BVH( Primitive* newPrimitives, PrimitiveType* newTypes, int N ) : N( N )
{
	primitiveData = newPrimitives;
	primitiveTypes = newTypes;
	bvhNode = new BVHNode[N * 2];
	primIdx = new uint[N];
	BuildBVH();
}

BVH::~BVH()
{
	delete[] primitiveData;
	primitiveData = nullptr;

	delete[] primitiveTypes;
	primitiveTypes = nullptr;

	delete[] primIdx;
	primIdx = nullptr;

	delete[] bvhNode;
	bvhNode = nullptr;
}


void BVH::BuildBVH()
{
	for ( int i = 0; i < N; i++ ) primIdx[i] = i;
	BVHNode& root = bvhNode[rootNodeIdx];
	root.firstPrimIdx = 0, root.primCount = N;
	UpdateNodeBounds( rootNodeIdx );
	Subdivide( rootNodeIdx );
}

void BVH::UpdateNodeBounds( uint nodeIdx )
{
	BVHNode& node = bvhNode[nodeIdx];
	node.aabbMin = float3( 1e30f );
	node.aabbMax = float3( -1e30f );
	for ( uint first = node.firstPrimIdx, i = 0; i < node.primCount; i++ )
	{
		uint leafPrimIdx = primIdx[first + i];

		if ( primitiveTypes[leafPrimIdx] == SPHERE )
		{
			Sphere& leafSphere = primitiveData[leafPrimIdx].sphere;
			node.aabbMin = fminf( node.aabbMin, leafSphere.pos - leafSphere.radius );
			node.aabbMax = fmaxf( node.aabbMax, leafSphere.pos + leafSphere.radius );
		}
		else if ( primitiveTypes[leafPrimIdx] == CUBE )
		{
			Cube& cube = primitiveData[leafPrimIdx].cube;
			node.aabbMin = fminf( node.aabbMin, cube.b[0] );
			node.aabbMax = fmaxf( node.aabbMax, cube.b[1] );
		}
	}
}

void BVH::Subdivide( uint nodeIdx )
{
	// terminate recursion
	BVHNode& node = bvhNode[nodeIdx];
	if ( node.primCount <= 2 ) return;
	// determine split axis and position
	float3 extent = node.aabbMax - node.aabbMin;
	int axis = 0;
	if ( extent.y > extent.x ) axis = 1;
	if ( extent.z > extent[axis] ) axis = 2;
	float splitPos = node.aabbMin[axis] + extent[axis] * 0.5f;
	// in-place partition
	int i = node.firstPrimIdx;
	int j = i + node.primCount - 1;
	while ( i <= j )
	{
		if ( primitiveTypes[primIdx[i]] == SPHERE )
		{
			if ( primitiveData[primIdx[i]].sphere.pos[axis] < splitPos )
				i++;
			else
				swap( primIdx[i], primIdx[j--] );
		}
		else if ( primitiveTypes[primIdx[i]] == CUBE )
		{
			if ( primitiveData[primIdx[i]].cube.b[0][axis] < splitPos )
				i++;
			else
				swap( primIdx[i], primIdx[j--] );
		}
	}
	// abort split if one of the sides is empty
	int leftCount = i - node.firstPrimIdx;
	if ( leftCount == 0 || static_cast<uint>(leftCount) == node.primCount ) return;
	// create child nodes
	int leftChildIdx = nodesUsed++;
	int rightChildIdx = nodesUsed++;
	bvhNode[leftChildIdx].firstPrimIdx = node.firstPrimIdx;
	bvhNode[leftChildIdx].primCount = leftCount;
	bvhNode[rightChildIdx].firstPrimIdx = i;
	bvhNode[rightChildIdx].primCount = node.primCount - leftCount;
	node.leftNode = leftChildIdx;
	node.primCount = 0;
	UpdateNodeBounds( leftChildIdx );
	UpdateNodeBounds( rightChildIdx );
	// recurse
	Subdivide( leftChildIdx );
	Subdivide( rightChildIdx );
}

void BVH::Add( PrimitiveType type, Primitive primitive )
{
	int newSize = N + 1;

	Primitive* newPrimitives = new Primitive[newSize];
	PrimitiveType* newPrimitiveTypes = new PrimitiveType[newSize];

	for ( int i = 0; i < N; i++ )
	{
		newPrimitives[i] = primitiveData[i];
		newPrimitiveTypes[i] = primitiveTypes[i];
	}

	BVHNode* newBVHNode = new BVHNode[newSize * 2];
	uint* newPrimIdx = new uint[newSize];

	for ( int i = 0; i < N; i++ )
	{
		newPrimIdx[i] = primIdx[i];
	}

	newPrimitives[N] = primitive;
	newPrimitiveTypes[N] = type;

	delete[] primitiveData;
	delete[] primitiveTypes;
	delete[] bvhNode;
	delete[] primIdx;

	primitiveData = newPrimitives;
	primitiveTypes = newPrimitiveTypes;
	bvhNode = newBVHNode;
	primIdx = newPrimIdx;

	N = newSize;

	BuildBVH();
}

void BVH::SetPosition( uint idx, float3 pos )
{
	if ( primitiveTypes[idx] == SPHERE )
	{
		primitiveData[idx].sphere.pos = pos;
	}
	else if ( primitiveTypes[idx] == CUBE )
	{
		float3 size = primitiveData[idx].cube.b[1] - primitiveData[idx].cube.b[0];
		primitiveData[idx].cube.b[0] = pos;
		primitiveData[idx].cube.b[1] = pos + size;
	}
	for ( int i = 0; i < N; i++ )
	{
		UpdateNodeBounds( i );
	}
}

Scene::Scene()
{
#ifdef USEBVH
	int size = 3;
	Primitive* newPrimitives = new Primitive[size];
	PrimitiveType* newPrimitiveTypes = new PrimitiveType[size];

	newPrimitives[0].sphere = Sphere( float3( 0.45f, 0.055f, 0.4f ),  0.0225f, float3( 1.0f, 0.3f, 0.3f ), 1 );
	newPrimitiveTypes[0] = SPHERE;
	newPrimitives[1].sphere = Sphere( float3( 0.33f, 0.055f, 0.1f ), 0.0225f, float3( 0.3f, 1.0f, 0.3f ), 1 );
	newPrimitiveTypes[1] = SPHERE;
	newPrimitives[2].sphere = Sphere( float3( 0.45f, 0.055f, 0.6f ), 0.0225f, float3( 0.3f, 0.3f, 1.0f ), 1 );
	newPrimitiveTypes[2] = SPHERE;


	bvh = new BVH( newPrimitives, newPrimitiveTypes, size );
	// the voxel world sits in a 1x1x1 cube
	cube = Cube( float3( 0, 0, 0 ), float3( 1, 1, 1 ) );

	// initialize the scene using Perlin noise, parallel over z
	grid = (uint*)MALLOC64( GRIDSIZE3 * sizeof( uint ) );
	memset( grid, 0, GRIDSIZE3 * sizeof( uint ) );


	
	
#endif
	
#ifndef SPHERES
	// the voxel world sits in a 1x1x1 cube
	cube = Cube( float3( 0, 0, 0 ), float3( 1, 1, 1 ) );
	// initialize the scene using Perlin noise, parallel over z
	grid = (uint*)MALLOC64( GRIDSIZE3 * sizeof( uint ) );
	memset( grid, 0, GRIDSIZE3 * sizeof( uint ) );

	for ( int z = 0; z < WORLDSIZE; z++ )
	{
		const float fz = (float)z / WORLDSIZE;
		for ( int y = 0; y < WORLDSIZE; y++ )
		{
			const float fy = (float)y / WORLDSIZE;
			float fx = 0;
			for ( int x = 0; x < WORLDSIZE; x++, fx += 1.0f / WORLDSIZE )
			{
				const float n = noise3D( fx, fy, fz );
				Set( x, y, z, n > 0.09f ? float3( 1.0f, 1.0f, 1.0f ) : float3( 0.0f ), 0 );
			}
		}
	}

	int size = 8;
	for ( int i = 0; i < size; i++ )
	{
		for ( int j = 0; j < size; j++ )
		{
			Set( 64 + i, 90, 64 + j, float3( 1.0f, 0.0f, 1.0f ), 0 ); // floor
			Set( 64 + i, 90 + size, 64 + j, float3( 1.0f, 1.0f, 1.0f ), 1 ); // roof

			Set( 64 + i, 91 + j, 64, float3( 1.0f, 1.0f, 1.0f ), 1 ); // wall
			Set( 64, 91 + j, 64 + i, float3( 1.0f, 1.0f, 1.0f ), 1 ); // wall
			Set( 64 + i, 91 + j, 64 + size, float3( 1.0f, 1.0f, 1.0f ), 1 ); // wall
			Set( 64 + size, 91 + j, 64 + i, float3( 1.0f, 1.0f, 1.0f ), 1 ); // wall
		}
	}
#endif // !SPHERES
}

Scene::~Scene()
{
	// Delete dynamically allocated grid
	delete[] grid;
	grid = nullptr;

	// No explicit deletion for cube, assuming it doesn't need dynamic allocation

	// Delete dynamically allocated BVH
	delete bvh;
	bvh = nullptr;
}



void Scene::Set( const uint x, const uint y, const uint z, const float3& c, const uint m )
{
	int index = x + y * GRIDSIZE + z * GRIDSIZE2;
	uint value;
	if ( c != float3( 0.0f ) )
	{
		uint r = static_cast<uint>( c.x * 255.0f ) & 0xFF;
		uint g = static_cast<uint>( c.y * 255.0f ) & 0xFF;
		uint b = static_cast<uint>( c.z * 255.0f ) & 0xFF;

		value = ( r << 24 ) | ( g << 16 ) | ( b << 8 ) | m;
	}
	else
	{
		value = 0;
	}

	grid[index] = value;
}


bool Scene::Setup3DDDA( const Ray& ray, DDAState& state ) const
{
	// if ray is not inside the world: advance until it is
	state.t = 0;
	if ( !cube.Contains( ray.origin ) )
	{
		state.t = cube.Intersect( ray );
		if ( state.t > 1e33f ) return false; // ray misses voxel data entirely
	}
	// setup amanatides & woo - assume world is 1x1x1, from (0,0,0) to (1,1,1)
	static const float cellSize = 1.0f / GRIDSIZE;
	state.step = make_int3( 1 - ray.Dsign * 2 );
	const float3 posInGrid = GRIDSIZE * ( ray.origin + ( state.t + 0.00005f ) * ray.direction );
	const float3 gridPlanes = ( ceilf( posInGrid ) - ray.Dsign ) * cellSize;
	const int3 P = clamp( make_int3( posInGrid ), 0, GRIDSIZE - 1 );
	state.X = P.x, state.Y = P.y, state.Z = P.z;
	state.tdelta = cellSize * float3( state.step ) * ray.rD;
	state.tmax = ( gridPlanes - ray.origin ) * ray.rD;
	// proceed with traversal
	return true;
}

void Scene::IntersectBVH( Ray& ray, const uint nodeIdx )
{
	BVHNode& node = bvh->bvhNode[nodeIdx];
	if ( !IntersectAABB( ray, node.aabbMin, node.aabbMax ) ) return;
	if ( node.isLeaf() )
	{
		for ( uint i = 0; i < node.primCount; i++ )
		{
			uint primIdx = bvh->primIdx[node.firstPrimIdx + i];
			PrimitiveType type = bvh->primitiveTypes[primIdx];
			if ( type == SPHERE )
			{
				Sphere& hitSphere = bvh->primitiveData[primIdx].sphere;
				float t;
				float3 n;
				if ( hitSphere.Intersect( ray, t, n ) )
				{
					if ( t < ray.t )
					{
						ray.t = t;
						ray.sphereIdx = primIdx;
						float3 c = hitSphere.color;
						uint r = static_cast<uint>( c.x * 255.0f ) & 0xFF;
						uint g = static_cast<uint>( c.y * 255.0f ) & 0xFF;
						uint b = static_cast<uint>( c.z * 255.0f ) & 0xFF;

						ray.hit = ( r << 24 ) | ( g << 16 ) | ( b << 8 ) | hitSphere.material;

						ray.N = n;
						ray.type = 0;
					}
				}
			}
			else if ( type == CUBE )
			{
				Cube& hitCube = bvh->primitiveData[primIdx].cube;
				float t = hitCube.Intersect( ray );
				if ( t  < ray.t )
				{
					ray.t = t;
					float3 c = hitCube.color;
					uint r = static_cast<uint>( c.x * 255.0f ) & 0xFF;
					uint g = static_cast<uint>( c.y * 255.0f ) & 0xFF;
					uint b = static_cast<uint>( c.z * 255.0f ) & 0xFF;

					ray.hit = ( r << 24 ) | ( g << 16 ) | ( b << 8 ) | hitCube.material;

					ray.N = hitCube.CalculateNormal( ray.origin + ray.t * ray.direction );
					ray.type = 1;
				}
			}
		}
	}
	else
	{
		IntersectBVH( ray, node.leftNode );
		IntersectBVH( ray, node.leftNode + 1 );
	}
}

void Scene::FindNearest( Ray& ray )
{
	Ray r = ray;
	IntersectBVH( r, bvh->rootNodeIdx );


	// setup Amanatides & Woo grid traversal
	DDAState s, bs;
	if ( !Setup3DDDA( ray, s ) ) {
		ray.t = r.t;
		ray.hit = r.hit;
		ray.N = r.N;
		ray.sphereIdx = r.sphereIdx;
		ray.type = r.type;
		return;
	}
	// start stepping
	while ( 1 ) {
		const uint cell = grid[s.X + s.Y * GRIDSIZE + s.Z * GRIDSIZE2];
		if ( cell && s.t < r.t ) {
			ray.t = s.t;
			ray.hit = cell;
			ray.N = ray.GetNormal();
			ray.type = 1;
			return;
		}
		if ( s.tmax.x < s.tmax.y ) {
			if ( s.tmax.x < s.tmax.z ) {
				s.t = s.tmax.x;
				s.X += s.step.x;
				if ( s.X >= GRIDSIZE )
					break;
				s.tmax.x += s.tdelta.x;
			} else {
				s.t = s.tmax.z;
				s.Z += s.step.z;
				if ( s.Z >= GRIDSIZE )
					break;
				s.tmax.z += s.tdelta.z;
			}
		} else {
			if ( s.tmax.y < s.tmax.z ) {
				s.t = s.tmax.y;
				s.Y += s.step.y;
				if ( s.Y >= GRIDSIZE )
					break;
				s.tmax.y += s.tdelta.y;
			} else {
				s.t = s.tmax.z;
				s.Z += s.step.z;
				if ( s.Z >= GRIDSIZE )
					break;
				s.tmax.z += s.tdelta.z;
			}
		}
	}

	ray.t = r.t;
	ray.hit = r.hit;
	ray.N = r.N;
	ray.sphereIdx = r.sphereIdx;
	ray.type = r.type;
}

bool Scene::IsOccluded( const Ray& ray )
{
	Ray r = ray;
	IntersectBVH( r, bvh->rootNodeIdx );

	// setup Amanatides & Woo grid traversal
	DDAState s, bs;
	if ( !Setup3DDDA( ray, s ) ) return false;
	// start stepping
	while ( s.t < ray.t ) {
		const uint cell = grid[s.X + s.Y * GRIDSIZE + s.Z * GRIDSIZE2];
		if ( cell && s.t < r.t ) {
			return true;
		}
		if ( s.tmax.x < s.tmax.y ) {
			if ( s.tmax.x < s.tmax.z ) { if ( ( s.X += s.step.x ) >= GRIDSIZE ) return false; s.t = s.tmax.x, s.tmax.x += s.tdelta.x; } else { if ( ( s.Z += s.step.z ) >= GRIDSIZE ) return false; s.t = s.tmax.z, s.tmax.z += s.tdelta.z; }
		} else {
			if ( s.tmax.y < s.tmax.z ) { if ( ( s.Y += s.step.y ) >= GRIDSIZE ) return false; s.t = s.tmax.y, s.tmax.y += s.tdelta.y; } else { if ( ( s.Z += s.step.z ) >= GRIDSIZE ) return false; s.t = s.tmax.z, s.tmax.z += s.tdelta.z; }
		}
	}

	if ( r.t < ray.t )
		return true;

	return false;
}

bool Scene::IntersectAABB( const Ray& ray, const float3 bmin, const float3 bmax )
{
	float tx1 = ( bmin.x - ray.origin.x ) / ray.direction.x, tx2 = ( bmax.x - ray.origin.x ) / ray.direction.x;
	float tmin = min( tx1, tx2 ), tmax = max( tx1, tx2 );
	float ty1 = ( bmin.y - ray.origin.y ) / ray.direction.y, ty2 = ( bmax.y - ray.origin.y ) / ray.direction.y;
	tmin = max( tmin, min( ty1, ty2 ) ), tmax = min( tmax, max( ty1, ty2 ) );
	float tz1 = ( bmin.z - ray.origin.z ) / ray.direction.z, tz2 = ( bmax.z - ray.origin.z ) / ray.direction.z;
	tmin = max( tmin, min( tz1, tz2 ) ), tmax = min( tmax, max( tz1, tz2 ) );
	return tmax >= tmin && tmin < ray.t && tmax > 0;
}