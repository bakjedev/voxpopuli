#pragma once

// high level settings
// #define TWOLEVEL
#define WORLDSIZE	32 // power of 2. Warning: max 512 for a 512x512x512x4 bytes = 512MB world!

#define SPHERES
#define USEBVH

#define MATERIALS 4

// low-level / derived
#define WORLDSIZE2	(WORLDSIZE*WORLDSIZE)
#ifdef TWOLEVEL
// feel free to replace with whatever suits your two-level implementation,
// should you chose this challenge.
#define BRICKSIZE	1
#define BRICKSIZE2	(BRICKSIZE*BRICKSIZE)
#define BRICKSIZE3	(BRICKSIZE*BRICKSIZE*BRICKSIZE)
#define GRIDSIZE	(WORLDSIZE/BRICKSIZE)
#define VOXELSIZE	(1.0f/WORLDSIZE)
#else
#define GRIDSIZE	WORLDSIZE
#endif
#define GRIDSIZE2	(GRIDSIZE*GRIDSIZE)
#define GRIDSIZE3	(GRIDSIZE*GRIDSIZE*GRIDSIZE)

namespace Tmpl8
{

	class Ray
	{
	public:
		Ray() = default;
		Ray( const float3 origin, const float3 direction, const float rayLength = 1e34f )
			: origin( origin ), direction( direction ), t( rayLength )
		{
			// calculate reciprocal ray direction for triangles and AABBs
			// TODO: prevent NaNs - or don't
			rD = float3( 1 / direction.x, 1 / direction.y, 1 / direction.z );


			uint xsign = *(uint*)&direction.x >> 31;
			uint ysign = *(uint*)&direction.y >> 31;
			uint zsign = *(uint*)&direction.z >> 31;
			Dsign = ( float3( (float)xsign * 2 - 1, (float)ysign * 2 - 1, (float)zsign * 2 - 1 ) + 1 ) * 0.5f;

		}
		float3 IntersectionPoint() const { return origin + t * direction; }
		float3 GetNormal() const;
		float2 GetUV() const;
		// ray data
		float3 origin;					// ray origin
		float3 rD;					// reciprocal ray direction
		float3 direction = float3( 0 );		// ray direction
		float t = 1e34f;			// ray length
		float3 Dsign = float3( 1 );	// inverted ray direction signs, -1 or 1
		uint hit = 0;				// 32-bit RGBM (M = material) of the hit (previously called voxel)
		// i would have liked for these not to be here, but i need this info in the trace function...
		float3 N;					// normal
		char type;					// 0 = sphere, 1 = cube/voxel
		uint sphereIdx;				// index of the sphere that was hit

		float3 GetColor()
		{
			return { ( hit >> 24u & 0xFF ) / 255.0f,
					 ( hit >> 16u & 0xFF ) / 255.0f,
					 ( hit >> 8u & 0xFF ) / 255.0f };
		}

		unsigned int GetMaterial()
		{
			if ( ( hit & 0xFF ) < MATERIALS )
				return hit & 0xFF;
			else
				return 0;
		}

	private:
		// min3 is used in normal reconstruction.
		__inline static float3 min3( const float3& a, const float3& b )
		{
			return float3( min( a.x, b.x ), min( a.y, b.y ), min( a.z, b.z ) );
		}
	};

	class Cube
	{
	public:
		Cube() = default;
		Cube( const float3 pos, const float3 size, float3 color = float3(0.0f), uint material = 0);
		float Intersect( const Ray& ray ) const;
		bool Contains( const float3& pos ) const;
		float3 b[2];
		float3 color;
		int material;

		// Inside Cube class
		float3 CalculateNormal( const float3& hitPoint ) const
		{
			float epsilon = 0.00001f; // Small epsilon value to handle floating-point inaccuracies

			// Check which face of the cube was hit
			if ( abs( hitPoint.x - b[0].x ) < epsilon )
				return float3( -1.0f, 0.0f, 0.0f ); // Left face
			else if ( abs( hitPoint.x - b[1].x ) < epsilon )
				return float3( 1.0f, 0.0f, 0.0f ); // Right face
			else if ( abs( hitPoint.y - b[0].y ) < epsilon )
				return float3( 0.0f, -1.0f, 0.0f ); // Bottom face
			else if ( abs( hitPoint.y - b[1].y ) < epsilon )
				return float3( 0.0f, 1.0f, 0.0f ); // Top face
			else if ( abs( hitPoint.z - b[0].z ) < epsilon )
				return float3( 0.0f, 0.0f, -1.0f ); // Front face
			else
				return float3( 0.0f, 0.0f, 1.0f ); // Back face
		}


	};

	class Sphere
	{
	public:
		Sphere() = default;
		Sphere( const float3 pos, const float radius, float3 color, uint material );
		bool Intersect( const Ray& ray, float& t, float3& n ) const;
		float3 pos;
		float radius;
		float3 color;
		int material;
	};

	enum PrimitiveType { SPHERE, CUBE };

	union Primitive
	{
		Sphere sphere;
		Cube cube;
	};

	struct BVHNode
	{
		float3 aabbMin, aabbMax;
		uint leftNode, firstPrimIdx, primCount;
		bool isLeaf() { return primCount > 0; }
	};

	class BVH
	{
	public:
		BVH() = default;
		BVH( Primitive* newPrimitives, PrimitiveType* newTypes, int N );
		~BVH();
		void BuildBVH();
		void UpdateNodeBounds( uint nodeIdx );
		void Subdivide( uint nodeIdx );

		void Add( PrimitiveType type, Primitive primitive );
		void SetPosition( uint idx, float3 pos );

		Primitive* primitiveData; // Array of primitive data
		PrimitiveType* primitiveTypes; // Array of primitive types
		uint* primIdx;

		uint rootNodeIdx = 0, nodesUsed = 1;

		BVHNode* bvhNode;

		int N;
	};

	class Scene
	{
	public:
		struct DDAState
		{
			int3 step;				// 16 bytes
			uint X, Y, Z;			// 12 bytes
			float t;				// 4 bytes
			float3 tdelta;
			float dummy1 = 0;		// 16 bytes
			float3 tmax;
			float dummy2 = 0;		// 16 bytes, 64 bytes in total
		};
		Scene();
		~Scene();
		void FindNearest( Ray& ray );
		bool IsOccluded( const Ray& ray );
		void Set( const uint x, const uint y, const uint z, const float3& c, const uint m );
		void IntersectBVH( Ray& ray, const uint nodeIdx );
		unsigned int* grid;
		Cube cube;
		BVH* bvh;
	private:
		bool Setup3DDDA( const Ray& ray, DDAState& state ) const;
		bool IntersectAABB( const Ray& ray, const float3 bmin, const float3 bmax );
	};

}