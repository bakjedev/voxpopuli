#pragma once
// default screen resolution
#define SCRWIDTH	512
#define SCRHEIGHT	320
// #define FULLSCREEN
#define DOUBLESIZE

namespace Tmpl8
{
	struct Plane
	{
		float3 normal = { 0.0f, 1.0f, 0.0f };
		float distance = 0.0f;
		float4 plane;

		Plane() = default;
		Plane( const float3& p, const float3& n ) : normal( normalize( n ) ), distance( -dot( normal, p ) )
		{
			plane = float4( normal, distance );
		}
	};

	struct Frustum
	{
		Plane right, left, top, bottom;
	};
	class Camera
	{
	public:
		Camera()
		{
			Initialize( float3( 0.375f, 0.85f, 0.3751f ), float3( 0.3751, -1.f, 0.3751f) );
		}

		void Initialize(float3 newCamPos, float3 newCamTarget)
		{
			aspect = (float)SCRWIDTH / (float)SCRHEIGHT;
			// setup a basic view frustum
			camPos = newCamPos;
			camTarget = newCamTarget;
			topLeft = float3( -aspect, 1, 0 );
			topRight = float3( aspect, 1, 0 );
			bottomLeft = float3( -aspect, -1, 0 );
		}

		Ray GetPrimaryRay( const float x, const float y )
		{
			const float u = x * ( 1.0f / SCRWIDTH );
			const float v = y * ( 1.0f / SCRHEIGHT );
			const float3 P = topLeft + u * ( topRight - topLeft ) + v * ( bottomLeft - topLeft );
			return Ray( camPos, normalize( P - camPos ) );
		}

		// credit: robin
		bool WorldSpaceToScreenSpace( const float3 worldPos, float2& screenPos )
		{
			float3 toPoint = worldPos - camPos;
			const float3 toPointNrm = normalize( toPoint );

			const float3 cameraDir = normalize( camTarget - camPos );
			const float3 viewRight = normalize( cross( float3( 0.0, 1.0, 0.0 ), cameraDir ) );
			const float3 viewUp = cross( cameraDir, viewRight );
			const float3 forward = cameraDir * 2.0f;

			float d = max( dot( cameraDir, toPointNrm ), 0.0f );
			if ( d < 0.01 )
				return false;

			d = 2.0f / d;

			toPoint = toPointNrm * d - forward;

			screenPos.x = dot( toPoint, viewRight );
			screenPos.y = dot( toPoint, viewUp );

			screenPos.x /= aspect;
			screenPos = screenPos * 0.5f + 0.5f;
			screenPos.y = 1 - screenPos.y;

			screenPos.x *= SCRWIDTH;
			screenPos.y *= SCRHEIGHT;

			if ( screenPos.x < 0 || screenPos.x > SCRWIDTH || screenPos.y < 0 || screenPos.y > SCRHEIGHT )
				return false;

			return true;
		}

		float3 RandomInUnitDisk()
		{
			float3 p;
			do
			{
				p = 2.0f * float3( RandomFloat(), RandomFloat(), 0 ) - float3( 1, 1, 0 );
			} while ( dot( p, p ) >= 1.0f );
			return p;
		}

		// depth of field - credit sven
		Ray GetSpecialPrimaryRay( const float x, const float y )
		{
			const float u = x * ( 1.0f / SCRWIDTH );
			const float v = y * ( 1.0f / SCRHEIGHT );

			const float3 pixelPosition = topLeft + u * ( topRight - topLeft ) + v * ( bottomLeft - topLeft );

			float3 initialDir = normalize( pixelPosition - camPos );

			float3 focalPlaneNormal = normalize( camTarget - camPos );
			float intersectionDistance = focalDistance / dot( focalPlaneNormal, initialDir );
			float3 focalPoint = camPos + initialDir * intersectionDistance;

			float3 lensPoint = camPos + RandomInUnitDisk() * lensRadius;

			float3 newDir = normalize( focalPoint - lensPoint );

			return Ray( lensPoint, newDir );
		}

		bool HandleInput( const float t )
		{
			if ( !WindowHasFocus() ) return false;
			float moveSpeed = 0.00025f * t;
			float lookSpeed = 0.0025f * t;
			float3 ahead = normalize( camTarget - camPos );
			float3 tmpUp( 0, 1, 0 );
			right = normalize( cross( tmpUp, ahead ) );
			up = normalize( cross( ahead, right ) );
			bool changed = false;
			if ( controls )
			{
				if ( IsKeyDown( GLFW_KEY_UP ) ) camTarget += lookSpeed * up, changed = true;
				if ( IsKeyDown( GLFW_KEY_DOWN ) ) camTarget -= lookSpeed * up, changed = true;
				if ( IsKeyDown( GLFW_KEY_LEFT ) ) camTarget -= lookSpeed * right, changed = true;
				if ( IsKeyDown( GLFW_KEY_RIGHT ) ) camTarget += lookSpeed * right, changed = true;

				ahead = normalize( camTarget - camPos );
				right = normalize( cross( tmpUp, ahead ) );
				up = normalize( cross( ahead, right ) );

				if ( IsKeyDown( GLFW_KEY_A ) ) camPos -= moveSpeed * right, changed = true;
				if ( IsKeyDown( GLFW_KEY_D ) ) camPos += moveSpeed * right, changed = true;
				if ( IsKeyDown( GLFW_KEY_W ) ) camPos += moveSpeed * ahead, changed = true;
				if ( IsKeyDown( GLFW_KEY_S ) ) camPos -= moveSpeed * ahead, changed = true;
				if ( IsKeyDown( GLFW_KEY_SPACE ) ) camPos += moveSpeed * up, changed = true;
				if ( IsKeyDown( GLFW_KEY_LEFT_CONTROL ) ) camPos -= moveSpeed * up, changed = true;
			}
			camTarget = camPos + ahead;
			ahead = normalize( camTarget - camPos );
			up = normalize( cross( ahead, right ) );
			right = normalize( cross( up, ahead ) );
			topLeft = camPos + 2 * ahead - aspect * right + up;
			topRight = camPos + 2 * ahead + aspect * right + up;
			bottomLeft = camPos + 2 * ahead - aspect * right - up;
			if ( !changed ) return false;
			return true;
		}

		// credit robin
		void CalculateFrustum()
		{
			float3 camDir = normalize( camTarget - camPos );
			// Calculate half vertical side and half horizontal side of the frustum
			const float halfVSide = 0.5f; // should be fov * 0.5f but fov is 1.0f
			const float halfHSide = halfVSide * aspect;

			// Initialize the frustum planes
			oldFrustum.right = { camPos, cross( camDir - right * halfHSide, up ) };
			oldFrustum.left = { camPos, cross( up, camDir + right * halfHSide ) };
			oldFrustum.top = { camPos, cross( right, camDir - up * halfVSide ) };
			oldFrustum.bottom = { camPos, cross( camDir + up * halfVSide, right ) };
		}

		// credit robin
		int2 PyramidSpace( float4 p )
		{
			// Calculate the distances d1, d2, d3, and d4 to the planes of the field of view pyramid
			float d1 = dot( p, oldFrustum.left.plane );
			float d2 = dot( p, oldFrustum.right.plane );
			float d3 = dot( p, oldFrustum.top.plane );
			float d4 = dot( p, oldFrustum.bottom.plane );

			// Calculate the screen space coordinates x and y using the formulas
			float x = d2 / ( d1 + d2 );
			float y = d4 / ( d3 + d4 );

			// Convert x and y to the range [0, 1]
			x = clamp( x, 0.0f, 1.0f );
			y = clamp( y, 0.0f, 1.0f );

			// Convert x and y to pixel positions
			float pixelX = x * SCRWIDTH;
			float pixelY = y * SCRHEIGHT;

			// Round pixelX and pixelY to the nearest integers

			int pixelPosX = min( static_cast<int>( std::round( pixelX ) ), SCRWIDTH - 1 );
			int pixelPosY = min( static_cast<int>( std::round( pixelY ) ), SCRHEIGHT - 1 );

			// Return the screen space coordinates as a int2
			return int2( pixelPosX, pixelPosY );
		}


		float aspect;
		float3 camPos, camTarget;
		float3 topLeft, topRight, bottomLeft;
		bool controls = false;

		float3 right, up;

		Frustum oldFrustum;

		float focalLength = 0.1f;
		float focalDistance = 0.35f;
		float lensRadius = 0.0f;
	};
}
