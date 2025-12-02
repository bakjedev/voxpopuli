#pragma once

namespace Tmpl8
{
	class Light
	{
	public:
		enum Type { Directional, Point, Spot };

		void SetType( const Type t ) { type = t; }
		Type GetType() const { return type; }
		void SetPosition( const float3 p ) { position = p; }
		float3 GetPosition() const { return position; }
		void SetDirection( const float3 d ) { direction = d; }
		float3 GetDirection() const { return direction; }
		void SetAngle( const float a ) { angle = a; }
		float GetAngle() const { return angle; }

		void SetColor( const float3 rgb )
		{
			color = ( (uint)( rgb.x * 255 ) << 16 ) | ( (uint)( rgb.y * 255 ) << 8 ) | (uint)( rgb.z * 255 );
		}

		float3 GetColor() const
		{
			float r = ( ( color >> 16 ) & 0xFF ) / 255.0f;
			float g = ( ( color >> 8 ) & 0xFF ) / 255.0f;
			float b = ( color & 0xFF ) / 255.0f;
			return float3( r, g, b );
		}

	private:
		float3 position; // all
		float3 direction; // directional and spot
		float angle; // spot
		uint color; // all
		Type type; // all
	};
}