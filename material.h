#pragma once

namespace Tmpl8
{
	class Material
	{
	public:
		Material( const char* texturePath = "" )
		{
			if ( texturePath && *texturePath != '\0' ) {
				texture = new Surface( texturePath );
				if ( !texture->pixels ) {
					printf( "Failed to load image\n" );
				}
			}
		}
		~Material()
		{
			delete texture;
			texture = nullptr;
		}

		Surface* texture{ nullptr };
		float reflectivity{ 0.0f };
		float refractivity{ 0.0f };
		float roughness{ 0.0f };
		float emission{ 0.0f };
		float ior{ 1.0f };
	};


}