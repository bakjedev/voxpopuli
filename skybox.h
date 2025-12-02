#pragma once
namespace Tmpl8 {
class Skybox
{
public:
	Skybox( const char* filepath );
	~Skybox();

	float3 render( const Ray& ray ) const;
private:
	SurfaceFloat* skydome;
};
}