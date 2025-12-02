#pragma once


#define EPSILON 0.0001f

#define MAX_LIGHTS 30

namespace Tmpl8
{
	class Skybox;
	class Game;
	class Material;
	class Light;

	struct DebugLine
	{
		DebugLine( float3 start, float3 end, int color ) : start( start ), end( end ), color( color ) {}
		float3 start, end;
		int color;
	};

	class Renderer : public TheApp
	{
	public:
		// game flow methods
		void Init();
		void Tick( const float deltaTime );
		void UI();
		void Shutdown();
		void MouseMove( int x, int y ) { mousePos.x = x; mousePos.y = y; }
		void ResetAccumulator() { reprojectAlpha = 0.9f; }
		Scene scene;
		Light* lights[MAX_LIGHTS] = { nullptr };
		Camera camera;
		int2 mousePos{ 0, 0 };
		int lightCount{ 0 };
		bool debug{ false };
	private:

		//float fps;
		float avg{ 10 }, alpha{ 1 }, fps;

		Skybox* skybox;
		float4* newAccumulator;
		float4* oldAccumulator;
		Material* materials[MATERIALS] = { nullptr };
		vector<DebugLine*> debugLines;

		Game* game;

		// settings
		bool toneMapping{ true };
		bool reprojection{ true };
		float reprojectAlpha{ 0.9f };
		uint activeMaterial{ 1 };

		// other methods
		float3 AcesApproximated( float3 C ) const;
		void AcesApproximatedSSE( const float3* C, const int numC, float3* result );


		// light transport
		int maxRayDepth{ 4 };

		// ui methods
		void PerformanceUI() const;
		void SettingsUI();

		//light transport methods
		float3 Trace( Ray& ray, int depth, const bool draw = false );
		float3 RandomDirectionOnHemisphere( const float3& normal ) const;
		float3 ComputeDirectLighting( const Ray& ray, const float3& normal );
		int PickLight( const float3& intersectionPoint, float* weights ) const;
	};

}