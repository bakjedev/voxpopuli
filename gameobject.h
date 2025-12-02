#pragma once
namespace Tmpl8
{
	class GameObject
	{
	protected:
		float3 position;
	public:
		GameObject( float3 newPosition );
		virtual void Init() = 0;
		virtual void Update() = 0;
		virtual void Shutdown() = 0;

		void SetPosition( const float3 position );
		float3 GetPosition() const;
	};

}