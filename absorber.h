#pragma once
#include "gameobject.h"
namespace Tmpl8
{
    class Absorber : public GameObject
    {
    private:
        float rate;
        uint primitiveIdx;
        uint color;
    public:
        Absorber( float3 position, uint primitiveIdx, float rate, uint color );
        void Init() override;
        void Update() override;
        void Shutdown() override;

        void SetPrimitiveIdx( const uint idx );
        uint GetPrimitiveIdx() const;

        void SetRate( const float newRate );
        float GetRate() const;

        void SetColor( const uint newColor );
        uint GetColor() const;

    };
}

