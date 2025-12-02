#include "precomp.h"
#include "game.h"
#include "light.h"
#include "gameobject.h"
#include "absorber.h"

void Game::Init()
{
	for ( int i = 0; i < 24; i++ )
	{
		for ( int j = 0; j < 24; j++ )
		{
			renderer->scene.Set( i, 64, j, float3( 1.0f, 0.8f, 1.0f ), 0 );
			if ( i == 0 || j == 0 || i == 23 || j == 23 )
			{
				renderer->scene.Set( i, 65, j, float3( 1.0f, 1.0f, 0.8f ), 0 );
				renderer->scene.Set( i, 66, j, float3( 1.0f, 1.0f, 0.8f ), 0 );
			}
		}
	}


	gameObjects[0] = new Absorber( float3( 0.3f, 0.05f, 0.3f ), 0, 0.2f, 0 );
	renderer->scene.bvh->SetPosition( 0, gameObjects[0]->GetPosition() );
	gameObjects[1] = new Absorber( float3( 0.3f, 0.05f, 0.6f ), 1, 0.2f, 1 );
	renderer->scene.bvh->SetPosition( 1, gameObjects[1]->GetPosition() );
	gameObjects[2] = new Absorber( float3( 0.1f, 0.05f, 0.2f ), 2, 0.2f, 2 );
	renderer->scene.bvh->SetPosition( 2, gameObjects[2]->GetPosition() );

	LoadLevel( 0 );
}

void Game::Update( const float deltaTime )
{
	// make light 3 direction go in circles around its position
	static float angle = 0.0f;
	float3 newDirection( 0.0f );
	float radius = 1.0f;
	float3 spotLightPos = renderer->lights[3]->GetPosition();
	newDirection.x = spotLightPos.x + radius * cos( angle );
	newDirection.z = spotLightPos.z + radius * sin( angle );
	renderer->lights[3]->SetDirection( newDirection );

	angle += 0.001f * deltaTime;
	if ( angle > 6.28319f ) angle = 0.0f;
	if ( testScene ) renderer->ResetAccumulator();


	if ( !playing ) return;
	float distances[MAX_LIGHTS];
	for ( int i = 0; i < MAX_LIGHTS; i++ )
	{
		if ( i >= renderer->lightCount ) break;

		float2 mousePos = renderer->mousePos;
		Ray primaryRay = renderer->camera.GetPrimaryRay( mousePos.x, mousePos.y );
		renderer->scene.FindNearest( primaryRay );
		float3 intersectionPoint = primaryRay.IntersectionPoint();

		float3 toNewPoint = intersectionPoint - renderer->lights[i]->GetPosition();

		distances[i] = length( toNewPoint );
	}

	// lowest distance is selected
	selectedLight = 0;
	for ( int i = 1; i < MAX_LIGHTS; i++ )
	{
		if ( i >= renderer->lightCount ) break;

		if ( distances[i] < distances[selectedLight] )
		{
			selectedLight = i;
		}
	}


	for ( int i = 0; i < MAX_LIGHTS; i++ )
	{
		if ( i == selectedLight )
		{
			// mouse pos to world pos
			if ( IsKeyDown( GLFW_KEY_SPACE ) )
			{
				float2 mousePos = renderer->mousePos;
				Ray primaryRay = renderer->camera.GetPrimaryRay( mousePos.x, mousePos.y );
				renderer->scene.FindNearest( primaryRay );
				float3 intersectionPoint = primaryRay.IntersectionPoint();

				float3 toNewPoint = intersectionPoint - renderer->lights[i]->GetPosition();
				float3 direction = normalize( toNewPoint );

				Ray moveRay( intersectionPoint, -direction, length( toNewPoint ) );


				if ( !renderer->scene.IsOccluded( moveRay ) && length( toNewPoint ) < 0.1f )
				{
					renderer->lights[i]->SetPosition( float3( intersectionPoint.x, renderer->lights[i]->GetPosition().y, intersectionPoint.z ) );
					renderer->ResetAccumulator();
				}
			}
		}


		if ( i >= renderer->lightCount || renderer->lights[i]->GetType() == Light::Directional || renderer->lights[i]->GetType() == Light::Spot ) break;

		for ( int j = 0; j < 10; j++ )
		{
			if ( gameObjects[j] == nullptr ) continue;
			Absorber* absorber = dynamic_cast<Absorber*>( gameObjects[j] );
			if ( absorber == nullptr ) continue;
			uint idx = absorber->GetPrimitiveIdx();
			switch ( renderer->scene.bvh->primitiveTypes[idx] )
			{
			case SPHERE:
			{


				float3 spherePos = renderer->scene.bvh->primitiveData[idx].sphere.pos;
				float3 sphereRadius = renderer->scene.bvh->primitiveData[idx].sphere.radius;
				float3 lightPos = renderer->lights[i]->GetPosition();

				float3 toLight = lightPos - spherePos;
				float3 direction = normalize( toLight );

				float3 up = float3( 0, 1, 0 );
				float3 right = cross( up, direction );
				up = cross( direction, right );

				float3 sphereEdges[3] = {
					spherePos - right * sphereRadius,
					spherePos + right * sphereRadius,
					spherePos + up * sphereRadius
				};
				for ( int k = 0; k < 3; k++ )
				{
					float3 intersectionPoint = sphereEdges[k];

					float3 toEdge = intersectionPoint - lightPos;
					float3 edgeDirection = normalize( toEdge );

					Ray occlusionRay( lightPos, edgeDirection, length( toEdge ) - 0.02f );

					//debug draw it

					bool occluded = renderer->scene.IsOccluded( occlusionRay );

					if ( !occluded )
					{
						float3 lightColor = renderer->lights[i]->GetColor();

						if ( lightColor.x != 0.0f ) lightColor.x = 1.0f;
						if ( lightColor.y != 0.0f ) lightColor.y = 1.0f;
						if ( lightColor.z != 0.0f ) lightColor.z = 1.0f;


						float chargingRate = absorber->GetRate();
						// if the absorber is red and the light is red it will increase redValue, if the absorber is red but the light is green it will decrease redValue
						if ( renderer->scene.bvh->primitiveData[idx].sphere.color == float3( 1.0f, 0.3f, 0.3f ) )
						{
							if ( lightColor == float3( 1.0f, 0.0f, 0.0f ) )
							{
								if ( redValue < maxValue )
								{
									redValue += chargingRate * deltaTime;
								}
							}
							else
							{
								if ( redValue > 0.0f )
								{
									redValue -= drainingRate * deltaTime;
								}

							}
						}
						else if ( renderer->scene.bvh->primitiveData[idx].sphere.color == float3( 0.3f, 1.0f, 0.3f ) )
						{
							if ( lightColor == float3( 0.0f, 1.0f, 0.0f ) )
							{
								if ( greenValue < maxValue )
								{
									greenValue += chargingRate * deltaTime;
								}
							}
							else
							{
								if ( greenValue > 0.0f )
								{
									greenValue -= drainingRate * deltaTime;
								}
							}
						}
						else if ( renderer->scene.bvh->primitiveData[idx].sphere.color == float3( 0.3f, 0.3f, 1.0f ) )
						{
							if ( lightColor == float3( 0.0f, 0.0f, 1.0f ) )
							{
								if ( blueValue < maxValue )
								{
									blueValue += chargingRate * deltaTime;
								}
							}
							else
							{
								if ( blueValue > 0.0f )
								{
									blueValue -= drainingRate * deltaTime;
								}
							}
						}
					}
				}
			}
			}
		}

	}
	DrawPointHUD();

	float3 distanceFromGoal( abs( maxValue - redValue ), abs( maxValue - greenValue ), abs( maxValue - blueValue ) );
	if ( length( distanceFromGoal ) < 5.0f )
	{
		//win
		win = true;
		playing = false;
		LoadLevel( 1 );
	}

}

void Game::Shutdown()
{
}

void Game::DrawPointHUD() const
{
	//background
	renderer->screen->Bar( 0, 0, 100, 40, 0xb9bec7 );

	//bars
	float redPercentage = redValue / maxValue;
	if ( redPercentage > 1.0f ) redPercentage = 1.0f;
	renderer->screen->Bar( 0, 0, static_cast<int>( 100.0f * redPercentage ), 10, 0xFF0000 );

	float greenPercentage = greenValue / maxValue;
	if ( greenPercentage > 1.0f ) greenPercentage = 1.0f;
	renderer->screen->Bar( 0, 15, static_cast <int>( 100.0f * greenPercentage ), 25, 0x00FF00 );

	float bluePercentage = blueValue / maxValue;
	if ( bluePercentage > 1.0f ) bluePercentage = 1.0f;
	renderer->screen->Bar( 0, 30, static_cast <int>( 100.0f * bluePercentage ), 40, 0x0000FF );
}

void Game::LoadLevel( const uint level )
{
	switch ( level )
	{
	case 0:
		for ( int i = 1; i < 23; i++ )
		{
			for ( int j = 3; j < 25; j++ )
			{
				renderer->scene.Set( i, 1, j, float3( 0.0f ), 0 );
			}
		}
		for ( int i = 0; i < 5; i++ )
		{
			renderer->scene.Set( 3 + i, 1, 11, float3( 1.0f ), 0 );
		}
		for ( int i = 0; i < 6; i++ )
		{
			renderer->scene.Set( 3, 1, 11 + i, float3( 1.0f ), 0 );
		}
		for ( int i = 0; i < 5; i++ )
		{
			renderer->scene.Set( 16 + i, 1, 4, float3( 1.0f ), 0 );
		}
		for ( int i = 0; i < 3; i++ )
		{
			renderer->scene.Set( 16, 1, 4 + i, float3( 1.0f ), 0 );
		}
		for ( int i = 0; i < 7; i++ )
		{
			renderer->scene.Set( 16 - i, 1, 7, float3( 1.0f ), 0 );
		}
		for ( int i = 0; i < 11; i++ )
		{
			renderer->scene.Set( 16, 1, 14 + i, float3( 1.0f ), 0 );
		}
		renderer->ResetAccumulator();
		break;
	case 1:
		for ( int i = 1; i < 23; i++ )
		{
			for ( int j = 3; j < 25; j++ )
			{
				renderer->scene.Set( i, 1, j, float3( 0.0f ), 0 );
			}
		}
		renderer->camera.camPos = float3(0.07f, 0.06f, 0.71f);
		renderer->camera.camTarget = float3(0.13f, 0.06f, 0.65f);
		renderer->scene.bvh->primitiveData[0].sphere.material = 2;
		renderer->scene.bvh->primitiveData[1].sphere.material = 3;
		renderer->scene.bvh->primitiveData[2].sphere.material = 2;

		renderer->scene.bvh->primitiveData[1].sphere.color = float3( 1.0f );

		renderer->scene.bvh->SetPosition( 0, float3( 0.16f, 0.055f, 0.51f ) );
		renderer->scene.bvh->SetPosition( 1, renderer->scene.bvh->primitiveData[1].sphere.pos + float3(-0.1f, 0.005f, 0.0f) );
		renderer->scene.bvh->SetPosition( 2, float3( 0.13f, 0.055f, 0.45f ) );
		for ( int i = 0; i < 4; i++ )
		{
			if ( i == 3 )
			{
				renderer->lights[i]->SetPosition( float3( 0.37f, 0.06f, 0.60f ) );
			}
			else
			{
				renderer->lights[i]->SetPosition( float3( 1.0f ) );
			}
		}
		testScene = true;
		renderer->ResetAccumulator();
		break;
	default:
		printf("Invalid level");
	}
	
}
