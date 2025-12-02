#include "precomp.h"
#include "skybox.h"
#include "game.h"
#include "material.h"
#include "light.h"

void Renderer::Init()
{
	newAccumulator = (float4*)MALLOC64( SCRWIDTH * SCRHEIGHT * 16 );
	memset( newAccumulator, 0, SCRWIDTH * SCRHEIGHT * 16 );

	oldAccumulator = (float4*)MALLOC64( SCRWIDTH * SCRHEIGHT * 16 );
	memset( oldAccumulator, 0, SCRWIDTH * SCRHEIGHT * 16 );

	skybox = new Skybox( "assets/quarry_cloudy_4k.hdr" ); // https://polyhaven.com/a/quarry_cloudy

	//camera.Initialize( float3( 0.43757f, 0.923119f, 0.440113 ), float3( 0.43757f, -0.076881, 0.440113));

	// load materials
	materials[0] = new Material(); // just normal elbedo
	materials[1] = new Material( "assets/white.png" ); // just normal with texture
	materials[2] = new Material(); // mirror
	materials[3] = new Material(); // glass

	materials[2]->reflectivity = 1.0f;

	materials[3]->ior = 1.2f;
	materials[3]->refractivity = 1.0f;

	Light pointLight;
	pointLight.SetType( Light::Type::Point );
	pointLight.SetPosition( float3( 0.2f, 0.0325f, 0.3f ) );
	pointLight.SetColor( float3( 0.1f, 0.0f, 0.0f ) );
	lights[lightCount++] = new Light( pointLight );

	Light pointLight2;
	pointLight2.SetType( Light::Type::Point );
	pointLight2.SetPosition( float3( 0.4f, 0.0325f, 0.5f ) );
	pointLight2.SetColor( float3( 0.0f, 0.1f, 0.0f ) );
	lights[lightCount++] = new Light( pointLight2 );

	Light pointLight3;
	pointLight3.SetType( Light::Type::Point );
	pointLight3.SetPosition( float3( 0.6f, 0.0325f, 0.23f ) );
	pointLight3.SetColor( float3( 0.0f, 0.0f, 0.1f ) );
	lights[lightCount++] = new Light( pointLight3 );

	Light spotLight;
	spotLight.SetType( Light::Type::Spot );
	spotLight.SetPosition(1.0f);
	spotLight.SetDirection( float3( 1.0f, 0.0f, 0.0f ) );
	spotLight.SetColor( float3( 0.0f, 0.0f, 0.1f ) );
	spotLight.SetAngle( 1.2f );
	lights[lightCount++] = new Light( spotLight );

	Light directionalLight;
	directionalLight.SetType( Light::Type::Directional );
	directionalLight.SetDirection( 1.0f );
	directionalLight.SetColor( float3( 1.0f ) );
	lights[lightCount++] = new Light( directionalLight );

	game = new Game();
	game->renderer = this;
	game->Init();
}

void Renderer::Tick( const float deltaTime )
{
	Timer t;
#pragma omp parallel for schedule(dynamic)
	for ( int y = 0; y < SCRHEIGHT; y++ )
	{
		for ( int x = 0; x < SCRWIDTH; x++ )
		{
			float4 pixelColor = float4( 0.0f );

			Ray primaryRay = camera.GetPrimaryRay( static_cast<float>( x ), static_cast<float>( y ) );
			float4 newColor = float4( Trace( primaryRay, 0 ), primaryRay.t );

			int index = x + y * SCRWIDTH;


			if ( primaryRay.t == 1e34f )
			{
				pixelColor = newColor;
				goto set; // lol goto, Jacco used this idk why i'm using it.
			}

			if ( reprojection )
			{
				float3 intersectionPoint = primaryRay.IntersectionPoint();
				int2 oldColorPos = camera.PyramidSpace( float4( intersectionPoint, 1.0f ) );
				float4 oldColor = oldAccumulator[oldColorPos.x + oldColorPos.y * SCRWIDTH];

				float oldDepth = oldColor.w;
				float depthDifference = abs( primaryRay.t - oldDepth );
				float adjustedAlpha = max( 0.0f, reprojectAlpha - depthDifference * 16 );

				// Vavg = aVprev + (1-a)Vnew
				pixelColor = adjustedAlpha * oldColor + ( 1 - adjustedAlpha ) * newColor;
				pixelColor.w = primaryRay.t;
			}
			else
			{
				pixelColor = newColor;
			}
		set:
			newAccumulator[index] = pixelColor;
			if ( toneMapping ) pixelColor = float4( AcesApproximated( float3( pixelColor.x, pixelColor.y, pixelColor.z ) ), primaryRay.t );
			screen->pixels[index] = RGBF32_to_RGB8( &pixelColor );
		}
	}
	memcpy( oldAccumulator, newAccumulator, SCRWIDTH * SCRHEIGHT * 16 );


	reprojectAlpha += 0.00005f * deltaTime;
	if ( reprojectAlpha > 0.999f ) reprojectAlpha = 0.999f;

	// handle user input
	if ( IsKeyPressed( GLFW_KEY_R ) )
	{
		//debugLines.clear();
		Ray ray = camera.GetPrimaryRay( static_cast<float>( mousePos.x ), static_cast<float>( mousePos.y ) );
		Trace( ray, 0, true );
		debugLines.push_back( new DebugLine( ray.origin, ray.IntersectionPoint(), 0xFF0000 ) );
	}

	if ( IsKeyPressed( GLFW_KEY_T ) ) debugLines.clear();


	if ( IsKeyPressed( GLFW_KEY_Y ) )
	{
		Ray ray = camera.GetPrimaryRay( static_cast<float>( mousePos.x ), static_cast<float>( mousePos.y ) );
		scene.FindNearest( ray );
		if ( ray.t != 1e34f )
		{
			activeMaterial = ray.GetMaterial();
		}
	}

	if ( IsKeyPressed( GLFW_KEY_U ) )
	{
		Ray ray = camera.GetPrimaryRay( static_cast<float>( mousePos.x ), static_cast<float>( mousePos.y ) );
		float3 result = Trace( ray, maxRayDepth );
		cout << result.x << " " << result.y << " " << result.z << endl;
	}

	if ( IsKeyPressed( GLFW_KEY_I ) )
	{
		debug = !debug;
	}

	for ( DebugLine* line : debugLines )
	{
		float2 screenPos1, screenPos2 = float2( 0.0f );
		bool succeeded = camera.WorldSpaceToScreenSpace( line->start, screenPos1 ) && camera.WorldSpaceToScreenSpace( line->end, screenPos2 );

		if ( succeeded )screen->Line( screenPos1.x, screenPos1.y, screenPos2.x, screenPos2.y, line->color );
	}
	game->Update( deltaTime );

	camera.CalculateFrustum(); // for reprojection
	if ( camera.HandleInput( deltaTime ) ) reprojectAlpha = 0.9f;

	avg = ( 1 - alpha ) * avg + alpha * t.elapsed() * 1000;
	if ( alpha > 0.05f ) alpha *= 0.5f;

	float deltaTimeSeconds = deltaTime / 1000.0f;
	fps = 1.0f / deltaTimeSeconds;
}

void Renderer::UI()
{
	if ( debug )
	{
		ImGui::Begin( "RuhTracer" );
		ImGui::BeginChild( "Scroll" );
		PerformanceUI();
		ImGui::NewLine();
		SettingsUI();
		ImGui::EndChild();
		ImGui::End();
	}

	if ( !game->playing )
	{
		int2 pauseWindowSize = int2( 400, 200 );
		if ( game->win )
		{
			ImGui::SetNextWindowPos( ImVec2( static_cast<float>( SCRWIDTH ) - pauseWindowSize.x / 2, static_cast<float>( SCRHEIGHT ) - pauseWindowSize.y / 2 - 240) );
		}
		else
		{
			ImGui::SetNextWindowPos( ImVec2( static_cast<float>( SCRWIDTH ) - pauseWindowSize.x / 2, static_cast<float>( SCRHEIGHT ) - pauseWindowSize.y / 2 ) );
		}
		ImGui::SetNextWindowSize( ImVec2( static_cast<float>( pauseWindowSize.x ), static_cast<float>( pauseWindowSize.y ) ) );
		ImGui::Begin( "Game", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove );
		const char* text = game->win ? "You won!" : "Absorb";
		ImVec2 textSize = ImGui::CalcTextSize( text );
		ImVec2 textPosition( ( ImGui::GetWindowSize().x - textSize.x ) * 0.5f, ( ImGui::GetWindowSize().y - textSize.y ) * 0.5f );
		ImGui::SetCursorPos( textPosition );
		ImGui::Text( text );

		ImVec2 buttonSize( 100, 25 );
		ImVec2 buttonPosition( ( ImGui::GetWindowSize().x - buttonSize.x ) * 0.5f, ( ImGui::GetWindowSize().y - buttonSize.y ) * 0.5f + 50 );
		ImGui::SetCursorPos( buttonPosition );
		if ( ImGui::Button( game->win ? "Quit" : "Play", buttonSize ) )
		{
			if ( !game->win )
			{
				game->playing = true;
			}
			else
			{
				shouldQuit = true;
			}
		}

		text = "Check the README.txt for the controls";
		textSize = ImGui::CalcTextSize( text );
		textPosition = ImVec2( ( ImGui::GetWindowSize().x - textSize.x ) * 0.5f, ( ImGui::GetWindowSize().y - textSize.y ) * 0.5f + 90);
		ImGui::SetCursorPos( textPosition );
		ImGui::Text( text );
		ImGui::End();
	}
}

void Renderer::Shutdown()
{
	game->Shutdown();
	delete skybox;
	skybox = nullptr;
	delete game;
	game = nullptr;

	for ( Light* light : lights )
	{
		delete light;
		light = nullptr;
	}
	for ( DebugLine* line : debugLines )
	{
		delete line;
		line = nullptr;
	}
	for ( int i = 0; i < MATERIALS; i++ )
	{
		delete materials[i];
		materials[i] = nullptr;
	}
}

float3 Renderer::Trace( Ray& ray, int depth, const bool draw )
{
	if ( depth > maxRayDepth ) return float3( 0.0f ); // Terminate recursion

	scene.FindNearest( ray );
	if ( ray.t == 1e34f ) return skybox->render( ray ); // if no voxel hit, return skybox color

	Material* material = materials[ray.GetMaterial()];
	float3 normal = ray.N;
	float3 albedo = ray.GetColor();

	if ( material->texture )
	{
		float3 newColor;
		uint pixel;
		switch ( ray.type )
		{
		case 0: //sphere
		{
			Sphere& sphere = scene.bvh->primitiveData[ray.sphereIdx].sphere;
			float3 sphereCenter = sphere.pos;
			float3 point = normalize( ray.IntersectionPoint() - sphereCenter );

			float theta = atan2( point.z, point.x );
			float phi = acos( point.y );

			float u = ( theta + PI ) / ( 2 * PI );
			float v = phi / PI;

			pixel = material->texture->PixelAt( static_cast<uint>( u * material->texture->width ), static_cast<uint>( v * material->texture->height ) );

			newColor.x = ( ( pixel >> 16 ) & 0xFF ) / 255.0f;
			newColor.y = ( ( pixel >> 8 ) & 0xFF ) / 255.0f;
			newColor.z = ( pixel & 0xFF ) / 255.0f;

			albedo *= newColor;
			break;
		}
		case 1: // cube/voxel
		{
			float2 uv = ray.GetUV();
			int x = static_cast<int>( uv.x * material->texture->width );
			int y = static_cast<int>( uv.y * material->texture->height );

			pixel = material->texture->PixelAt( x, y );

			newColor.x = ( ( pixel >> 16 ) & 0xFF ) / 255.0f;
			newColor.y = ( ( pixel >> 8 ) & 0xFF ) / 255.0f;
			newColor.z = ( pixel & 0xFF ) / 255.0f;

			albedo *= newColor;
			break;
		}
		}

	}

	float3 outColor = float3( 0.0f );

	float reflectivity = material->reflectivity;
	float refractivity = material->refractivity;
	float roughness = material->roughness;
	float diffuse = 1 - ( reflectivity + refractivity );

	float3 intersectionPoint = ray.IntersectionPoint();

	if ( reflectivity > 0 )
	{
		float3 reflection = reflect( ray.direction, normal );

		float3 roughReflection = mix( reflection, RandomDirectionOnHemisphere( normal ), roughness );

		Ray reflectedRay( intersectionPoint + roughReflection * EPSILON, roughReflection );

		if ( material->reflectivity > RandomFloat() )
			outColor += albedo * Trace( reflectedRay, depth + 1, draw );
		else
			diffuse = 1 - refractivity;

		if ( draw )
		{
			debugLines.push_back( new DebugLine( intersectionPoint, intersectionPoint + roughReflection * ( reflectedRay.t == 1e34f ? 0.5f : reflectedRay.t ), 0x00FF00 ) );
			debugLines.push_back( new DebugLine( intersectionPoint, intersectionPoint + ray.N * 0.025f, 0xFFFF00 ) );
		}
	}

	if ( refractivity > 0 )
	{
		bool tir;
		float3 refractedDirection;
		tir = refract( ray.direction, normal, material->ior, refractedDirection );
		Ray refractedRay( intersectionPoint + refractedDirection * EPSILON, refractedDirection );
		if ( !tir )
		{
			outColor += albedo * Trace( refractedRay, depth + 1, draw );
			if ( draw )
			{
				debugLines.push_back( new DebugLine( refractedRay.origin, intersectionPoint + refractedDirection * ( refractedRay.t == 1e34f ? 0.25f : refractedRay.t ), refractedRay.t == 1e34f ? 0x00FF00 : 0x0000FF ) );
				debugLines.push_back( new DebugLine( intersectionPoint, intersectionPoint + ray.N * 0.025f, 0xFFFF00 ) );
			}
		}
	}

	if ( diffuse > 0 )
	{
		float3 directLight = ComputeDirectLighting( ray, normal );

		float3 indirectLight = float3( 0.0f );

		float3 randomDirection = RandomDirectionOnHemisphere( normal );
		Ray indirectRay( ray.IntersectionPoint() + randomDirection * EPSILON, randomDirection );
		indirectLight = Trace( indirectRay, depth + 1, draw ) + ( albedo * material->emission );

		float3 tempSolve = outColor + diffuse * albedo * ( directLight + indirectLight );

		//https://github.com/jbikker/lighthouse2/blob/e61e65444d8ed3074775003f7aa7d60cb0d4792e/lib/rendercore_vulkan_rt/shaders/tools.glsl#L169
		float survivalProbability = min( 1.0f, max( max( tempSolve.x, tempSolve.y ), tempSolve.z ) ); // i think this is pretty fucky on performance but idk
		if ( RandomFloat() < survivalProbability )
		{
			outColor = tempSolve / survivalProbability;
		}
	}

	return outColor;
}

float3 Renderer::ComputeDirectLighting( const Ray& ray, const float3& normal )
{
	bool notEmpty = false;
	for ( int i = 0; i < MAX_LIGHTS; i++ ) { if ( lights[i] != nullptr ) { notEmpty = true; break; } }
	if ( !notEmpty ) return float3( 0.0f );

	float3 result = float3( 0.0f );

	float3 intersectionPoint = ray.IntersectionPoint();

	float weights[MAX_LIGHTS];
	int lightIndex = PickLight( intersectionPoint, weights );
	Light* light = lights[lightIndex];

	float3 lightPosition = light->GetPosition();
	float3 lightDirection = light->GetDirection();
	float3 lightColor = light->GetColor();


	switch ( light->GetType() )
	{
	case Light::Type::Directional:
	{
		const float3 direction = normalize( lightDirection );
		const float cosTheta = dot( normal, direction );

		Ray shadowRay( intersectionPoint + direction * EPSILON, direction - EPSILON );
		if ( !scene.IsOccluded( shadowRay ) )
		{
			result += max( 0.0f, cosTheta ) * ( lightColor / weights[lightIndex] );
		}
		break;
	}
	case Light::Type::Point:
	{

		const float3 toLight = lightPosition - intersectionPoint;
		const float3 direction = normalize( toLight );
		const float distance = length( toLight );

		const float cosTheta = dot( normal, direction );
		if ( cosTheta < EPSILON ) return float3( 0.0f );

		const float attenuation = 1.0f / ( distance * distance );

		// shadows
		Ray shadowRay( intersectionPoint + direction * EPSILON, direction, distance - EPSILON );
		if ( !scene.IsOccluded( shadowRay ) )
		{
			result += ( lightColor / weights[lightIndex] ) * attenuation * cosTheta;
		}
		break;
	}
	case Light::Type::Spot:
	{
		const float3 toLight = lightPosition - intersectionPoint;
		const float3 direction = normalize( toLight );

		const float cosSpotAngle = dot( normalize( lightDirection ), -direction );

		const float cosAngle = cos( light->GetAngle() * 0.5f );

		if ( cosSpotAngle > cosAngle )
		{
			const float distance = length( lightPosition - intersectionPoint );

			const float attenuation = 1.0f / ( distance * distance );

			const float falloff = smoothstep( cosAngle, 1.0f, cosSpotAngle );

			Ray shadowRay( intersectionPoint + direction * EPSILON, direction, distance - EPSILON );
			if ( !scene.IsOccluded( shadowRay ) )
			{
				result += max( 0.0f, falloff ) * attenuation * ( lightColor / weights[lightIndex] );
			}
		}
		break;
	}
	}
	return result;
}

int Renderer::PickLight( const float3& intersectionPoint, float* weights ) const
{
	float invDistances[MAX_LIGHTS];
	float brightnesses[MAX_LIGHTS];

	float totalBrightness = 0.0f;
	float totalInvDistance = 0.0f;


	for ( size_t i = 0; i < MAX_LIGHTS; i++ )
	{
		if ( lights[i] == nullptr ) break;
		float3 lightPosition = lights[i]->GetPosition();
		float3 toLight = lightPosition - intersectionPoint;
		invDistances[i] = 1.0f / length( toLight );

		brightnesses[i] = length( lights[i]->GetColor() );

		totalBrightness += brightnesses[i];
		totalInvDistance += invDistances[i];

	}

	float totalWeight = 0.0f;
	for ( int i = 0; i < MAX_LIGHTS; i++ )
	{
		if ( lights[i] == nullptr ) break;
		float brightnessWeight = brightnesses[i] / totalBrightness;
		float distanceWeight = invDistances[i] / totalInvDistance;
		float weight = brightnessWeight * distanceWeight;
		weights[i] = weight;
		totalWeight += weights[i];
	}

	float accumulatedWeight = 0.0f;
	for ( int i = 0; i < MAX_LIGHTS; i++ )
	{
		if ( lights[i] == nullptr ) break;
		weights[i] /= totalWeight;

		accumulatedWeight += weights[i];
		if ( RandomFloat() < accumulatedWeight )
		{
			return i;
		}
	}
	return 0;
}


float3 Renderer::RandomDirectionOnHemisphere( const float3& normal ) const
{
	float theta = 2.0f * PI * RandomFloat();
	float phi = acos( sqrt( RandomFloat() ) );

	float3 randomDirection;
	randomDirection.x = sin( phi ) * cos( theta );
	randomDirection.y = sin( phi ) * sin( theta );
	randomDirection.z = cos( phi );

	if ( dot( randomDirection, normal ) < 0.0f )
		randomDirection = -randomDirection;

	return randomDirection;
}

float3 Renderer::AcesApproximated( float3 C ) const
{
	C = C * 0.6f;
	float a = 2.51f;
	float b = 0.03f;
	float c = 2.43f;
	float d = 0.59f;
	float e = 0.14f;
	C = ( C * ( a * C + b ) ) / ( C * ( c * C + d ) + e );
	return clamp( C, 0.0f, 1.0f );
}

void Renderer::AcesApproximatedSSE( const float3* C, const int numC, float3* result )
{
	__m128 a = _mm_set1_ps( 2.51f );
	__m128 b = _mm_set1_ps( 0.03f );
	__m128 c = _mm_set1_ps( 2.43f );
	__m128 d = _mm_set1_ps( 0.59f );
	__m128 e = _mm_set1_ps( 0.14f );
	__m128 multiplier = _mm_set1_ps( 0.6f );

	for ( int i = 0; i < numC; ++i )
	{
		__m128 C_simd = _mm_set_ps( 0, C[i].z, C[i].y, C[i].x );
		C_simd = _mm_mul_ps( C_simd, multiplier );

		// incredible naming scheme, I know
		__m128 aTimesCPlusB = _mm_add_ps( _mm_mul_ps( a, C_simd ), b );
		__m128 cTimesCPlusD = _mm_add_ps( _mm_mul_ps( c, C_simd ), d );

		__m128 cTimesATimesCPlusB = _mm_mul_ps( C_simd, aTimesCPlusB );
		__m128 cTimesCTimesCPlusDPlusE = _mm_add_ps( _mm_mul_ps( C_simd, cTimesCPlusD ), e );

		__m128 division = _mm_div_ps( cTimesATimesCPlusB, cTimesCTimesCPlusDPlusE );
		// clamp
		division = _mm_max_ps( _mm_min_ps( division, _mm_set1_ps( 1.0f ) ), _mm_setzero_ps() );
		// funky conversion lol
		_mm_store_ps( reinterpret_cast<float*>( &result[i] ), division );
	}
}

void Renderer::PerformanceUI() const
{
	ImGui::Text( "fps: %.1ffps", fps );
	ImGui::Text( "ray time: %5.2fms", avg );
	ImGui::Text( "%f %f %f", camera.camPos.x, camera.camPos.y, camera.camPos.z );
	ImGui::Text( "%i %i %i", static_cast<int>( camera.camPos.x * WORLDSIZE ), static_cast<int>( camera.camPos.y * WORLDSIZE ), static_cast<int>( camera.camPos.z * WORLDSIZE ) );
	ImGui::Text( "x%i y%i", mousePos.x, mousePos.y );
	ImGui::Text( "camera target: %f %f %f", camera.camTarget.x, camera.camTarget.y, camera.camTarget.z );
}

void Renderer::SettingsUI()
{
	// Lighting settings
	ImGui::Text( "Lighting settings" );
	ImGui::Separator();
	ImGui::Indent( 16.0f );
	ImGui::SliderInt( "maxRayDepth", &maxRayDepth, 0, 100 );
	ImGui::Unindent( 16.0f );


	// Material settings
	ImGui::Text( "Material settings" );
	ImGui::Separator();
	ImGui::Indent( 16.0f );
	ImGui::Text( "Active material: %i", activeMaterial );
	ImGui::SliderFloat( "reflectivity", &materials[activeMaterial]->reflectivity, 0.f, 1.f );
	ImGui::SliderFloat( "refractivity", &materials[activeMaterial]->refractivity, 0.f, 1.f );
	ImGui::SliderFloat( "roughness", &materials[activeMaterial]->roughness, 0.f, 1.f );
	ImGui::SliderFloat( "ior", &materials[activeMaterial]->ior, 1.f, 2.0f );
	ImGui::SliderFloat( "emission", &materials[activeMaterial]->emission, 0.0f, 25.0f );
	ImGui::Unindent( 16.0f );

	// Misc settings
	ImGui::Text( "Misc settings" );
	ImGui::Separator();
	ImGui::Indent( 16.0f );
	ImGui::Checkbox( "toneMapping", &toneMapping );
	ImGui::Checkbox( "reprojection", &reprojection );
	ImGui::SliderFloat( "reprojectAlpha", &reprojectAlpha, 0.f, 1.f );

	ImGui::Unindent( 16.0f );

	// Spawning
	ImGui::Text( "Spawning" );
	ImGui::Separator();
	ImGui::Indent( 16.0f );
	if ( ImGui::Button( "Point light" ) )
	{
		if ( lightCount >= MAX_LIGHTS ) return;
		Light* pointLight = new Light;
		pointLight->SetType( Light::Type::Point );
		pointLight->SetPosition( camera.camPos );
		pointLight->SetColor( float3( 0.03f, 0.03f, 0.03f ) );
		lights[lightCount++] = pointLight;
	}
	ImGui::SameLine();
	if ( ImGui::Button( "Spot light" ) )
	{
		if ( lightCount >= MAX_LIGHTS ) return;
		Light* spotLight = new Light;
		spotLight->SetType( Light::Type::Spot );
		spotLight->SetPosition( camera.camPos );
		spotLight->SetColor( float3( 0.03f, 0.03f, 0.03f ) );
		spotLight->SetDirection( float3( 0.0f, 0.0f, -1.0f ) );
		spotLight->SetAngle( 30 * PI / 180.0 );
		lights[lightCount++] = spotLight;
	}
	ImGui::Unindent( 16.0f );

	//Lights
	ImGui::Text( "Lights" );
	ImGui::Separator();
	ImGui::Indent( 16.0f );
	for ( int i = static_cast<int>( MAX_LIGHTS ) - 1; i >= 0; i-- )
	{
		if ( lights[i] == nullptr ) continue;
		Light::Type type = lights[i]->GetType();

		// light title
		if ( type == Light::Type::Point )
			ImGui::Text( "Point light" );
		else if ( type == Light::Type::Spot )
			ImGui::Text( "Spot light" );
		else if ( type == Light::Type::Directional )
			ImGui::Text( "Directional light" );

		// delete
		if ( ImGui::Button( ( "Delete" + std::to_string( i ) ).c_str() ) )
		{
			delete lights[i];
			lights[i] = nullptr;
			for ( int j = i; j < MAX_LIGHTS - 1; ++j )
			{
				lights[j] = lights[j + 1];
			}

			--lightCount;
			ResetAccumulator();
			return;
		}

		// set pos to camera
		if ( ImGui::Button( ( "Set to camera" + std::to_string( i ) ).c_str() ) )
		{
			lights[i]->SetPosition( camera.camPos );
		}

		//pos input fields
		float3 position = lights[i]->GetPosition();
		if ( ImGui::InputFloat3( ( "Position" + std::to_string( i ) ).c_str(), &position.x ) )
		{
			lights[i]->SetPosition( position );
		}

		// color
		float3 color = lights[i]->GetColor();
		if ( ImGui::ColorEdit3( ( "Color" + std::to_string( i ) ).c_str(), &color.x ) )
		{
			lights[i]->SetColor( color );
		}

		// direction
		if ( type == Light::Type::Spot || type == Light::Type::Directional )
		{
			float3 direction = lights[i]->GetDirection();
			if ( ImGui::InputFloat3( ( "Direction" + std::to_string( i ) ).c_str(), &direction.x ) )
			{
				lights[i]->SetDirection( direction );
			}
		}

		// angle
		if ( type == Light::Type::Spot )
		{
			float angle = lights[i]->GetAngle();
			if ( ImGui::SliderFloat( ( "Angle" + std::to_string( i ) ).c_str(), &angle, 0.0f, PI ) )
			{
				lights[i]->SetAngle( angle );
			}
		}
		ImGui::Separator();
	}
	ImGui::Unindent( 16.0f );
}