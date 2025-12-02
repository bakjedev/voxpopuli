#pragma once

namespace Tmpl8
{
	class Renderer;
	class GameObject;


class Game
{
public:
	Game() = default;
	~Game() = default;

	void Init();
	void Update(const float deltaTime);
	void Shutdown();

	void LoadLevel(const uint level);

	Renderer* renderer;
	bool playing = false;
	bool win = false;
private:
	void DrawPointHUD() const;

	bool testScene = false;

	int selectedLight = 0;

	GameObject* gameObjects[10] = { nullptr };

	float maxValue = 1000.0f;
	float drainingRate = 0.4f;

	float redValue = 0.0f;
	float greenValue = 0.0f;
	float blueValue = 0.0f;
};

}