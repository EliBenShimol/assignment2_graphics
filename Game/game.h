#pragma once
#include "scene.h"

class Game : public Scene
{
public:
	
	Game();
	Game(float angle,float relationWH,float near, float far);
	void Init();
	void Update(const glm::mat4 &MVP,const glm::mat4 &Model,const int  shaderIndx);
	void ControlPointUpdate();
	void WhenRotate();
	void WhenTranslate();
	void Motion();
	bool hitSphere(glm::vec3 posPixel, glm::vec4 sphere, glm::vec3 vector, float ans[]);
	float magnitude(glm::vec3 vector);
	void solveForT(float* parts, float* ans);
	bool isHitBefore(glm::vec3 posVector, glm::vec4 objects, glm::vec3 pos, int isRef, int isTrans);
	float my_dist(glm::vec3 a, glm::vec3 b);
	float solveForP(glm::vec4 plane, glm::vec3 vec, glm::vec3 pos);
	glm::vec3 calcColor(glm::vec3 amb, int countLights, glm::vec4* lights, glm::vec4 obj, glm::vec3 hitPos, glm::vec3* colorLights,
		glm::vec4* spotLight, glm::vec4 objColor, glm::vec3 tmpVec, int countObj, float dist, glm::vec4* objects,
		glm::vec3 pos, glm::vec4* objectsColors, int isRef, int isTrans, int refPlace[], int transPlace[], glm::vec3 ambient,
		int countRef, glm::vec3 norm);
	void printVec(glm::vec3 a);
	void printVec(glm::vec4 a);
	glm::vec3 calcRef();
	~Game(void);
};

