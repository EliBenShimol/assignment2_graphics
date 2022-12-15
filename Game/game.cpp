#include "game.h"
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include <stb_image.h>
#include <fstream>
using namespace glm;
using namespace std;

static void printMat(const glm::mat4 mat)
{
	std::cout << " matrix:" << std::endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout << mat[j][i] << " ";
		std::cout << std::endl;
	}
}

Game::Game() : Scene()
{
}

Game::Game(float angle, float relationWH, float near1, float far1) : Scene(angle, relationWH, near1, far1)
{
}

void Game::Init()
{

	AddShader("../res/shaders/pickingShader");
	AddShader("../res/shaders/basicShader");

	int width, height, numComponents;
	std::string fileName = "../res/textures/plane.png";
	//unsigned char* data = stbi_load((fileName).c_str(), &width, &height, &numComponents, 4);
	unsigned char* data = (unsigned char*)malloc(4*800*800);
	for (int row = 0; row < 800; row++) {
		for (int count = 0; count < 4 * 800; count = count + 4) {
			data[row * 4 * 800 + count] = 0.0f;
			data[row * 4 * 800 + count + 1] = 0.0f;
			data[row * 4 * 800 + count + 2] = 0.0f;
			data[row * 4 * 800 + count + 3] = 1.0f;
		}
	}
	char readData[100];
	char* path = "../res/textures/scene5.txt";
	ifstream infile;
	infile.open(path);
	//getting cam position
	glm::vec3 cam = glm::vec3(0.0, 0.0, 0.0);
	glm::vec3 specular = glm::vec3(0.7, 0.7, 0.7);
	glm::vec3 ambient = glm::vec3(0.0, 0.0, 0.0);
	infile >> readData;
	for (int i = 0; i < 3; i++) {
		infile >> readData;
		cam[i] = stof(readData);
	}
	infile >> readData;

	//getting ambient
	infile >> readData;
	for (int i = 0; i < 3; i++) {
		infile >> readData;
		ambient[i] = stof(readData);
	}
	infile >> readData;

	//getting all the objects
	infile >> readData;
	int countObj = 0;
	vector <vec4> objectsTmp;
	vector <int> refPlace;
	vector <int> transPlace;
	while (readData[0] != 'c') {
		if (readData[0] == 'r')
			refPlace.push_back(countObj);
		else
			refPlace.push_back(-1);
		if (readData[0] == 't')
			transPlace.push_back(countObj);
		else
			transPlace.push_back(-1);
		countObj++;
		vec4 object = vec4(0.0, 0.0, 0.0, 0.0);
		for (int i = 0; i < 4; i++) {
			infile >> readData;
			object[i] = stof(readData);
		}
		//printVec(object);
		objectsTmp.push_back(object);
		infile >> readData;
	}

	glm::vec4* objects = new vec4[countObj];
	for (int i = 0; i < countObj; i++) {
		objects[i] = objectsTmp[i];
	}
	
	//objects colors
	glm::vec4* objectsColors = new vec4[countObj];
	int countTmp = 0;
	for (int i = 0; i < countObj; i++) {
		vec4 object = vec4(0.0, 0.0, 0.0, 0.0);
		for (int i = 0; i < 4; i++) {
			infile >> readData;
			object[i] = stof(readData);
		}
		objectsColors[i] = object;
		countTmp++;
		infile >> readData;
	}
	
	//lights
	vector <vec4> lightsTmp;
	int countLights = 0;
	while (readData[0] == 'd') {
		vec4 object = vec4(0.0, 0.0, 0.0, 0.0);
		for (int i = 0; i < 4; i++) {
			infile >> readData;
			object[i] = stof(readData);
		}
		lightsTmp.push_back(object);
		countLights++;
		infile >> readData;
	}
	glm::vec4* lights = new vec4[countLights];
	glm::vec4* spotLight = new vec4[countLights];
	for (int i = 0; i < countLights; i++) {
		lights[i] = lightsTmp[i];
		if (lights[i][3] == 1.0f) {
			vec4 object = vec4(0.0, 0.0, 0.0, 0.0);
			for (int i = 0; i < 4; i++) {
				infile >> readData;
				object[i] = stof(readData);
			}
			spotLight[i] = object;
			infile >> readData;
		}
	}
	glm::vec3* colorLights = new vec3[countLights];
	for (int i = 0; i < countLights; i++) {
		vec4 object = vec4(0.0, 0.0, 0.0, 0.0);
		for (int i = 0; i < 4; i++) {
			infile >> readData;
			object[i] = stof(readData);
		}
		colorLights[i] = vec3(object);
		infile >> readData;
	}
	double r = 1.0 / 400.0;
	vector<vector<double>>* pixelsPlace = new vector<vector<double>>();
	vector<vector<double>>* pixelsVector = new vector<vector<double>>();
	for (int row = 0; row < 800; row++) {
		vector<double> tmpPlace;
		vector<double> tmpVector;
		for (int count = 0; count < 800; count++) {
			double x = (count - 399) * r;
			double y = (row - 399) * r;
			double z = 0.0;
			glm::vec3 tmpVec = glm::vec3(x, y, z);
			tmpPlace.push_back(x);
			tmpPlace.push_back(y);
			tmpPlace.push_back(z);
			glm::vec3 normalVec = tmpVec - cam;
			normalVec = glm::normalize(normalVec);
			tmpVector.push_back(normalVec[0]);
			tmpVector.push_back(normalVec[1]);
			tmpVector.push_back(normalVec[2]);
		}
		pixelsPlace->push_back(tmpPlace);
		pixelsVector->push_back(tmpVector);
	}
	float ans[] = { -1.0, -1.0 };
	glm::vec3 pos = glm::vec3(0.0, 0.0, 0.0);
	glm::vec3 tmpVec = glm::vec3(0.0, 0.0, 0.0);
	int* tmpRef = new int[countObj];
	for (int i = 0; i < countObj; i++)
		tmpRef[i] = refPlace[i];
	int* tmpTrans = new int[countObj];
	for (int i = 0; i < countObj; i++)
		tmpTrans[i] = transPlace[i];
	for (int row = 0; row < 800; row++) {
		for (int count = 0; count < 800; count++) {
			ans[0] = -1.0;
			ans[1] = -1.0;
			pos[0] = (*pixelsPlace)[row][count * 3];
			pos[1] = (*pixelsPlace)[row][count * 3 + 1];
			pos[2] = (*pixelsPlace)[row][count * 3 + 2];
			tmpVec[0] = (*pixelsVector)[row][count * 3];
			tmpVec[1] = (*pixelsVector)[row][count * 3 + 1];
			tmpVec[2] = (*pixelsVector)[row][count * 3 + 2];
			glm::vec4 obj = objects[0];
			glm::vec4 objColor = objectsColors[0];
			glm::vec3 hitPos = pos;
			int isRef = -2;
			int isTrans = -2;
			float dist = -1;
			for (int i = 0; i < countObj; i++) {
				ans[0] = -1.0;
				ans[1] = -1.0;
				bool check = false;
				check = hitSphere(pos, objects[i], tmpVec, ans);
				if (check) {
					float tmpHit = ans[0];
					if (ans[1] < ans[0] && ans[1] > 0.0)
						tmpHit = ans[1];
					glm::vec3 tmpHitPos = pos + tmpVec * tmpHit;
					float tmpDist = my_dist(tmpHitPos, pos);
					if ((tmpDist < dist && tmpDist > 0) || (dist < 0 && tmpDist > 0)) {
						dist = tmpDist;
						obj = objects[i];
						objColor = objectsColors[i];
						hitPos = tmpHitPos;
						isRef = refPlace[i];
						isTrans = transPlace[i];
					}

				}

			}
			//printVec(glm::vec3(obj));
			if (dist > 0) {
				glm::vec3 norm = glm::vec3(obj);
				if (obj[3] > 0) {
					norm = glm::vec3(hitPos[0] - obj[0], hitPos[1] - obj[1], hitPos[2] - obj[2]);
				}
				norm = glm::normalize(norm);
				//ambient
				glm::vec3 amb = glm::vec3(objColor[0] * ambient[0], objColor[1] * ambient[1], objColor[2] * ambient[2]);
				vec3 color = calcColor(amb, countLights, lights, obj, hitPos, colorLights, spotLight, objColor, tmpVec, countObj,
					dist, objects, pos, objectsColors, isRef, isTrans, tmpRef, tmpTrans, ambient, 0, norm);
				if (color[0] > 1)
					color[0] = 1;
				if (color[1] > 1)
					color[1] = 1;
				if (color[2] > 1)
					color[2] = 1;
				data[row * 4 * 800 + 4*count] = color[0] * 255.0f;
				data[row * 4 * 800 + 4*count + 1] = color[1] * 255.0f;
				data[row * 4 * 800 + 4*count + 2] = color[2] * 255.0f;
			}
		}
	}
	AddTexture(800, 800, data);//greyScale image
	AddShape(Plane, -1, TRIANGLES);
	pickedShape = 0;
	
	SetShapeTex(0,0);
	MoveCamera(0,zTranslate,10);
	pickedShape = -1;
	
	//ReadPixel(); //uncomment when you are reading from the z-buffer
}

void Game::Update(const glm::mat4 &MVP,const glm::mat4 &Model,const int  shaderIndx)
{
	Shader *s = shaders[shaderIndx];
	int r = ((pickedShape+1) & 0x000000FF) >>  0;
	int g = ((pickedShape+1) & 0x0000FF00) >>  8;
	int b = ((pickedShape+1) & 0x00FF0000) >> 16;
	s->Bind();
	s->SetUniformMat4f("MVP", MVP);
	s->SetUniformMat4f("Normal",Model);
	s->SetUniform4f("lightDirection", 0.0f , 0.0f, -1.0f, 0.0f);
	if(shaderIndx == 0)
		s->SetUniform4f("lightColor",r/255.0f, g/255.0f, b/255.0f,1.0f);
	else 
		s->SetUniform4f("lightColor",0.7f,0.8f,0.1f,1.0f);
	s->Unbind();
}

void Game::WhenRotate()
{
}

void Game::WhenTranslate()
{
}

void Game::Motion()
{
	if(isActive)
	{
	}
}

bool Game::hitSphere(glm::vec3 posPixel, glm::vec4 sphere, glm::vec3 vector, float ans[])
{
	if (sphere[3] > 0) {
		float a = 1.0;
		glm::vec3 tmp = glm::vec3(posPixel[0] - sphere[0], posPixel[1] - sphere[1], posPixel[2] - sphere[2]);
		float b = 2 * dot(vector, tmp);
		float c = magnitude(tmp);
		c = c * c - sphere[3] * sphere[3];
		float parts[] = { a, b, c };
		solveForT(parts, ans);
		if (ans[0] <= 0.0 && ans[1] <= 0.0) {
			return false;
		}
		return true;
	}
	else {
		float ret = solveForP(sphere, vector, posPixel);
		if (ret > 0) {
			ans[0] = ret;
			return true;
		}
		return false;
	}
}

float Game::magnitude(glm::vec3 vector)
{
	float x = vector[0] * vector[0];
	float y = vector[1] * vector[1];
	float z = vector[2] * vector[2];
	float ans = x + y + z;
	ans = sqrt(ans);
	return ans;
}

void Game::solveForT(float* parts, float* ans)
{
	float part1 = parts[1] * parts[1];
	part1 = part1 - 4 * parts[0] * parts[2];
	if (part1 >= 0) {
		ans[0] = -1 * parts[1] + part1;
		ans[1] = -1 * parts[1] - part1;
		ans[0] = ans[0] / (2 * parts[0]);
		ans[1] = ans[1] / (2 * parts[0]);
	}
}

bool Game::isHitBefore(glm::vec3 posVector, glm::vec4 object, glm::vec3 pos, int isRef, int isTrans)
{
	float ans[] = { -1.0, -1.0 };
		if (hitSphere(pos, object, posVector, ans)) {
			float close = ans[0];
			if (ans[1] < ans[0] && ans[1] > 0)
				close = ans[1];
			if(close > 0.0001f && object[3] > 0 && isRef == -1 && isTrans == -1)
				return true;
		}
		return false;
}

float Game::solveForP(glm::vec4 plane, glm::vec3 vec, glm::vec3 pos)
{
	float noT = plane[0] * pos[0] + plane[1] * pos[1] + plane[2] * pos[2] + plane[3];
	float withT = plane[0] * vec[0] + plane[1] * vec[1] + plane[2] * vec[2];
	noT = -1 * noT;
	float t = noT / withT;
	if (t < 0.0001f)
		return -1.0f;
	if (t < 0.0f)
		t = 0.0;
	return t;
}

float Game::my_dist(glm::vec3 a, glm::vec3 b)
{
	float x = b[0] - a[0];
	float y = b[1] - a[1];
	float z = b[2] - a[2];
	float ans = x * x + y * y + z * z;
	ans = sqrt(ans);
	return ans;
}

glm::vec3 Game::calcColor(vec3 amb, int countLights, vec4* lights, vec4 obj, vec3 hitPos, vec3* colorLights,
	vec4* spotLight, vec4 objColor, vec3 tmpVec, int countObj, float dist, vec4* objects,
	vec3 pos, vec4* objectsColors, int isRef, int isTrans, int refPlace[], int transPlace[], glm::vec3 ambient, int countRef, vec3 norm)
{
	vec3 color = vec3(0.0, 0.0, 0.0);
	if (countRef > 5)
		return color;
	if (isRef >= 0) {
			vec3 refl = reflect(tmpVec, norm);
			refl = normalize(refl);
			glm::vec4 tmpObj = objects[0];
			glm::vec4 tmpObjColor = vec4(0.0, 0.0, 0.0, 0.0);
			glm::vec3 tmpHitPos = vec3(0.0, 0.0, 0.0);
			int isRef = -2;
			int isTrans = -2;
			float dist = -1;
			float ans[] = { -1.0, -1.0 };
			for (int i = 0; i < countObj; i++) {
				ans[0] = -1.0;
				ans[1] = -1.0;
				bool check = false;
				check = hitSphere(hitPos, objects[i], refl, ans);
				if (check) {
					float tmpHit = ans[0];
					if (ans[1] < ans[0] && ans[1] > 0.0)
						tmpHit = ans[1];
					glm::vec3 tmpHitPos2 = hitPos + refl * tmpHit;
					float tmpDist = my_dist(tmpHitPos, hitPos);
					if ((tmpDist < dist && tmpDist > 0) || (dist < 0 && tmpDist > 0)) {
						dist = tmpDist;
						tmpObj = objects[i];
						tmpObjColor = objectsColors[i];
						tmpHitPos = tmpHitPos2;
						isRef = refPlace[i];
						isTrans = transPlace[i];
					}


				}
			}
			if (dist > 0) {
				vec3 newNorm = glm::vec3(tmpObj);
				if (tmpObj[3] > 0) {
					newNorm = glm::vec3(tmpHitPos[0] - tmpObj[0], tmpHitPos[1] - tmpObj[1], tmpHitPos[2] - tmpObj[2]);
				}
				newNorm = glm::normalize(newNorm);
				vec3 newAmb = glm::vec3(tmpObjColor[0] * ambient[0], tmpObjColor[1] * ambient[1], tmpObjColor[2] * ambient[2]);
				color = amb + calcColor(newAmb, countLights, lights, tmpObj, tmpHitPos, colorLights, spotLight,
						tmpObjColor, refl, countObj, dist, objects, hitPos, objectsColors, isRef, isTrans, refPlace,
						transPlace, ambient, countRef + 1, newNorm);
			}
	}
	else if (isTrans >= 0) {
		float pi = 3.1415926535;
		float cosVal = dot(norm, tmpVec);
		float ang = std::acos(cosVal) * 180 / pi;
		float ang2 = std::asin((1.0f / 1.5f) * std::sin(ang)) * 180 / pi;
		vec3 inT = ((1.0f / 1.5f) * cosVal - std::cos(ang2)) * norm - ((1.0f / 1.5f) * tmpVec);
		inT = normalize(inT);

		float ans[] = { -1.0, -1.0 };
		hitSphere(hitPos, obj, inT, ans);
		float sec = ans[0];
		if (ans[1] > ans[0])
			sec = ans[1];
		if (sec > 0) {
			vec3 newHit = hitPos + inT * sec;
			vec3 newNorm = glm::vec3(obj);
			if (obj[3] > 0) {
				newNorm = glm::vec3(newHit[0] - obj[0], newHit[1] - obj[1], newHit[2] - obj[2]);
			}
			newNorm = glm::normalize(newNorm);

			cosVal = dot(newNorm, inT);
			ang = std::acos(cosVal) * 180 / pi;
			//cout << 1.5f * std::sin(ang) << endl;
			float checkAng = 1.5f * std::sin(ang);
			if (checkAng <= 1 && checkAng >= -1) {
				ang2 = std::asin(1.5f * std::sin(ang)) * 180 / pi;
				vec3 out = (1.5f * cosVal - std::cos(ang2)) * newNorm - 1.5f * inT;
				out = normalize(out);

				glm::vec4 tmpObj = objects[0];
				glm::vec4 tmpObjColor = vec4(0.0, 0.0, 0.0, 0.0);
				glm::vec3 tmpHitPos = vec3(0.0, 0.0, 0.0);
				int isRef = -2;
				int isTrans = -2;
				float dist = -1;
				ans[0] = -1.0;
				ans[1] = -1.0;
				for (int i = 0; i < countObj; i++) {
					ans[0] = -1.0;
					ans[1] = -1.0;
					bool check = false;
					check = hitSphere(newHit, objects[i], out, ans);
					if (check) {
						float tmpHit = ans[0];
						if (ans[1] < ans[0] && ans[1] > 0.0)
							tmpHit = ans[1];
						glm::vec3 tmpHitPos2 = newHit + out * tmpHit;
						float tmpDist = my_dist(tmpHitPos, newHit);
						if ((tmpDist < dist && tmpDist > 0) || (dist < 0 && tmpDist > 0)) {
							dist = tmpDist;
							tmpObj = objects[i];
							tmpObjColor = objectsColors[i];
							tmpHitPos = tmpHitPos2;
							isRef = refPlace[i];
							isTrans = transPlace[i];
						}


					}
				}
				if (dist > 0) {
					newNorm = glm::vec3(tmpObj);
					if (tmpObj[3] > 0) {
						newNorm = glm::vec3(tmpHitPos[0] - tmpObj[0], tmpHitPos[1] - tmpObj[1], tmpHitPos[2] - tmpObj[2]);
					}
					newNorm = glm::normalize(newNorm);
					vec3 newAmb = glm::vec3(tmpObjColor[0] * ambient[0], tmpObjColor[1] * ambient[1], tmpObjColor[2] * ambient[2]);
					color = calcColor(newAmb, countLights, lights, tmpObj, tmpHitPos, colorLights, spotLight,
						tmpObjColor, out, countObj, dist, objects, newHit, objectsColors, isRef, isTrans, refPlace,
						transPlace, ambient, countRef, newNorm) + color;
				}
			}
		}


	}
	else {
		color = amb;
		vec3 specular = vec3(0.7, 0.7, 0.7);
		glm::vec3 diffuse = glm::vec3(0.0, 0.0, 0.0);
		glm::vec3 specul = glm::vec3(0.0, 0.0, 0.0);
		for (int i = 0; i < countLights; i++) {

			glm::vec3 lightVec = glm::vec3(lights[i]);
			lightVec = glm::normalize(lightVec);
			glm::vec3 norm = glm::vec3(obj);
			if (obj[3] > 0) {
				norm = glm::vec3(hitPos[0] - obj[0], hitPos[1] - obj[1], hitPos[2] - obj[2]);
			}
			norm = glm::normalize(norm);
			glm::vec3 intensity = colorLights[i];

			if (lights[i][3] == 1.0f) {
				glm::vec3 position = glm::vec3(spotLight[i]);
				glm::vec3 dir = position - hitPos;
				dir = glm::normalize(dir);
				float spotAng = glm::max(glm::dot(-dir, lightVec), 0.0f);
				float dist2 = my_dist(hitPos, position);
				if (spotAng <= spotLight[i][3])
					intensity = glm::vec3(0.0f, 0.0f, 0.0f);
				else {
					float power = 5.0f * spotAng / dist2;
					intensity = intensity * power;
				}
			}

			//diffuse
			float d = glm::max(glm::dot(norm, lightVec), 0.0f);
			if (obj[3] > 0)
				d = glm::max(glm::dot(norm, -lightVec), 0.0f);
			glm::vec3 diff = glm::vec3(objColor[0], objColor[1], objColor[2]) * d * intensity;
			if ((int(1.5f * hitPos[0]) - 2.0f * glm::floor(int(1.5f * hitPos[0]) / 2.0f) ==
				int(1.5f * hitPos[1]) - 2.0f * glm::floor(int(1.5f * hitPos[1]) / 2.0f)) && obj[3] <= 0) {
				diff = 0.5f * diff;
			}


			//specular
			glm::vec3 ref1 = glm::reflect(lightVec, norm);
			ref1 = glm::normalize(ref1);
			float d3 = glm::max(glm::dot(ref1, -tmpVec), 0.0f);
			float tmp = d3;
			for (int i = 0; i < objColor[3]; i++) {
				d3 = d3 * tmp;
			}

			glm::vec3 spec = specular * d3 * intensity;

			color = color + diff + spec;
			diffuse = diffuse + diff;
			specul = specul + spec;
			for (int i = 0; i < countObj; i++)
				if (isHitBefore(-lightVec, objects[i], hitPos, isRef, isTrans) && obj[3] <= 0 && objects[i][3] > 0 && refPlace[i] == -1
					&& transPlace[i] == -1) {
					color = color*0.5f;
				}
		}
	}
	return color;
}


void Game::printVec(glm::vec3 a)
{
	cout << a[0];
	cout << " ";
	cout << a[1];
	cout << " ";
	cout << a[2] << endl;
}

glm::vec3 Game::calcRef()
{
	return glm::vec3();
}

void Game::printVec(glm::vec4 a)
{
	cout << a[0];
	cout << " ";
	cout << a[1];
	cout << " ";
	cout << a[2];
	cout << " ";
	cout << a[3] << endl;
}

Game::~Game(void)
{
}


