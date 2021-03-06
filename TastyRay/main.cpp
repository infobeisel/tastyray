#include <vector>
#include <iostream>
#include "GLFunctions.h"
#include "GLFW/glfw3.h"
#include "Context.h"
#include "RenderFunctions.h"
#include "util.h"
#include "Camera.h"
#include <tuple>
#include "SHBasis.h"
#define GLM_ENABLE_EXPERIMENTAL
#include "glm/gtx/transform.hpp"
#include <omp.h>
#include <glm/gtx/string_cast.hpp>
std::string TastyQuad::RenderFunctions::pathToShaderDeclarations = "./shaders/shaderDeclarations.glsl";
std::string TastyQuad::RenderFunctions::pathToShaderFolder = "./";

int pixelCountMin = 32;
int pixelCountMax = 512;

int pixelCountX = pixelCountMin;
int defaultFboWidth, defaultFboHeight;
GLFWwindow* window;
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	//defaultFboWidth = width;
	//defaultFboHeight = height;
	//glViewport(0, 0, width, height);
}

void initWindow() {
	int success = glfwInit();
	MYREQUIRE(success);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
 	glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GLFW_TRUE);
	window = glfwCreateWindow(1024, 1024, "TastyRay", NULL, NULL);
	
	glfwMakeContextCurrent(window);
	MYREQUIRE(window);
	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	}
	TastyQuad::util::setupGLDebugCallback();
}

struct RayInstance {
	TastyQuad::Material * resMat = nullptr;
	TastyQuad::TexturePowerOf2 * resAlbedo = nullptr;
	glm::vec3 position, dir, basisX, basisY;
	std::vector<glm::vec4> resColors;
	//sh "buffer"
	std::vector<float> shEvals;
	std::vector<float> shEvalsDX;
	std::vector<float> shEvalsDY;
	std::vector<float> shEvalsDZ;

	
};


const float paraboloidFar = 1000.0f;
const float paraboloidNear = 0.0f;

glm::vec2 paraboloidSample(glm::vec3 direction) {
	float dz = glm::length(direction);    
    glm::vec3 localParaboloidPos = glm::normalize(direction);    
    localParaboloidPos[2] = (dz - paraboloidNear) / (paraboloidFar - paraboloidNear);
    localParaboloidPos[0] /= (1.0 + localParaboloidPos[2]);
    localParaboloidPos[1] /= (1.0 + localParaboloidPos[2]);
    return glm::vec2(localParaboloidPos[0] , localParaboloidPos[1]);

}

int main(void) {
	initWindow();
	//load data
	TastyQuad::GPULessContext context(std::string("./Data/TinyForest/assets.db"));
	TastyQuad::GLTextureLoader glTexLoader;
	auto arena = context.load("shtestscene");

	

	//create render texture
	auto newTex = context.getTextureLoader().createEmptyTexture(pixelCountMax);
	glTexLoader.loadToGl(*newTex);
	GLuint fbo = 0;
	glGenFramebuffers(1, &fbo);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo);
	newTex->attachToFbo(fbo, GL_COLOR_ATTACHMENT0, 0);
	GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (!(status == GL_FRAMEBUFFER_COMPLETE)) {
		std::cout << "RenderFunctions::setupGBuffer(): incomplete gbuffer status " << status << std::endl;
	}

	float const fov = (50.0f / 180.0f) * glm::pi<float>();
	float const near = 0.5f;
	float const far = 1000.0f;
	float const aspectRatio = 1.0f;
	int const maxMarchSteps = 2;
	float const constStepSize = 1.0f / (float)maxMarchSteps;
	float const nextProbe = 0.0001f;
	auto camT = context.findTransformByName("Camera");
	int bounces = 1;
	float tolerance = 0.2f;

	

	glm::quat correction = glm::rotate(glm::quat(1, 0, 0, 0), -glm::pi<float>() / 2.0f, glm::vec3(1.0f, 0.0f, 0.0f));
	auto cameraCorrectionTransform = context.createTransformation(camT.get());
	cameraCorrectionTransform->modify([&](auto & pos, auto & rot, auto & scale) {
		pos = glm::vec3(0);
		rot = correction;
		scale = glm::vec3(1);
	});

	glm::mat4 V;
	glm::mat4 P;
	glm::mat4 W = camT->calculateWorldMatrix(context);
	glm::vec3 forward;
	
	context.updateDirty();

	//ugly cam
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	double oldCursorPosX = xpos;
	double oldCursorPosY = ypos;
	float cameraXAxisAngle = 0.0f;
	float cameraYAxisAngle = 0.0f;
	glm::mat4 newRx;
	glm::mat4 newRy;
	glm::mat4 id = glm::mat4(1);
	//----------

	//initialize ray buffer
	std::vector< RayInstance> rays;
	rays.resize(pixelCountX * pixelCountX);

	context.acceptVisitorGivePointers([&](auto * name, auto * scObj, TastyQuad::ITransformation * t, TastyQuad::Material * mat, auto * albedo, auto * normal, auto * mr, TastyQuad::Mesh * mesh) {
		glm::mat4 worldT = t->calculateWorldMatrix(context);
		glm::mat4 iWorldT = t->calculateWorldInverse(context);
		bool consider = false;
		auto extends = glm::vec3(0);
		auto shScale = 0.0f;
		auto shs = context.getSHs(*t);
		if (shs != nullptr) {//it's a sh object
			std::cout << "read a sh object at " << glm::to_string(iWorldT) << std::endl;
		}
	});

	double frameTimePoint = glfwGetTime();
	bool iPressed, kPressed;
	int showShIndex = 20;
	while (!glfwWindowShouldClose(window))
	{

		//ugly cam
		glfwGetCursorPos(window, &xpos, &ypos);
		
		//total camera rotation angles
		float newXAxisAngle = cameraXAxisAngle - 0.01f * (float)(ypos - oldCursorPosY);
		cameraYAxisAngle = cameraYAxisAngle - 0.01f * (float)(xpos - oldCursorPosX);
		//check for "salto" and > 360\B0
		//need to include the camera correction in checks
		cameraXAxisAngle = newXAxisAngle > 0.5f * glm::pi<float>() + (glm::pi<float>() / 2.0f) ? cameraXAxisAngle :
			(newXAxisAngle < -glm::pi<float>() * 0.5f + (glm::pi<float>() / 2.0f) ? cameraXAxisAngle : newXAxisAngle);
		cameraYAxisAngle = cameraYAxisAngle > glm::pi<float>() ? cameraYAxisAngle - (glm::pi<float>() * 2.0f) :
			(cameraYAxisAngle < -glm::pi<float>() ? cameraYAxisAngle + (glm::pi<float>() * 2.0f) : cameraYAxisAngle);
		//total camera rotation
		newRx = glm::rotate(id, cameraXAxisAngle, glm::vec3(1, 0, 0));
		newRy = glm::rotate(id, cameraYAxisAngle, glm::vec3(0, 1, 0));

		//calculate movement in local camera coordinates 

		auto cameraPosDelta = glm::vec4(0, 0, 0, 1);
		int state = glfwGetKey(window, GLFW_KEY_SPACE);
		int stateShift = glfwGetKey(window, GLFW_KEY_LEFT_SHIFT);
		float speed = state == GLFW_PRESS && stateShift == GLFW_PRESS ? 0.04f : (state == GLFW_PRESS ? 10.0f : 1.0f);
		state = glfwGetKey(window, GLFW_KEY_W);
		if (state == GLFW_PRESS) 
			cameraPosDelta += glm::vec4(0, 0, speed  * -1.0f, 0);
		state = glfwGetKey(window, GLFW_KEY_A);
		if (state == GLFW_PRESS)
			cameraPosDelta += glm::vec4(speed  * -1.0f, 0, 0, 0);
		state = glfwGetKey(window, GLFW_KEY_D);
		if (state == GLFW_PRESS)
			cameraPosDelta += glm::vec4(speed  * 1.0f, 0, 0, 0);
		state = glfwGetKey(window, GLFW_KEY_S);
		if (state == GLFW_PRESS)
			cameraPosDelta += glm::vec4(0, 0, speed  * 1.0f, 0);
		state = glfwGetKey(window, GLFW_KEY_I);
		if (state == GLFW_PRESS && !iPressed) {
			showShIndex = (showShIndex + 1) % 36;
			std::cout << std::to_string(showShIndex) << std::endl;
		}
		iPressed = state == GLFW_PRESS;

		state = glfwGetKey(window, GLFW_KEY_K);
		if (state == GLFW_PRESS && !kPressed) {
			showShIndex = (showShIndex - 1) % 36;
			std::cout << std::to_string(showShIndex) << std::endl;
		}
		kPressed = state == GLFW_PRESS;
		

		camT->modify([&](glm::vec3 & pos, glm::quat & rot) {
		//	rot = glm::quat_cast(newRy * newRx);
	//		pos += glm::vec3(newRy * newRx  * glm::mat4_cast(correction) * cameraPosDelta);
		});
		glm::vec2 movementDelta = glm::vec2((float)(xpos - oldCursorPosX), (float)(ypos - oldCursorPosY)) 
			+ glm::vec2(glm::length(glm::vec3(cameraPosDelta)));
		oldCursorPosX = xpos;
		oldCursorPosY = ypos;
		//----------


		float delta =   static_cast<float>(glfwGetTime() - frameTimePoint);
		if (delta > 1.0f && glm::length(movementDelta) < 0.01) {
			if(pixelCountX < pixelCountMax)
				pixelCountX *= 2;
			rays.resize(pixelCountX * pixelCountX);
			frameTimePoint = glfwGetTime();
		}
		else if (glm::length(movementDelta) > 0.01) {
			pixelCountX = pixelCountMin;
			frameTimePoint = glfwGetTime();
		}

		

		W = cameraCorrectionTransform->calculateWorldMatrix(context);
		glm::mat4 iW = cameraCorrectionTransform->calculateWorldInverse(context);
		forward = glm::vec3(W * glm::vec4(0, 0, -1, 1));
		TastyQuad::Camera::constructPerspViewProjection(fov, near, far, aspectRatio, glm::vec3(W[3]), glm::vec3(W[3]) + forward, V, P);
		
		#pragma omp parallel
		{
		#pragma omp for collapse(2)
		for (int i = 0; i < pixelCountX; i++) {
			for (int j = 0; j < pixelCountX; j++) {
				RayInstance & ray = rays[i*pixelCountX + j];
				

				float viewPlaneWorldWidth = (glm::tan(fov / 2.0f) * near) * 2.0f;
				glm::vec2 normalizedScreenCoord = glm::vec2(static_cast<float>(j) / static_cast<float>(pixelCountX),
					static_cast<float>(i) / static_cast<float>(pixelCountX));
				glm::vec2 normalizedScreenCoordNeighbor = glm::vec2(static_cast<float>(j+1) / static_cast<float>(pixelCountX),
					static_cast<float>(i+1) / static_cast<float>(pixelCountX));
				ray.position = glm::vec3(-0.5f * viewPlaneWorldWidth * aspectRatio + normalizedScreenCoord.x * viewPlaneWorldWidth * aspectRatio,
					-0.5f * viewPlaneWorldWidth * aspectRatio + normalizedScreenCoord.y * viewPlaneWorldWidth,
					near);
				ray.dir = glm::normalize(ray.position);
				ray.position = glm::vec3(W * glm::vec4(ray.position, 1));
				ray.dir.z = -ray.dir.z;
				ray.dir = glm::mat3(W) * ray.dir;
				auto neverColinear = ray.dir;
				neverColinear = glm::vec3(neverColinear.y,neverColinear.z,-neverColinear.x);
				ray.basisX = glm::normalize(glm::cross(neverColinear, ray.dir));
				ray.basisY =  glm::normalize(glm::cross( ray.dir, ray.basisX));
				ray.resColors.resize(0);
				ray.resColors.resize(bounces,glm::vec4(0.0f));
				
			}
		}

		#pragma omp barrier
		
		for(int bounce = 0; bounce < bounces; bounce++) {
			#pragma omp for collapse(2)
			for (int y = 0; y < pixelCountX; y++) {
				for (int x = 0; x < pixelCountX; x++) {
					RayInstance & ray = rays[y*pixelCountX + x];
					ray.shEvals.clear();
					ray.shEvals.resize(36,0.0f);
					auto dir = ray.dir;
					glm::mat3 rayBasis = glm::mat3 (
								ray.basisX,
								ray.basisY,
								ray.dir);
					float maxDist = far;
					context.acceptVisitorGivePointers([&](auto * name, auto * scObj, TastyQuad::ITransformation * transform, TastyQuad::Material * mat, auto * albedo, auto * normal, auto * mr, TastyQuad::Mesh * mesh) {
						glm::mat4 worldT = transform->calculateWorldMatrix(context);
						glm::mat4 iWorldT = transform->calculateWorldInverse(context);
						bool consider = false;
						auto extends = glm::vec3(0);
						auto shScale = 0.0f;
						auto shs = context.getSHs(*transform);
						if (shs != nullptr) {//it's a sh object
							shs->readRef([&](auto const & params, float const & radius) {
								transform->read([&](auto p, auto r, glm::vec3 scale) {
									shScale = scale.x;
									extends = glm::vec3(radius);
									extends = glm::abs(glm::mat3(worldT) * (extends  / shScale)); //divide by the sh local scale to get world extends
									consider = true;
								});
							});
						}
						glm::vec3 position = ray.position;
						auto s = glm::vec3(worldT[3]) - position;
						s = s * rayBasis;
						if (s.z - glm::sign(s.z) * glm::length(extends) > 0 && consider) { //sphere in front of ray
							float distToRay = glm::length(glm::vec2(s));
							if(distToRay < glm::length(extends)) { //ray intersects sphere. 
//-------------------------------- until now it was ray-boundingsphere-intersectiontest. from now on it is ray-sh-intersectiontest ----------------------------------------------------------------------
								float dist = glm::max(0.0f,s.z - glm::length(extends)); //start in front of sphere, but at least at 0.0 (already inside the sphere)
								glm::vec3 intersectionDir;
								glm::vec3 marchPos = position + dir * dist; 
								glm::vec3 testCol = glm::vec3(1);
								float toMinimize = std::numeric_limits<float>::max();
								glm::vec3 localRayDir = glm::mat3(iWorldT) * dir;
								localRayDir = glm::normalize(localRayDir);
								glm::vec3 localPos = glm::vec3(iWorldT * glm::vec4(marchPos,1));
								float t = 0.0f;
								//find good start position.
								auto kt = localPos;// + localRayDir * t;
								auto ktN = glm::normalize(kt);
								auto f = glm::length(kt);
								SHFunctions::SHEval6(ktN.x, ktN.y, ktN.z, &ray.shEvals[0]);
								float f_star = 0.0f;
								shs->readRef([&](auto & params, auto) {
									for(int i  = 0; i < 36; i++)
											f_star +=  ray.shEvals[i] * params[i];
									//f_star +=  ray.shEvals[showShIndex] ;//* params[i];
								});
								auto gt = glm::abs(f_star);
								auto kt_to_surface = ktN * gt - kt;
								float conservative = 1.0f / 1.0f;
								t = 0.77f * conservative * glm::dot(localRayDir, kt_to_surface); // project surface point on ray dir
								//std::cout << std::to_string(t) << std::endl;
								toMinimize = f - gt;

								for (int iteration = 0; glm::abs(toMinimize) > tolerance && iteration < maxMarchSteps && dist < maxDist; iteration++) {
									
									kt = localPos + localRayDir * t;
									ktN = glm::normalize(kt);
									f = glm::length(kt);
									auto df = glm::dot( kt / f , localRayDir);
									float te = t + nextProbe;
									auto kte = localPos + localRayDir * te;
									auto kteN = glm::normalize(kte);
									
									SHFunctions::SHEval6(ktN.x, ktN.y, ktN.z, &ray.shEvals[0]);
									f_star = 0.0f;
									shs->readRef([&](auto & params, auto) {
										for(int i  = 0; i < 36; i++)
											f_star +=  ray.shEvals[i] * params[i];
										//f_star +=  ray.shEvals[showShIndex] ;//* params[i];
									});

									SHFunctions::SHEval6(kteN.x, kteN.y, kteN.z, &ray.shEvals[0]);
									float f_starNext = 0.0f;
									shs->readRef([&](auto & params, auto) {
										for(int i  = 0; i < 36; i++)
											f_starNext +=  ray.shEvals[i] * params[i];
										//f_starNext +=  ray.shEvals[showShIndex] ;//* params[i];

									});

									gt = glm::abs(f_star);
									auto gte = glm::abs(f_starNext);
									auto dg = gte - gt;

									toMinimize = f - gt;

									//if(glm::abs(toMinimize) < tolerance) { //found intersection

									//} else {
										auto newtonDiv = df - dg;
										//if (glm::abs(newtonDiv) < std::numeric_limits<float>::epsilon()) {
											//no intersection, TODO opt out this case
											//t += constStepSize;
											
										//} else {
											t = t - conservative * (toMinimize / newtonDiv);
										//}
									//}
									intersectionDir = ktN;
								}
								if(glm::max(toMinimize,0.0f) < tolerance) {
									//std::cout << "hit  " << std::to_string(distToCenter) << std::endl;
									ray.resMat = mat;
									auto uv = paraboloidSample(intersectionDir);
										uv = uv * 0.5f + glm::vec2(0.5f);
									uv = glm::min(glm::vec2(1.0),glm::max(glm::vec2(0.0),uv));
									//ray.resAlbedo = albedo->sample(uv.x,uv.y);
									if(albedo && normal && mr) { //sh bake textures
										auto bakedNormal = intersectionDir.z > 0.0 ? albedo->sample(uv.x,uv.y) : mr->sample(uv.x,uv.y);
										bakedNormal = bakedNormal * 2.0f - glm::vec4(1.0f);
										auto bakedUv4 = normal->sample(uv.x,uv.y);
										auto bakedUv2 = intersectionDir.z > 0.0 ? glm::vec2(bakedUv4.x,bakedUv4.y) : glm::vec2(bakedUv4.z,bakedUv4.w);
										context.visitParent( [&] (auto * p) {
											auto * parentMat = context.getMaterial(*p);
											context.findCorrespondingTextures([&] (auto * a,auto * n,auto * mr ) {
												glm::vec4 sampledAlbedo = glm::vec4(1);

												if(a )
													sampledAlbedo = a->sample(bakedUv2.x,bakedUv2.y);
												if(parentMat) {
													parentMat->read([&] (auto aa,auto,auto,auto) {sampledAlbedo *= aa;});
												}
												//float l = glm::length(testCol);
												ray.resColors[bounce] = sampledAlbedo;//glm::vec4(bakedUv2,0,0);//sampledAlbedo;		

											}, *p);
										}, *transform);

										//prepare next trace
										
										ray.dir = glm::reflect(ray.dir,glm::normalize(glm::mat3(worldT) * glm::vec3(bakedNormal)));
										auto neverColinear = glm::vec3(ray.dir.y,ray.dir.z,-ray.dir.x);
										ray.basisX = glm::normalize(glm::cross(neverColinear, ray.dir));
										ray.basisY =  glm::normalize(glm::cross( ray.dir, ray.basisX));
										
									}
								}

								
								
							}
						} 
					});
				}
			}
		}
		#pragma omp barrier
		

		

		//paint phase
		newTex->modify([&] (unsigned int & edgeLength, std::vector<unsigned char> & bytes) {
			//for each ray
			#pragma omp for collapse(2)
			for (int y = 0; y < pixelCountX ; y++) {
				for (int x = 0; x < pixelCountX; x++) {
					glm::vec4 color = glm::vec4(0.0f);
					//if (rays[y*pixelCountX + x].hit) {
					for(auto c : rays[y*pixelCountX + x].resColors) color += c; 

					//}
					//else {
					//	color = glm::vec4(0);
					//}
					color = color * 255.0f;
					//write into color buffer
					bytes[edgeLength * y * 4 + 4 * x + 0] = static_cast<unsigned char>(color.r);// 0xff;
					bytes[edgeLength * y * 4 + 4 * x + 1] = static_cast<unsigned char>(color.g);
					bytes[edgeLength * y * 4 + 4 * x + 2] = static_cast<unsigned char>(color.b);
					bytes[edgeLength * y * 4 + 4 * x + 3] = static_cast<unsigned char>(color.a);
				}
			}
		});

		}

		glTexLoader.reloadToGl(*newTex);


		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
		glReadBuffer(GL_COLOR_ATTACHMENT0);
		glBlitFramebuffer(0, 0, pixelCountX, pixelCountX,
			0, 0, 1024, 1024, GL_COLOR_BUFFER_BIT, GL_NEAREST);
		glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);


		glfwSwapBuffers(window);

		glfwPollEvents();



	}

	glfwTerminate();
	return 0;

}