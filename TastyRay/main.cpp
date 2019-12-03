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
	bool sh = false;
	TastyQuad::Material * resMat = nullptr;
	TastyQuad::TexturePowerOf2 * resAlbedo = nullptr;
	glm::vec3 position, dir, basisX, basisY;
	bool hit = false;
	//sh "buffer"
	std::vector<float> shEvals;
	
};

int main(void) {
	initWindow();
	//load data
	TastyQuad::GPULessContext context(std::string("./Data/TinyForest/assets.db"));
	TastyQuad::GLTextureLoader glTexLoader;
	auto arena = context.load("monkey");
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
	int const maxMarchSteps = 30;
	float const sufficientDist = 0.001f;
	float const constStepSize = 0.033f;
	auto camT = context.findTransformByName("Camera");
	

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

	

	double frameTimePoint = glfwGetTime();
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
		state = glfwGetKey(window, GLFW_KEY_F);

		camT->modify([&](glm::vec3 & pos, glm::quat & rot) {
			rot = glm::quat_cast(newRy * newRx);
			pos += glm::vec3(newRy * newRx  * glm::mat4_cast(correction) * cameraPosDelta);
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
		
		//std::cout << "cam pos " <<
		#pragma omp parallel
		{

		//std::cout << "num threads" <<  omp_get_num_threads() << std::endl;
		//std::cout << std::endl << "rays:" ;
		#pragma omp for collapse(2)
		for (int i = 0; i < pixelCountX; i++) {
			//std::cout << std::endl;
			for (int j = 0; j < pixelCountX; j++) {
				

				rays[i*pixelCountX + j].sh = false;
				rays[i*pixelCountX + j].hit = false;
				rays[i*pixelCountX + j].shEvals.resize(36,0.0f);
				
				float viewPlaneWorldWidth = (glm::tan(fov / 2.0f) * near) * 2.0f;
				glm::vec2 normalizedScreenCoord = glm::vec2(static_cast<float>(j) / static_cast<float>(pixelCountX),
					static_cast<float>(i) / static_cast<float>(pixelCountX));
				glm::vec2 normalizedScreenCoordNeighbor = glm::vec2(static_cast<float>(j+1) / static_cast<float>(pixelCountX),
					static_cast<float>(i+1) / static_cast<float>(pixelCountX));
				rays[i*pixelCountX + j].position = glm::vec3(-0.5f * viewPlaneWorldWidth * aspectRatio + normalizedScreenCoord.x * viewPlaneWorldWidth * aspectRatio,
					-0.5f * viewPlaneWorldWidth * aspectRatio + normalizedScreenCoord.y * viewPlaneWorldWidth,
					near);
				
				rays[i*pixelCountX + j].dir = glm::normalize(rays[i*pixelCountX + j].position);
				rays[i*pixelCountX + j].position = glm::vec3(W * glm::vec4(rays[i*pixelCountX + j].position, 1));
				rays[i*pixelCountX + j].dir.z = -rays[i*pixelCountX + j].dir.z;
				rays[i*pixelCountX + j].dir = glm::mat3(W) * rays[i*pixelCountX + j].dir;
				auto neverColinear = rays[i*pixelCountX + j].dir;
				neverColinear = glm::vec3(neverColinear.y,neverColinear.z,-neverColinear.x);
				rays[i*pixelCountX + j].basisX = glm::normalize(glm::cross(neverColinear, rays[i*pixelCountX + j].dir));
				rays[i*pixelCountX + j].basisY =  glm::normalize(glm::cross( rays[i*pixelCountX + j].dir, rays[i*pixelCountX + j].basisX));
				//std::cout << std::to_string(rays[i*pixelCountX + j].dir.z) << "  ";
			}
		}

		#pragma omp barrier
		
		#pragma omp for collapse(2)
		//for each ray
		for (int y = 0; y < pixelCountX; y++) {
			for (int x = 0; x < pixelCountX; x++) {
				
				//information per ray, local variables
				auto dir = rays[y*pixelCountX + x].dir;
				glm::mat3 rayBasis = glm::mat3 (
							rays[y*pixelCountX + x].basisX,
							rays[y*pixelCountX + x].basisY,
							rays[y*pixelCountX + x].dir);
				float maxDist = far;
				rays[y*pixelCountX + x].hit = false;
								
				//visit the scene
				context.acceptVisitorGivePointers([&](auto * name, auto * scObj, TastyQuad::ITransformation * t, TastyQuad::Material * mat, auto * albedo, auto * normal, auto * mr, TastyQuad::Mesh * mesh) {
					glm::mat4 worldT = t->calculateWorldMatrix(context);
					bool consider = false;
					auto extends = glm::vec3(0);
					auto shScale = 0.0f;
					auto shs = context.getSHs(*t);
					if (shs != nullptr) {//it's a sh object
						shs->readRef([&](auto const & params, float const & radius) {
							t->read([&](auto p, auto r, glm::vec3 scale) {
								shScale = scale.x;
								extends = glm::vec3(radius);
								extends = glm::mat3(worldT) * (extends  / shScale); //divide by the sh local scale to get world extends
								consider = true;
							});
						});
					}
					//project sphere into ray basis
					glm::vec3 position = rays[y*pixelCountX + x].position;
					auto s = glm::vec3(worldT[3]) - position;
					s = s * rayBasis;
					if (s.z - glm::sign(s.z) * extends.x > 0 && consider) { //sphere in front of ray
						
						float distToRay = glm::length(glm::vec2(s));
						if(distToRay < extends.x) { //ray intersects sphere

							float dist = glm::max(0.0f,s.z - extends.x); //start in front of sphere, but at least at 0.0 (already inside the sphere)
							bool hit = false;
							//maxDist = dist;
							rays[y*pixelCountX + x].resMat = mat;
							rays[y*pixelCountX + x].resAlbedo = albedo;
							rays[y*pixelCountX + x].hit = true;
							//march n times
							//hit = true;

							for (int iteration = 0; iteration < maxMarchSteps && dist < maxDist && !hit; iteration++) {
								glm::mat4 iWorldT = t->calculateWorldInverse(context);
								glm::vec3 marchPos = position + dir * dist; 
								glm::vec3 unitLocalDir = glm::normalize(glm::vec3(iWorldT * glm::vec4(marchPos, 1)));
								SHFunctions::SHEval6(unitLocalDir.x, unitLocalDir.y, unitLocalDir.z, &rays[y*pixelCountX + x].shEvals[0]);
								float f_star = 0.0f;
								shs->readRef([&](auto & params, auto) {
									for (int i = 0; i < 36; i++) {
										f_star += rays[y*pixelCountX + x].shEvals[i] * params[i];
									}
								});
								auto distToCenter = glm::length(marchPos - glm::vec3(worldT[3]));
								//std::cout << std::to_string(glm::abs(f_star) * shScale) << " VS " << std::to_string(distToCenter) << std::endl;
								if (glm::abs(f_star) * shScale  > distToCenter) {
									rays[y*pixelCountX + x].sh = true;
									rays[y*pixelCountX + x].resMat = mat;
									rays[y*pixelCountX + x].resAlbedo = albedo;
									rays[y*pixelCountX + x].hit = true;
									hit = true;
								}
								dist += constStepSize * extends.x;
							}
							if(hit) {
								maxDist = dist;
							}
							
						}
					} 
				});
			}
		}

		#pragma omp barrier
		

		

		//paint phase
		newTex->modify([&] (unsigned int & edgeLength, std::vector<unsigned char> & bytes) {
			//for each ray
			#pragma omp for collapse(2)

			//std::cout << std::endl << "paint! dists:" ;
			for (int y = 0; y < pixelCountX ; y++) {
				//std::cout << std::endl;
				for (int x = 0; x < pixelCountX; x++) {
					glm::vec4 color;
					//std::cout << std::to_string(minDist) << "  ";
					//found a close triangle ?
					if (rays[y*pixelCountX + x].hit) {
			
						glm::vec2 uv = glm::vec2(0);
						uv.x = uv.x < 0.0f ? uv.x + 1.0f : uv.x;
						uv.y = uv.y < 0.0f ? uv.y + 1.0f : uv.y;
						if (rays[y*pixelCountX + x].resMat != nullptr) {
							rays[y*pixelCountX + x].resMat->read([&](auto alb, auto e, auto m, auto r) {
								color = (rays[y*pixelCountX + x].resAlbedo != nullptr ? rays[y*pixelCountX + x].resAlbedo->sample(uv.x, uv.y) : glm::vec4(1.0f)) * alb;
							});
						}
						else
							color = glm::vec4(1);
						color = glm::vec4(1);
					}
					else {
						color = glm::vec4(0);
					}
					color *= rays[y*pixelCountX + x].sh ? glm::vec4(1, 0, 0, 1) : glm::vec4(1);
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