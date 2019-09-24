#include <iostream>
#include "GLFunctions.h"
#include "GLFW/glfw3.h"
#include "Context.h"
#include "RenderFunctions.h"
#include "util.h"
#include "Camera.h"
#include <tuple>
std::string TastyQuad::RenderFunctions::pathToShaderDeclarations = "./shaders/shaderDeclarations.glsl";
std::string TastyQuad::RenderFunctions::pathToShaderFolder = "./";

int const pixelCountX = 8;
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
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac, from imgui example
	glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GLFW_TRUE);
	window = glfwCreateWindow(256, 256, "TastyRay", NULL, NULL);
	
	glfwMakeContextCurrent(window);
	MYREQUIRE(window);
	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	}
	TastyQuad::util::setupGLDebugCallback();
}

struct RayMarchInstance {
	float minDist;
	float coveredDist;
	glm::vec2 uvA, uvB, uvC;
	TastyQuad::Material * resMat = nullptr;
	TastyQuad::TexturePowerOf2 * resAlbedo = nullptr;
	glm::vec3 position, dir;
	float incircleRadius;
	glm::vec3 incircleCtr;
};

int main(void) {
	initWindow();
	//load data
	TastyQuad::GPULessContext context(std::string("./Data/TinyForest/assets.db"));
	TastyQuad::GLTextureLoader glTexLoader;
	auto arena = context.load("monkeyLightCam");
	//create render texture
	auto newTex = context.getTextureLoader().createEmptyTexture(pixelCountX);
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
	int const maxMarchSteps = 4;
	float const sufficientDist = 0.1f;
	auto camT = context.findTransformByName("Camera");
	glm::mat4 V;
	glm::mat4 P;
	glm::mat4 W = camT->calculateWorldMatrix(context);
	glm::vec3 forward;

	glm::quat correction = glm::rotate(glm::quat(1, 0, 0, 0), -glm::pi<float>() / 2.0f, glm::vec3(1.0f, 0.0f, 0.0f));
	auto cameraCorrectionTransform = context.createTransformation(camT.get());
	cameraCorrectionTransform->modify([&](auto & pos, auto & rot, auto & scale) {
		pos = glm::vec3(0);
		rot = correction;
		scale = glm::vec3(1);
	});

	context.updateDirty();

	//initialize ray buffer
	std::vector< RayMarchInstance> rays;
	rays.resize(pixelCountX * pixelCountX);

	while (!glfwWindowShouldClose(window))
	{
		W = cameraCorrectionTransform->calculateWorldMatrix(context);
		forward = glm::vec3(W * glm::vec4(0, 0, -1, 1));
		TastyQuad::Camera::constructPerspViewProjection(fov, near, far, aspectRatio, glm::vec3(W[3]), glm::vec3(W[3]) + forward, V, P);
		

		std::cout << std::endl << "rays:" ;

		for (int i = 0; i < pixelCountX; i++) {
			std::cout << std::endl;
			for (int j = 0; j < pixelCountX; j++) {
				rays[i*pixelCountX + j].minDist = far;
				rays[i*pixelCountX + j].coveredDist = 0.0f;
				float viewPlaneWorldWidth = (glm::tan(fov / 2.0f) * near) * 2.0f;
				glm::vec2 normalizedScreenCoord = glm::vec2(static_cast<float>(j) / static_cast<float>(pixelCountX),
					static_cast<float>(i) / static_cast<float>(pixelCountX));
				rays[i*pixelCountX + j].position = glm::vec3(-0.5f * viewPlaneWorldWidth * aspectRatio + normalizedScreenCoord.x * viewPlaneWorldWidth * aspectRatio,
					-0.5f * viewPlaneWorldWidth * aspectRatio + normalizedScreenCoord.y * viewPlaneWorldWidth,
					near);
				rays[i*pixelCountX + j].dir = glm::normalize(rays[i*pixelCountX + j].position);
				rays[i*pixelCountX + j].position = glm::vec3(W * glm::vec4(rays[i*pixelCountX + j].position, 1));
				rays[i*pixelCountX + j].dir.z = -rays[i*pixelCountX + j].dir.z;
				rays[i*pixelCountX + j].dir = glm::mat3(W) * rays[i*pixelCountX + j].dir;
				std::cout << std::to_string(rays[i*pixelCountX + j].dir.z) << "  ";
			}
		}
		//march n times
		for (int iteration = 0; iteration < maxMarchSteps ; iteration++) {
			//visit the scene
			context.acceptVisitorGivePointers([&](auto * name, auto * scObj, TastyQuad::ITransformation * t, TastyQuad::Material * mat, auto * albedo, auto * normal, auto * mr, TastyQuad::Mesh * mesh) {
				if (mesh != nullptr) {
					//visit the mesh
					mesh->modify([&](std::pair<glm::vec3, glm::vec3> & boundingBox, std::vector<float> & vertices,
						std::vector<float> & normals, std::vector<float> & uvs, std::vector<unsigned int>& indices) {
						//for each ray
						for (int y = 0; y < pixelCountX; y++) {
							for (int x = 0; x < pixelCountX; x++) {
								if (rays[y*pixelCountX + x].coveredDist < far) { //only march rays that are still interesting
									//information per ray, local variables
									glm::vec3 position = rays[y*pixelCountX + x].position;
									glm::vec3 a, b, c, A, B, C, D, H, ctr; //triangle points, edges and new base vectors
									float r,area2, al,bl,cl,ls, dist; //transformed triangle point coordinates, squared lengths of base vectors
									dist = far;
									glm::mat4 worldT = t->calculateWorldMatrix(context);
									for (int i = 0; i < indices.size(); i = i + 3) { //
										a = glm::vec3(vertices[indices[i + 0] * 4 + 0], vertices[indices[i + 0] * 4 + 1], vertices[indices[i + 0] * 4 + 2]);
										b = glm::vec3(vertices[indices[i + 1] * 4 + 0], vertices[indices[i + 1] * 4 + 1], vertices[indices[i + 1] * 4 + 2]);
										c = glm::vec3(vertices[indices[i + 2] * 4 + 0], vertices[indices[i + 2] * 4 + 1], vertices[indices[i + 2] * 4 + 2]);
										a = glm::vec3(worldT * glm::vec4(a, 1));
										b = glm::vec3(worldT * glm::vec4(b, 1));
										c = glm::vec3(worldT * glm::vec4(c, 1));
										area2 = glm::cross(b - a, c - a).length();
										al = a.length();
										bl = b.length();
										cl = c.length();
										ls = al + bl + cl;
										r = area2 / ls;
										ctr = glm::mat3(a, b, c) * glm::vec3(al, bl, cl) / ls;
										//sphere dist
										dist = glm::distance(ctr,  position) - r;
										if (dist < rays[y*pixelCountX + x].minDist) { //found new result
											rays[y*pixelCountX + x].uvA = glm::vec2(uvs[indices[i + 0] * 4 + 0], uvs[indices[i + 0] * 4 + 1]);
											rays[y*pixelCountX + x].uvB = glm::vec2(uvs[indices[i + 1] * 4 + 0], uvs[indices[i + 1] * 4 + 1]);
											rays[y*pixelCountX + x].uvC = glm::vec2(uvs[indices[i + 2] * 4 + 0], uvs[indices[i + 2] * 4 + 1]);
											rays[y*pixelCountX + x].resMat = mat;
											rays[y*pixelCountX + x].resAlbedo = albedo;
											rays[y*pixelCountX + x].minDist = dist;
											rays[y*pixelCountX + x].incircleCtr = ctr;
											rays[y*pixelCountX + x].incircleRadius = r;
										}
									}
								}
							}
						}
					});
				}

			});
			std::cout << "march phase " << std::to_string(iteration) << std::endl;
			//for each ray
			for (int y = 0; y < pixelCountX ; y++) {
				for (int x = 0; x < pixelCountX; x++) {
					//march!
					rays[y*pixelCountX+x].position += rays[y*pixelCountX + x].dir * rays[y*pixelCountX + x].minDist;
					rays[y*pixelCountX + x].coveredDist += rays[y*pixelCountX + x].minDist;
					if(iteration < maxMarchSteps - 1)//if not last iteration, reset for next iteration
						rays[y*pixelCountX + x].minDist = far;
				}
			}
		}

		//paint phase
		newTex->modify([&](unsigned int & edgeLength, std::vector<unsigned char> & bytes) {
			//for each ray
			std::cout << std::endl << "dists:" ;
			for (int y = 0; y < pixelCountX ; y++) {
				std::cout << std::endl;
				for (int x = 0; x < pixelCountX; x++) {
					glm::vec4 color;
					float minDist = rays[y*pixelCountX + x].minDist;
					std::cout << std::to_string(minDist) << "  ";
					//found a close triangle ?
					if (minDist < sufficientDist) {
					
				/*pd = glm::clamp(rays[y*pixelCountX + x].pd / glm::sqrt(rays[y*pixelCountX + x].Ds), 0.0f, 1.0f);
						pc = glm::clamp(rays[y*pixelCountX + x].pc / glm::sqrt(rays[y*pixelCountX + x].Cs), 0.0f, 1.0f);
						pb = glm::clamp(rays[y*pixelCountX + x].pb / glm::sqrt(rays[y*pixelCountX + x].Bs), 0.0f, 1.0f);
						glm::vec2 uv = (rays[y*pixelCountX + x].uvB * pb + rays[y*pixelCountX + x].uvA * (1.0f - pb)) * (1.0f - pc) + rays[y*pixelCountX + x].uvC * pc;
						// uv - glm::vec2(glm::ivec2(uv));*/
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
					color = color * 255.0f;
					//write into color buffer
					bytes[edgeLength * y * 4 + 4 * x + 0] = static_cast<unsigned char>(color.r);// 0xff;
					bytes[edgeLength * y * 4 + 4 * x + 1] = static_cast<unsigned char>(color.g);
					bytes[edgeLength * y * 4 + 4 * x + 2] = static_cast<unsigned char>(color.b);
					bytes[edgeLength * y * 4 + 4 * x + 3] = static_cast<unsigned char>(color.a);

				}
			}
		});

		glTexLoader.reloadToGl(*newTex);


		//Switching between menu mode and game mode
		int state = glfwGetKey(window, GLFW_KEY_LEFT_ALT);
		if (state == GLFW_PRESS) {
			
		}

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
		glReadBuffer(GL_COLOR_ATTACHMENT0);
		glBlitFramebuffer(0, 0, pixelCountX, pixelCountX,
			0, 0, 256, 256, GL_COLOR_BUFFER_BIT, GL_NEAREST);
		glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);


		glfwSwapBuffers(window);

		glfwPollEvents();



	}

	glfwTerminate();
	return 0;

}