#include <SDL.h>

#include <iostream>

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/ext/matrix_projection.hpp>

glm::vec2 raySphere(glm::vec3 rayOrigin, glm::vec3 rayDir, glm::vec3 sphereCentre, float sphereRadius) {
    glm::vec3 offset = rayOrigin - sphereCentre;
    float a = 1.0; // Set to dot(rayDir, rayDir) if rayDir might be not normalized
    float b = 2.0 * dot(offset, rayDir);
    float c = dot(offset, offset) - sphereRadius * sphereRadius;
    float d = b * b - 4.0 * a * c; // Discriminant from quadratic formula
    // Number of intersections: 0 when d < 0; 1 when d = 0; 2 when d > 0
    if(d > 0.0){
        float s = sqrt(d);
        float dstToSphereNear = glm::max(0.0, (-b - s) / (2.0 * a));
        float dstToSphereFar = (-b + s) / (2.0 * a);
        // Ignore intersections that occur behind the ray
        if(dstToSphereFar >= 0.0) {
            return glm::vec2(dstToSphereNear, dstToSphereFar);
        }
    }
    // Ray did not intersect sphere
    return glm::vec2(-1.0, 0);
}

glm::vec3 rayTriangle(glm::vec3 ray_orig, glm::vec3 ray_dir, glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, float tmin, float tmax)
{
	glm::vec3 e1 = v0 - v1;
	glm::vec3 e2 = v2 - v0;
	glm::vec3 n = cross(e1, e2);
	
	glm::vec3 c = v0 - ray_orig;
	glm::vec3 r = glm::cross(ray_dir, c);
	float inv_det = 1.0 / glm::dot(n, ray_dir);

	float u = glm::dot(r, e2) * inv_det;
	float v = glm::dot(r, e1) * inv_det;
	float w = 1.0 - u - v;

	// These comparisons are designed to return false
	// when one of t, u, or v is a NaN
	if (u >= 0 && v >= 0 && w >= 0) {
		float t = glm::dot(n, c) * inv_det;
		if (t >= tmin && t <= tmax)
		//if (t >= 0.0)
			return glm::vec3(t, u, v);
	}

	return glm::vec3(FLT_MAX, 0.0, 0.0);
}

glm::vec3 UnprojectNDC(glm::vec3 ndc, glm::mat4 invmvp)
{
	glm::vec4 world = invmvp * glm::vec4(ndc, 1.0);
	return glm::vec3(world) / world.w;
}
 
int main(int argc, char ** argv)
{
    bool quit = false;
	int win_width = 640, win_height = 480;
 
    SDL_Init(SDL_INIT_VIDEO);	
 
    SDL_Window * window = SDL_CreateWindow("Thing", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, win_width, win_height, 0);
	
	SDL_Renderer * renderer = SDL_CreateRenderer(window, -1, 0);
	SDL_Texture * texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STATIC, win_width, win_height);
	
	SDL_PixelFormat* pxfmt = SDL_AllocFormat(SDL_PIXELFORMAT_ARGB8888);	
	
	glm::vec4* color_buffer = new glm::vec4[win_width * win_height];
	
	Uint32 * pixels = new Uint32[win_width * win_height];
	memset(pixels, 255, win_width * win_height * sizeof(Uint32));
	
	bool leftMouseButtonDown = false;
 
    while (!quit)
    {
		SDL_Event event;
		SDL_UpdateTexture(texture, NULL, pixels, win_width * sizeof(Uint32));
        SDL_PollEvent(&event);
		
		static glm::vec3 translation = glm::vec3(0, 0, 3);
 
        switch (event.type)
        {
            case SDL_QUIT: quit = true; break;
			case SDL_KEYDOWN:
                /* Check the SDLKey values and move change the coords */
                switch( event.key.keysym.sym ){
                    case SDLK_LEFT:
                        translation.x -= 1;
                        break;
                    case SDLK_RIGHT:
                        translation.x += 1;
                        break;
                    case SDLK_UP:
                        translation.z -= 1;
                        break;
                    case SDLK_DOWN:
                        translation.z += 1;
                        break;
                    default:
                        break;
                }
        }
		
		//rendering code
		glm::mat4 view = glm::translate(translation);
		glm::mat4 proj = glm::perspective(glm::radians(90.0f), float(win_width)/float(win_height), 0.1f, 100.0f);
		glm::mat4 viewproj = proj * glm::inverse(view);
		glm::mat4 invviewproj = glm::inverse(viewproj);
		for(int y = 0; y < win_height; ++y)
		{
			for(int x = 0; x < win_width; ++x)
			{
				glm::vec4 color = glm::vec4(0, 0, 0, 1);
				
				glm::vec2 ndc_xy = glm::vec2(float(x)/float(win_width)*2.0-1.0, -(float(y)/float(win_height)*2.0-1.0));				
				glm::vec3 ray_orig = UnprojectNDC(glm::vec3(ndc_xy, 0.0f), invviewproj);
				glm::vec3 ray_dir = UnprojectNDC(glm::vec3(ndc_xy, 1.0f), invviewproj);
				ray_dir = glm::normalize(ray_dir - ray_orig);
				
				//glm::vec3 rt_result = rayTriangle(ray_orig, ray_dir, glm::vec3(1, 0, 0), glm::vec3(-1, 0, 0), glm::vec3(0, 1, 0), 0.0, FLT_MAX);
				//if(rt_result.x != FLT_MAX)
				glm::vec3 sphere_orig = glm::vec3(0.0, 0.0, 0.0);
				glm::vec2 rt_result = raySphere(ray_orig, ray_dir, sphere_orig, 1.0f);
				if(rt_result.x >= 0.0)
				{
					glm::vec3 normal = glm::normalize(ray_orig + ray_dir*rt_result.x - sphere_orig);
					color = glm::vec4((glm::vec3(normal.x, normal.y, normal.z)+1.0f)*0.5f, 1.0);
				}
				color_buffer[x+win_width*y] = color;
			}
		}
		
		for(int y = 0; y < win_height; ++y)
		{
			for(int x = 0; x < win_width; ++x)
			{
				int i = x+win_width*y;
				Uint32 pixel = SDL_MapRGBA(
					pxfmt,
					glm::clamp(color_buffer[i].r, 0.0f, 1.0f)*255,
					glm::clamp(color_buffer[i].g, 0.0f, 1.0f)*255,
					glm::clamp(color_buffer[i].b, 0.0f, 1.0f)*255,
					255
				);
				pixels[i] = pixel;
			}
		}
		
		SDL_RenderClear(renderer);
		SDL_RenderCopy(renderer, texture, NULL, NULL);
		SDL_RenderPresent(renderer);
    }
 
    SDL_DestroyWindow(window);
    SDL_Quit();
	
	delete[] pixels;
SDL_DestroyTexture(texture);
SDL_DestroyRenderer(renderer);
 
    return 0;
}