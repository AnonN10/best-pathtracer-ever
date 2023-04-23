#include <SDL.h>

#include <iostream>

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/ext/matrix_projection.hpp>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define PI 3.141592653589793238462643

//sky
glm::vec3 sky_color = glm::vec3(0.1, 0.4, 1.0);
glm::vec3 sun_color = glm::vec3(1.0, 1.0, 1.0);
float sun_intensity = 1.0f;
glm::vec3 sun_dir = glm::normalize(glm::vec3(0.0, -1.0, 0.0));
float sun_angle = 3.1415926f*0.2f;
glm::vec3 calc_sky_color(glm::vec3 dir)
{
	float dir_up_cosine = glm::dot(dir, glm::vec3(0, 1, 0));
	float dir_sun_angle = glm::acos(glm::dot(dir, -sun_dir));
	return (dir_sun_angle < sun_angle)?(sun_color*sun_intensity):(sky_color*(glm::clamp(1.0f-dir_up_cosine*0.7f, 0.3f, 1.0f)));
	//return vec3(1.0);//white furnace
	//return vec3(0.5);//gray furnace
	//return vec3(0.0);
}

//intersection code
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

//low discrepancy sequence
glm::uvec2 Sobol(glm::uint n) {
    glm::uvec2 p = glm::uvec2(0u);
    glm::uvec2 d = glm::uvec2(0x80000000u);

    for(; n != 0u; n >>= 1u) {
        if((n & 1u) != 0u)
            p ^= d;
        
        d.x >>= 1u; // 1st dimension Sobol matrix, is same as base 2 Van der Corput
        d.y ^= d.y >> 1u; // 2nd dimension Sobol matrix
    }
    
    return p;
}

glm::uint ReverseBits(glm::uint x) {
    x = ((x & 0xaaaaaaaau) >> 1) | ((x & 0x55555555u) << 1);
    x = ((x & 0xccccccccu) >> 2) | ((x & 0x33333333u) << 2);
    x = ((x & 0xf0f0f0f0u) >> 4) | ((x & 0x0f0f0f0fu) << 4);
    x = ((x & 0xff00ff00u) >> 8) | ((x & 0x00ff00ffu) << 8);
    return (x >> 16) | (x << 16);
}

glm::uint OwenHash(glm::uint x, glm::uint seed) { // works best with random seeds
    x ^= x * 0x3d20adeau;
    x += seed;
    x *= (seed >> 16) | 1u;
    x ^= x * 0x05526c56u;
    x ^= x * 0x53a22864u;
    return x;
}

glm::uint OwenScramble(glm::uint p, glm::uint seed) {
    p = ReverseBits(p);
    p = OwenHash(p, seed);
    return ReverseBits(p);
}

//sampling
glm::vec3 map_to_unit_sphere(glm::vec2 uv)
{
    float cos_theta = 2.0*uv.x-1.0;
    float phi = 2.0*PI*uv.y;
    float sin_theta = glm::sqrt(1.0-cos_theta*cos_theta);
    float sin_phi = glm::sin(phi);
    float cos_phi = glm::cos(phi);
    
    return glm::vec3(sin_theta*cos_phi, cos_theta, sin_theta*sin_phi);
}

glm::vec3 map_to_unit_hemisphere_cosine_weighted(glm::vec2 uv, glm::vec3 n)
{
    glm::vec3 p = map_to_unit_sphere(uv);
    return n+p;
}

//etc
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

	SDL_Window * window = SDL_CreateWindow("Behold", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, win_width, win_height, 0);

	SDL_Renderer * renderer = SDL_CreateRenderer(window, -1, 0);
	SDL_Texture * texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STATIC, win_width, win_height);

	SDL_PixelFormat* pxfmt = SDL_AllocFormat(SDL_PIXELFORMAT_ARGB8888);	

	glm::vec4* color_buffer = new glm::vec4[win_width * win_height];

	Uint32 * pixels = new Uint32[win_width * win_height];
	memset(pixels, 255, win_width * win_height * sizeof(Uint32));
	
	//load blue noise
	int blue_noise_width, blue_noise_height, blue_noise_components;
	const uint8_t* blue_noise_data = stbi_load("bluenoise256.png", &blue_noise_width, &blue_noise_height, &blue_noise_components, 4);
	std::cout << blue_noise_components << std::endl;
	glm::vec4* blue_noise_pixels = new glm::vec4[blue_noise_width * blue_noise_height];
	for(int i = 0; i < blue_noise_width * blue_noise_height; ++i)
		blue_noise_pixels[i] = glm::vec4(
			blue_noise_data[i*4+0],
			blue_noise_data[i*4+1],
			blue_noise_data[i*4+2],
			blue_noise_data[i*4+3]
		) / 255.0f;
		
	//rendering variables
	float sample_count = 0.0f;
 
    while (!quit)
    {
		SDL_Event event;
		SDL_UpdateTexture(texture, NULL, pixels, win_width * sizeof(Uint32));
        SDL_PollEvent(&event);
		
		static glm::vec3 translation = glm::vec3(0, 0, 3);
		
		bool reset_rendering = false;
        switch (event.type)
        {
            case SDL_QUIT: quit = true; break;
			case SDL_KEYDOWN:
                /* Check the SDLKey values and move change the coords */
                switch( event.key.keysym.sym ){
                    case SDLK_LEFT:
                        translation.x -= 1;
						reset_rendering = true;
                        break;
                    case SDLK_RIGHT:
                        translation.x += 1;
						reset_rendering = true;
                        break;
                    case SDLK_UP:
                        translation.z -= 1;
						reset_rendering = true;
                        break;
                    case SDLK_DOWN:
                        translation.z += 1;
						reset_rendering = true;
                        break;
                    default:
                        break;
                }
        }
		
		if(reset_rendering)
		{
			for(int i = 0; i < win_width * win_height; ++i) color_buffer[i] = glm::vec4(0.0, 0.0, 0.0, 1.0);
			sample_count = 0.0f;
		}
		
		//rendering code
		glm::mat4 view = glm::translate(translation);
		glm::mat4 proj = glm::perspective(glm::radians(90.0f), float(win_width)/float(win_height), 0.1f, 100.0f);
		glm::mat4 viewproj = proj * glm::inverse(view);
		glm::mat4 invviewproj = glm::inverse(viewproj);
		static glm::uint frame_idx = 0;
		for(int y = 0; y < win_height; ++y)
		{
			for(int x = 0; x < win_width; ++x)
			{
				//read previous color
				glm::vec3 color = glm::vec3(color_buffer[x+win_width*y]);
				
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
					//sample sky by bouncing off the sphere
					
					glm::uvec2 ip = Sobol(frame_idx);
					ip.x = OwenScramble(ip.x, 0xe7843fbau);
					ip.y = OwenScramble(ip.y, 0x8d8fb1e9u);
					glm::vec2 rand_uv = glm::vec2(ip) / float(0xffffffffu);
					rand_uv = glm::fract(rand_uv + glm::vec2(blue_noise_pixels[x%blue_noise_width + y%blue_noise_height*blue_noise_height]));
					glm::vec3 incident_dir = glm::normalize(map_to_unit_hemisphere_cosine_weighted(rand_uv, normal));
					//float cosine = glm::dot(normal, incident_dir);
															
					//pretty much all terms are canceled out by cosine weighted sampling with pdf=cosine/PI,
					//would have been this otherwise: cosine * radiance * albedo/PI / pdf
					color += calc_sky_color(incident_dir);
				}
				else
				{
					color += calc_sky_color(ray_dir);
				}
				
				//store new color
				color_buffer[x+win_width*y] = glm::vec4(color, 1.0);
			}
		}
		sample_count += 1.0f;
		++frame_idx;
		
		//write float color buffer to the texture pixels
		for(int y = 0; y < win_height; ++y)
		{
			for(int x = 0; x < win_width; ++x)
			{
				int i = x+win_width*y;
				//NOTE the sample count division, doing it here to not involve another intermediate buffer
				auto color = color_buffer[i] / sample_count;
				pixels[i] = SDL_MapRGBA(
					pxfmt,
					glm::clamp(color.r, 0.0f, 1.0f)*255,
					glm::clamp(color.g, 0.0f, 1.0f)*255,
					glm::clamp(color.b, 0.0f, 1.0f)*255,
					255
				);
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