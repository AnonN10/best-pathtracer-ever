#include <SDL.h>

#include <iostream>

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/ext/matrix_projection.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define PI 3.141592653589793238462643

#define ENABLE_NEE

//sky
glm::vec3 sky_color = glm::vec3(0.2, 0.55, 1.0);
glm::vec3 sun_color = glm::vec3(1.0, 0.6, 0.3);
float sun_intensity = 150.0f;
glm::vec3 sun_dir = glm::normalize(glm::vec3(2.0, -1.0, 1.5));
float sun_angle = 3.1415926f*0.03f;
glm::vec3 calc_sky_color(glm::vec3 dir)
{
	float dir_up_cosine = glm::dot(dir, glm::vec3(0, 1, 0));
	float dir_sun_angle = glm::acos(glm::dot(dir, -sun_dir));
	return (dir_sun_angle < sun_angle)?(sun_color*sun_intensity):(sky_color*(glm::clamp(1.0f-dir_up_cosine*0.7f, 0.3f, 1.0f)));
	//return vec3(1.0);//white furnace
	//return vec3(0.5);//gray furnace
	//return vec3(0.0);
}

//transformations
struct Transform
{
	glm::vec3 translation = glm::vec3(0);
	glm::vec3 scaling = glm::vec3(1);
	glm::quat rotation = glm::quat();
	
	glm::mat4 Compose() {return glm::translate(translation) * glm::toMat4(rotation) * glm::scale(scaling);}
};

//generates local coordinates given a normal vector.
//taken from "Building an Orthonormal Basis, Revisited", Duff et al.
glm::mat3 BuildLocalCoords(const glm::vec3& n)
{
    float s = (n.z < 0.0)?-1.0:1.0;
    float a = -1.0 / (s + n.z);
    float b = n.x * n.y * a;
    glm::vec3 t = glm::vec3(1.0 + s * n.x* n.x * a, s * b, -s * n.x);
    glm::vec3 bt = glm::vec3(b, s + n.y * n.y * a, -n.y);
    return glm::mat3(t, n, bt);
}

//materials
struct Material
{
	glm::vec3 albedo = glm::vec3(1.0);
	glm::vec3 emission = glm::vec3(0.0);
};

//intersection code
struct IntersectionData
{
	float t, u, v, w;
	glm::vec3 normal;
	Material* material;
	
	void reset(){t = -1.0f;}
	bool hit(){return t >= 0.0f;}
};

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

//scene
constexpr int num_spheres = 3;
glm::vec3 sphere_positions[num_spheres] = {glm::vec3(-2, 0, -5), glm::vec3(1, 1, -5), glm::vec3(0, -1001, -5)};
float sphere_radii[num_spheres] = {1.0f, 2.0f, 1000.0f};
Material sphere_materials[num_spheres] = {
	{glm::vec3(1.0, 0.1, 0.1), glm::vec3(0.0)},
	{glm::vec3(0.1, 1.0, 0.1), glm::vec3(0.0)},
	{glm::vec3(1.0, 1.0, 1.0), glm::vec3(0.0)}
};

IntersectionData TraceScene(const glm::vec3& ray_orig, const glm::vec3& ray_dir)
{
	IntersectionData isect_data;
	isect_data.t = -1.0f;
	for(int i = 0; i < num_spheres; ++i)
	{
		glm::vec2 rs_result = raySphere(ray_orig, ray_dir, sphere_positions[i], sphere_radii[i]);
		if(rs_result.x >= 0.0 && (!isect_data.hit() || rs_result.x < isect_data.t))
		{
			isect_data.t = rs_result.x;
			//no barycentrics on the sphere, so u, v, w are ignored
			isect_data.normal = glm::normalize(ray_orig + ray_dir*rs_result.x - sphere_positions[i]);
			isect_data.material = &sphere_materials[i];
		}
	}
	
	return isect_data;
}

//hash
glm::uint WangHash(glm::uint seed)
{
    seed = (seed ^ 61u) ^ (seed >> 16u);
    seed *= 9u;
    seed = seed ^ (seed >> 4u);
    seed *= 0x27d4eb2du;
    seed = seed ^ (seed >> 15u);
    return 1u+seed;
}

//low discrepancy sequence
//n-dimensional R-sequence by Martin Roberts,
//from "The Unreasonable Effectiveness of Quasirandom Sequences"
float RLDS(float phi_d, float d, float n) {
	float a = 1.0f/glm::pow(phi_d, d);
	return glm::fract(5.0f + a*n);
}

float RLDSComputePhi(glm::uint d_max) {
	float phi_d = 2.0f;
	for(int i = 0; i < 15; ++i)
	{
		phi_d = glm::pow(1.0f + phi_d, 1.0f/float(phi_d+1));
	}
	return phi_d;
}

//owen scrambling
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

float OwenScrambleFloat(float p, glm::uint seed) {
    return float(OwenScramble(glm::uint(p*float(0xffffffffu)), seed)) / float(0xffffffffu);
}

//sampling
glm::vec3 MapToUnitSphere(glm::vec2 uv)
{
    float cos_theta = 2.0*uv.x-1.0;
    float phi = 2.0*PI*uv.y;
    float sin_theta = glm::sqrt(1.0-cos_theta*cos_theta);
    float sin_phi = glm::sin(phi);
    float cos_phi = glm::cos(phi);
    
    return glm::vec3(sin_theta*cos_phi, cos_theta, sin_theta*sin_phi);
}

glm::vec3 MapToUnitHemisphereCosineWeighted(glm::vec2 uv, glm::vec3 n)
{
    glm::vec3 p = MapToUnitSphere(uv);
    return n+p;
}

float MapToUnitHemisphereCosineWeightedPDF(float cos_theta)
{
	return cos_theta / PI;
}

glm::vec3 MapToSolidAngle(const glm::vec2& uv, float theta_max)
{
	float phi = 2.0 * PI * uv.x;
	float cos_theta = 1.0 - uv.y * (1.0 - glm::cos(theta_max));
	float sin_theta = glm::sqrt(1.0 - cos_theta * cos_theta);
	return glm::vec3(glm::cos(phi) * sin_theta, cos_theta, glm::sin(phi) * sin_theta);
}

float MapToSolidAnglePDF(float theta_max)
{
	return 1.0f/(2.0f * PI * (1.0f - glm::cos(theta_max)));
}

glm::vec3 SampleSkySun(const glm::vec2& uv, const glm::vec3& normal, const glm::mat3& sun_transform, float& cosine, float& PDF)
{	
	glm::vec3 dir = sun_transform * MapToSolidAngle(uv, sun_angle);
	cosine = glm::clamp(glm::dot(normal, dir), 0.0f, 1.0f);
	PDF = MapToSolidAnglePDF(sun_angle);
    return dir;
}

float SampleSkySunPDF(glm::vec3 dir)
{
	float dir_sun_angle = glm::acos(glm::clamp(glm::dot(dir, -sun_dir), 0.0f, 1.0f));
	return (dir_sun_angle < sun_angle)?MapToSolidAnglePDF(sun_angle):0.0f;
}

//MIS
float BalanceHeuristic(float p_f, float p_g)
{
	return p_f/(p_f+p_g);
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
	glm::vec4* blue_noise_pixels = new glm::vec4[blue_noise_width * blue_noise_height];
	for(int i = 0; i < blue_noise_width * blue_noise_height; ++i)
		blue_noise_pixels[i] = glm::vec4(
			blue_noise_data[i*4+0],
			blue_noise_data[i*4+1],
			blue_noise_data[i*4+2],
			blue_noise_data[i*4+3]
		) / 255.0f;
		
	//rendering variables
	Transform camera_transform;
	glm::vec2 rotation_angles = glm::vec2(0, 0);
	constexpr int max_bounces = 10;
	float sample_count = 0.0f;
	glm::mat3 sun_transform = BuildLocalCoords(-sun_dir);
	
	//R LDS
	constexpr glm::uint RLDS_d_max = max_bounces*2;
	float RLDS_phi_d = RLDSComputePhi(RLDS_d_max);
	//std::cout << "d_max = " << RLDS_d_max << std::endl;
	//std::cout << "phi_d = " << RLDS_phi_d << std::endl;
	
	//Owen scrambling
	glm::uint lds_scramble_seeds[max_bounces*2];
	for(glm::uint i = 0; i < max_bounces*2; ++i)
		lds_scramble_seeds[i] = WangHash(i);
 
    while (!quit)
    {
		SDL_Event event;
		SDL_UpdateTexture(texture, NULL, pixels, win_width * sizeof(Uint32));
        SDL_PollEvent(&event);		
		
		bool reset_rendering = false;
        switch (event.type)
        {
            case SDL_QUIT: quit = true; break;
        }
		
		const Uint8* state = SDL_GetKeyboardState(nullptr);
		if (state[SDL_SCANCODE_KP_PLUS]) {
			camera_transform.scaling *= 2.0;
			reset_rendering = true;
		}
		if (state[SDL_SCANCODE_KP_MINUS]) {
			camera_transform.scaling /= 2.0;
			reset_rendering = true;
		}
		if (state[SDL_SCANCODE_LEFT]) {
			rotation_angles.x += 0.1;
			reset_rendering = true;
		}
		if (state[SDL_SCANCODE_RIGHT]) {
			rotation_angles.x -= 0.1;
			reset_rendering = true;
		}
		if (state[SDL_SCANCODE_DOWN]) {
			rotation_angles.y -= 0.1;
			reset_rendering = true;
		}
		if (state[SDL_SCANCODE_UP]) {
			rotation_angles.y += 0.1;
			reset_rendering = true;
		}
		camera_transform.rotation = glm::angleAxis(rotation_angles.x, glm::vec3(0, 1, 0)) * glm::angleAxis(rotation_angles.y, glm::vec3(1, 0, 0));		
		if (state[SDL_SCANCODE_A]) {
			camera_transform.translation += camera_transform.rotation * (glm::vec3(-1, 0, 0) * camera_transform.scaling);
			reset_rendering = true;
		}
		if (state[SDL_SCANCODE_D]) {
			camera_transform.translation += camera_transform.rotation * (glm::vec3(1, 0, 0) * camera_transform.scaling);
			reset_rendering = true;
		}
		if (state[SDL_SCANCODE_Q]) {
			camera_transform.translation += camera_transform.rotation * (glm::vec3(0, -1, 0) * camera_transform.scaling);
			reset_rendering = true;
		}
		if (state[SDL_SCANCODE_E]) {
			camera_transform.translation += camera_transform.rotation * (glm::vec3(0, 1, 0) * camera_transform.scaling);
			reset_rendering = true;
		}
		if (state[SDL_SCANCODE_W]) {
			camera_transform.translation += camera_transform.rotation * (glm::vec3(0, 0, -1) * camera_transform.scaling);
			reset_rendering = true;
		}
		if (state[SDL_SCANCODE_S]) {
			camera_transform.translation += camera_transform.rotation * (glm::vec3(0, 0, 1) * camera_transform.scaling);
			reset_rendering = true;
		}
		//just reset rendering by holding
		if (state[SDL_SCANCODE_R]) {reset_rendering = true;}
		//got tired of rendering
		if (state[SDL_SCANCODE_ESCAPE]) {quit = true; break;}
		
		if(reset_rendering)
		{
			for(int i = 0; i < win_width * win_height; ++i) color_buffer[i] = glm::vec4(0.0, 0.0, 0.0, 1.0);
			sample_count = 0.0f;
		}
		
		//rendering code
		static glm::uint frame_idx = 0;
		
		glm::mat4 view = camera_transform.Compose();
		glm::mat4 proj = glm::perspective(glm::radians(90.0f), float(win_width)/float(win_height), 0.1f, 100.0f);
		glm::mat4 viewproj = proj * glm::inverse(view);
		glm::mat4 invviewproj = glm::inverse(viewproj);
		
		constexpr int lds_seq_length = 1;
		float lds_points[RLDS_d_max][lds_seq_length];
		for(int d = 0; d < RLDS_d_max; ++d)
			for(int i = 0; i < lds_seq_length; ++i)
			{
				lds_points[d][i] = RLDS(RLDS_phi_d, d+1, frame_idx * lds_seq_length + i);
				//lds_points[d][i] = OwenScrambleFloat(lds_points[d][i], lds_scramble_seeds[d]);
			}
		
		/*sample_count = 1.0f;
		glm::ivec2 lds_coord = glm::ivec2(glm::vec2(lds_points[2][0], lds_points[3][0]) * glm::vec2(win_width, win_height));
		if(lds_coord.x < win_width && lds_coord.x >= 0 && lds_coord.y < win_height && lds_coord.y >= 0)
		{
			color_buffer[lds_coord.x+win_width*lds_coord.y] += glm::vec4(glm::vec3(0.5f), 1.0f);
		}*/
		
		for(int y = 0; y < win_height; ++y)
		{
			for(int x = 0; x < win_width; ++x)
			{
				//read previous color
				glm::vec3 color = glm::vec3(color_buffer[x+win_width*y]);
				
				//portion of radiance that passes from the last vertex of the path to the sensor, a R[0; 1]^3 vector value
				glm::vec3 throughput = glm::vec3(1);
				
				//obtain NDC space screen coordinate
				glm::vec2 ndc_xy = glm::vec2(float(x)/float(win_width)*2.0-1.0, -(float(y)/float(win_height)*2.0-1.0));
				//add jitter for anti-aliasing
				ndc_xy += glm::vec2(lds_points[0][0], lds_points[1][0])/glm::vec2(win_width, win_height)*2.0f;
				
				//primary ray
				glm::vec3 ray_orig = UnprojectNDC(glm::vec3(ndc_xy, 0.0f), invviewproj);
				glm::vec3 ray_dir = UnprojectNDC(glm::vec3(ndc_xy, 1.0f), invviewproj);
				ray_dir = glm::normalize(ray_dir - ray_orig);
				
				
				for(int bounce = 0; bounce < max_bounces; ++bounce)
				{
					IntersectionData isect_data = TraceScene(ray_orig, ray_dir);
					float cosine;
					if(isect_data.hit())
					{
						//add contribution from the surface emission
						color += throughput * isect_data.material->emission;
						
						//fetch blue noise mask value for this bounce
						glm::vec2 blue_noise = glm::vec2(blue_noise_pixels[(x+1324*bounce)%blue_noise_width + (y+764352*bounce)%blue_noise_height*blue_noise_height]);
						//do the blue noise dithering of the low discrepancy sequence
						glm::vec2 rand_uv = glm::fract(glm::vec2(lds_points[bounce*2+0][0], lds_points[bounce*2+1][0]) + blue_noise);
						
						//make new ray for a higher dimensional integral of the BRDF estimator
						ray_orig = ray_orig + ray_dir * isect_data.t;
						ray_orig += isect_data.normal * 0.0005f;//ray bias
						ray_dir = glm::normalize(MapToUnitHemisphereCosineWeighted(rand_uv, isect_data.normal));
						cosine = glm::clamp(glm::dot(isect_data.normal, ray_dir), 0.0f, 1.0f);
						
#ifdef ENABLE_NEE
						//sample the sun disk in the sky directly with a Next Event Estimator
						float nee_cosine;
						float nee_pdf;
						glm::vec3 nee_ray_dir = SampleSkySun(
							/*glm::vec2(
								OwenScrambleFloat(rand_uv.x, lds_scramble_seeds[bounce*2+0]),
								OwenScrambleFloat(rand_uv.y, lds_scramble_seeds[bounce*2+1])
							)*/rand_uv,
							isect_data.normal,
							sun_transform,
							nee_cosine,
							nee_pdf
						);
						bool is_light_visible = !TraceScene(ray_orig, nee_ray_dir).hit();
						//adding bounce!=max_bounces-1 condition makes the image match with the image without NEE,
						//otherwise NEE is always 1 bounce ahead of the renderer without it, making it impossible to validate
						if(is_light_visible && nee_cosine > 0 && bounce!=max_bounces-1)
						{
							//compute Multiple Importance Sampling weight and add the weighted contribution
							float nee_mis_weight = BalanceHeuristic(nee_pdf, MapToUnitHemisphereCosineWeightedPDF(nee_cosine));
							color += nee_mis_weight * throughput * (isect_data.material->albedo/float(PI)) * nee_cosine * calc_sky_color(nee_ray_dir) / nee_pdf;
						}
#endif
						
						//pretty much all the terms are canceled out by the cosine weighted hemisphere sampling with BRDF=albedo/PI and PDF=cosine/PI,
						//would have been this otherwise: cosine * BRDF / PDF
						throughput *= isect_data.material->albedo;
					}
					else
					{
						//the ray didn't hit any geometry, so fetch the sky color and multiply by the portion of how much light goes to the sensor
#ifdef ENABLE_NEE
						//must not ignore sun for primary rays, hence bounce!=0 check is added
						float mis_weight = bounce!=0? BalanceHeuristic(MapToUnitHemisphereCosineWeightedPDF(cosine), SampleSkySunPDF(ray_dir)): 1.0f;						
						color += mis_weight * throughput * calc_sky_color(ray_dir);
#else
						color += throughput * calc_sky_color(ray_dir);						
#endif
						break;
					}
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
				//approximate sRGB transfer
				color = glm::pow(color, glm::vec4(1.0f/2.2f));
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

	delete[] color_buffer;
	delete[] pixels;
	SDL_DestroyTexture(texture);
	SDL_DestroyRenderer(renderer);
 
    return 0;
}