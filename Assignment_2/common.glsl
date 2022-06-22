/**
 * common.glsl
 * Common types and functions used for ray tracing.
 */

const float pi = 3.14159265358979;
const float epsilon = 0.001;

struct Ray {
    vec3 o;     // origin
    vec3 d;     // direction - always set with normalized vector
    float t;    // time, for motion blur
};

Ray createRay(vec3 o, vec3 d, float t)
{
    Ray r;
    r.o = o;
    r.d = d;
    r.t = t;
    return r;
}

Ray createRay(vec3 o, vec3 d)
{
    return createRay(o, d, 0.0);
}

vec3 pointOnRay(Ray r, float t)
{
    return r.o + r.d * t;
}

float gSeed = 0.0;

uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    uvec2 rz = uvec2(n, n * 48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU)) / float(0x7fffffff);
}

vec3 hash3(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
    uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
    return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

float rand(vec2 v)
{
    return fract(sin(dot(v.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 toLinear(vec3 c)
{
    return pow(c, vec3(2.2));
}

vec3 toGamma(vec3 c)
{
    return pow(c, vec3(1.0 / 2.2));
}

vec2 randomInUnitDisk(inout float seed) {
    vec2 h = hash2(seed) * vec2(1.0, 6.28318530718);
    float phi = h.y;
    float r = sqrt(h.x);
	return r * vec2(sin(phi), cos(phi));
}

vec3 randomInUnitSphere(inout float seed)
{
    vec3 h = hash3(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0);
    float phi = h.y;
    float r = pow(h.z, 1.0/3.0);
	return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}

vec3 randomUnitVector(inout float seed) //to be used in diffuse reflections with distribution cosine
{
    return(normalize(randomInUnitSphere(seed)));
}

struct Camera
{
    vec3 eye;
    vec3 u, v, n;
    float width, height;
    float lensRadius;
    float planeDist, focusDist;
    float time0, time1;
};

Camera createCamera(
    vec3 eye,
    vec3 at,
    vec3 worldUp,
    float fovy,
    float aspect,
    float aperture,  //diametro em multiplos do pixel size
    float focusDist,  //focal ratio
    float time0,
    float time1)
{
    Camera cam;
    if(aperture == 0.0) cam.focusDist = 1.0; //pinhole camera then focus in on vis plane
    else cam.focusDist = focusDist;
    vec3 w = eye - at;
    cam.planeDist = length(w);
    cam.height = 2.0 * cam.planeDist * tan(fovy * pi / 180.0 * 0.5);
    cam.width = aspect * cam.height;

    cam.lensRadius = aperture * 0.5 * cam.width / iResolution.x;  //aperture ratio * pixel size; (1 pixel=lente raio 0.5)
    cam.eye = eye;
    cam.n = normalize(w);
    cam.u = normalize(cross(worldUp, cam.n));
    cam.v = cross(cam.n, cam.u);
    cam.time0 = time0;
    cam.time1 = time1;
    return cam;
}

Ray getRay(Camera cam, vec2 pixel_sample)  //rnd pixel_sample viewport coordinates
{
    vec2 ls = cam.lensRadius * randomInUnitDisk(gSeed);  //ls - lens sample for DOF
    float time = cam.time0 + hash1(gSeed) * (cam.time1 - cam.time0);

    // Powerpoint Distribution Ray-Tracing, Slide 39
    vec3 eye_offset;
    vec3 ray_dir; 

    vec3 ps;

    ps.x = cam.width * (pixel_sample.x / iResolution.x - 0.5);
    ps.y = cam.height * (pixel_sample.y / iResolution.y - 0.5);

    vec3 p;

    p.x = ps.x * cam.focusDist;
    p.y = ps.y * cam.focusDist;
    
    //d = normalize((px - lsx)u + (py - lsy)v - fz
    ray_dir = (cam.u * (p.x - ls.x) + cam.v * (p.y - ls.y) + cam.n * cam.focusDist * -cam.planeDist);
   
    //eye_offset = eye + lsx * u + lsy * v
    eye_offset = cam.eye + (cam.u * ls.x) + (cam.v * ls.y);


    return createRay(eye_offset, normalize(ray_dir), time);
}

// MT_ material type
#define MT_DIFFUSE 0
#define MT_METAL 1
#define MT_DIALECTRIC 2

struct Material
{
    int type;
    vec3 albedo;  //diffuse color
    vec3 specColor;  //the color tint for specular reflections. for metals and opaque dieletrics like coloured glossy plastic
    vec3 emissive; //
    float roughness; // controls roughness for metals. It can be used for rough refractions
    float refIdx; // index of refraction for dialectric
    vec3 refractColor; // absorption for beer's law
    float refractionRoughness; //roughness for refraction
};

Material createDiffuseMaterial(vec3 albedo)
{
    Material m;
    m.type = MT_DIFFUSE;
    m.albedo = albedo;
    m.specColor = vec3(0.0);
    m.roughness = 1.0;  //ser usado na iluminação direta
    m.refIdx = 1.0;
    m.refractColor = vec3(0.0);
    m.emissive = vec3(0.0);
    return m;
}

Material createMetalMaterial(vec3 specClr, float roughness)
{
    Material m;
    m.type = MT_METAL;
    m.albedo = vec3(0.0);
    m.specColor = specClr;
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

Material createDialectricMaterial(vec3 refractClr, float refIdx, float roughness, float refractionRoughness)
{
    Material m;
    m.type = MT_DIALECTRIC;
    m.albedo = vec3(0.0);
    m.specColor = vec3(0.04);
    m.refIdx = refIdx;
    m.refractColor = refractClr;  
    m.roughness = roughness;
    m.refractionRoughness = refractionRoughness;
    m.emissive = vec3(0.0);
    return m;
}

struct HitRecord
{
    vec3 pos;
    vec3 normal;
    float t;            // ray parameter
    Material material;
};


float schlick(float cosine, float ior_1, float ior_t)
{
    // Powerpoint -> Whitted Ray - tracing: Practice, slide 50 
    float R0_aux = (ior_1 - ior_t) / (ior_1 + ior_t);
    float R0 = pow(R0_aux, 2.0);
    
    float Kr = R0 + (1.0 - R0) * pow(1.0 - cosine, 5.0);
    return Kr;
}

bool scatter(Ray rIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{
    if(rec.material.type == MT_DIFFUSE)
    {
        // consider a unit radius sphere tangent to the hit point and
        // calculate point S on its surface
        vec3 S = rec.pos + rec.normal + randomInUnitSphere(gSeed);
        
        // pointing from the ray position to S
        vec3 rayDirection = normalize(S - rec.pos);
        
        // creates ray from the hit point with the new direction
        rScattered = createRay(rec.pos + rec.normal * epsilon, rayDirection, rIn.t);
        
        atten = rec.material.albedo * max(dot(rScattered.d, rec.normal), 0.0) / pi;
        return true;
    }
    if(rec.material.type == MT_METAL)
    {
        // calulate direction taking into account roughness
        // ray.direction = normalize(R + roughness * rand_in_unit_sphere())
        vec3 rayDirection = rIn.d - 2.0 * rec.normal * dot(rIn.d, rec.normal);
        rayDirection = normalize(rayDirection + randomInUnitSphere(gSeed) * rec.material.roughness);

        // creates ray from the hit point with the new direction
        rScattered = createRay(rec.pos + rec.normal * epsilon, rayDirection, rIn.t);
        
        atten = rec.material.specColor;
        return true;
    }
    if(rec.material.type == MT_DIALECTRIC)
    {
        atten = vec3(1.0);  
        vec3 outwardNormal;
        float niOverNt;
        float cosine;
        float ior_1;
        float ior_t;

        if(dot(rIn.d, rec.normal) >= 0.0) //hit inside
        {
            outwardNormal = -rec.normal;
            niOverNt = rec.material.refIdx;
            cosine = dot(rIn.d, rec.normal);
           
            ior_1 = rec.material.refIdx;
            ior_t = 1.0;

            atten = exp(-rec.material.refractColor * rec.t);  //beer's law
        }
        else  //hit from outside
        {
            outwardNormal = rec.normal;
            niOverNt = 1.0 / rec.material.refIdx;
            cosine = -dot(rIn.d, rec.normal); 

            ior_t = 1.0;
            ior_1 = rec.material.refIdx; 
        }

        //Use probabilistic math to decide if scatter a reflected ray or a refracted ray

        float reflectProb;

        //https://blog.demofox.org/2020/06/14/casual-shadertoy-path-tracing-3-fresnel-rough-refraction-absorption-orbit-camera/
        float discriminant = 1.0 - niOverNt * niOverNt * (1.0 - cosine * cosine);

        if (discriminant > 0.0){ //if no total reflection  reflectProb = schlick(cosine, rec.material.refIdx);  
            reflectProb = schlick(cosine, ior_1, ior_t);
        } else {
            reflectProb = 1.0;
        }

        if( hash1(gSeed) < reflectProb) { //Reflection

            // reflection vector = 2(L*n)n - L
            vec3 rayDirection = normalize(rIn.d - 2.0 * outwardNormal * dot(rIn.d, outwardNormal)); 
            // takes into account roughness to create the new direction
            rayDirection = normalize(rayDirection + randomInUnitSphere(gSeed) * rec.material.roughness);
            
            // creates ray from the hit point with the new direction
            rScattered = createRay(rec.pos + rec.normal * epsilon, rayDirection, rIn.t);

            //atten *= vec3(reflectProb); not necessary since we are only scattering reflectProb rays and not all reflected rays
        
        }else{  //Refraction 
            // as described on page 30 of "Raytracing in One Weekend" by Peter Shirley
            vec3 rayDirection = normalize(niOverNt * rIn.d + (niOverNt * cosine - sqrt(discriminant)) * outwardNormal);

            //https://blog.demofox.org/2020/06/14/casual-shadertoy-path-tracing-3-fresnel-rough-refraction-absorption-orbit-camera/
            rayDirection = normalize(mix(rayDirection, normalize(outwardNormal + randomInUnitSphere(gSeed)), rec.material.refractionRoughness * rec.material.refractionRoughness));
             
            // creates ray from the hit point with the new direction
            rScattered = createRay(rec.pos - outwardNormal * epsilon, rayDirection, rIn.t);
            //atten *= vec3(1.0 - reflectProb); not necessary since we are only scattering 1-reflectProb rays and not all refracted rays
        }

        return true;
    }
    return false;
}

struct Triangle {vec3 a; vec3 b; vec3 c; };

Triangle createTriangle(vec3 v0, vec3 v1, vec3 v2)
{
    Triangle t;
    t.a = v0; t.b = v1; t.c = v2;
    return t;
}

bool hit_triangle(Triangle t, Ray r, float tmin, float tmax, out HitRecord rec)
{
    // https://www.youtube.com/watch?v=fK1RPmF_zjQ 
    vec3 vert0 = t.a; 
    vec3 vert1 = t.b;
    vec3 vert2 = t.c;

    // get edges of the triangle (sharing v0)
	vec3 edge1 = vert1 - vert0;
	vec3 edge2 = vert2 - vert0;

    vec3 cross_rayDir_edge2 = cross(r.d, edge2);
    float det = dot(edge1, cross_rayDir_edge2);

    if (det > -0.0000001f && det < 0.0000001f) return false; //ray is parallel with triangle (n podemos dividir por 0)

    float inv_det = 1.0f / det;

    // distance from ray origin to v0 (t)   
    vec3 orig_minus_vert0 = r.o - vert0;

    //calculate baryU
	float u = dot(inv_det, dot(orig_minus_vert0, cross_rayDir_edge2));

    // check if the values are within the triangle (bary not bigger than 1 and not smaller than 0 if it is within)
	if (u < 0.0f || u > 1.0f) return false;

   vec3 cross_origMinusVert0_edge1 = cross(orig_minus_vert0, edge1);

    //calculate baryV
   float v =  dot(inv_det, dot(r.d, cross_origMinusVert0_edge1));

    // check if the values are within the triangle (bary not bigger than 1 and not smaller than 0 if it is within) 
	if (v < 0.0f || u + v > 1.0f) return false;

    float taux = dot(inv_det, dot(edge2, cross_origMinusVert0_edge1));
    vec3 normal = normalize(cross(edge1, edge2));
    
    //calculate a valid t and normal
    if(taux < tmax && taux > tmin)
    {
        rec.t = taux;
        rec.normal = normal;
        rec.pos = pointOnRay(r, rec.t);
        return true;
    }
    return false;
}



struct Sphere
{
    vec3 center;
    float radius;
};

Sphere createSphere(vec3 center, float radius)
{
    Sphere s;
    s.center = center;
    s.radius = radius;
    return s;
}


struct MovingSphere
{
    vec3 center0, center1;
    float radius;
    float time0, time1;
};

MovingSphere createMovingSphere(vec3 center0, vec3 center1, float radius, float time0, float time1)
{
    MovingSphere s;
    s.center0 = center0;
    s.center1 = center1;
    s.radius = radius;
    s.time0 = time0;
    s.time1 = time1;
    return s;
}

vec3 center(MovingSphere mvsphere, float time)
{
    return mvsphere.center0 + (mvsphere.center1 - mvsphere.center0) * ((time - mvsphere.time0) / (mvsphere.time1 - mvsphere.time0));
}


/*
 * The function naming convention changes with these functions to show that they implement a sort of interface for
 * the book's notion of "hittable". E.g. hit_<type>.
 */

bool hit_sphere(Sphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    vec3 OC = s.center - r.o; // center of sphere - origin of ray
    
    float b = dot(OC, r.d); // d * OC
    float c = dot(OC, OC) - s.radius * s.radius; //OC * OC - r^2
    float t = 0.0;
	
    if(c > 0.0) { // ray origin is outside
        if (b <= 0.0){ // sphere is behind of the ray and return false
            return false;
        }
    }

    float discriminant = b * b - c; // discriminant = b ^ 2 - c
    if (discriminant <= 0.0){ // if discriminant is less or equal to zero, there is no hit
        return false;
    }

    if(c > 0.0){ // if origin of ray is outside, compute the smallest root
        t = b - sqrt(discriminant);
    } else { // if origin of ray is inside, compute the positive root
        t = b + sqrt(discriminant);
    }

    if(t < tmax && t > tmin) {
        rec.t = t;
        rec.pos = pointOnRay(r, rec.t);

        if(s.radius >= 0.0){
            rec.normal = normalize(rec.pos - s.center);
        } else{
            rec.normal = normalize(s.center - rec.pos);
        }
        return true;
    }
    else return false;
}

bool hit_movingSphere(MovingSphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    vec3 OC = center(s, r.t) - r.o; // center of sphere - origin of ray
    
    float b = dot(OC, r.d); // d * OC
    float c = dot(OC, OC) - s.radius * s.radius; //OC * OC - r^2
    float t = 0.0;
	
    if(c > 0.0) { // ray origin is outside
        if (b <= 0.0){ // sphere is behind of the ray and return false
            return false;
        }
    }

    float discriminant = b * b - c; // discriminant = b ^ 2 - c
    if (discriminant <= 0.0){ // if discriminant is less or equal to zero, there is no hit
        return false;
    }

    if(c > 0.0){ // if origin of ray is outside, compute the smallest root
        t = b - sqrt(discriminant);
    } else { // if origin of ray is inside, compute the positive root
        t = b + sqrt(discriminant);
    }

    if(t < tmax && t > tmin) {
        rec.t = t;
        rec.pos = pointOnRay(r, rec.t);

        if(s.radius >= 0.0){
            rec.normal = normalize(rec.pos - center(s, r.t));
        } else{
            rec.normal = normalize(center(s, r.t) - rec.pos);
        }
        return true;
    }
    else return false;
}

struct pointLight {
    vec3 pos;
    vec3 color;
};

pointLight createPointLight(vec3 pos, vec3 color) 
{
    pointLight l;
    l.pos = pos;
    l.color = color;
    return l;
}