//
// Implementation for Yocto/Model
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_model.h"

#include <yocto/yocto_sampling.h>

#include "ext/perlin-noise/noise1234.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::vector;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EXAMPLE OF PROCEDURAL MODELING
// -----------------------------------------------------------------------------
namespace yocto {

float noise(const vec3f& p) { return ::noise3(p.x, p.y, p.z); }
vec2f noise2(const vec3f& p) {
  return {noise(p + vec3f{0, 0, 0}), noise(p + vec3f{3, 7, 11})};
}
vec3f noise3(const vec3f& p) {
  return {noise(p + vec3f{0, 0, 0}), noise(p + vec3f{3, 7, 11}),
      noise(p + vec3f{13, 17, 19})};
}
float fbm(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 1.0f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * fabs(noise(p * scale));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}
vec3f hash3(vec3f x) {
  auto  p = vec3f{x.x, x.y,x.z};
  auto  q = vec3f{dot(p, vec3f{127.1,269.1 ,311.7}), dot(p, vec3f{11.6, 237.1,311.7}),
      dot(p, vec3f{163.1,269.5 ,183.3})};
  vec3f z = vec3f{
      sin(q.x) * 43758.5453f, sin(q.y) * 43758.5453f, sin(q.z) * 43758.5453f};
  vec3f r = {z.x - floor(z.x), z.y - floor(z.y), z.z - floor(z.z)};
  return r;
}
vec3f hash2(vec2f x) {
  auto  p = vec2f{x.x, x.y};
  auto  q = vec2f{dot(p, vec2f{127.1,311.7}), dot(p, vec2f{269.5,183.3})};
  vec2f z = vec2f{sin(p.x) * 43758.5453f, sin(q.y) * 43758.5453f};
  vec3f r = {z.x - floor(z.x), z.y - floor(z.y)};
  return r;
}
float turbulence(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 1.0f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * fabs(noise(p * scale));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}
float ridge(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 0.5f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * (1 - fabs(noise(p * scale))) * (1 - fabs(noise(p * scale)));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}
bool Isvalid3d(vec3f pos, vec3f S, float s, float r, int*** g, vector<vec3f> p) {
  if (pos.x >= 0.f && pos.x <= S.x && pos.y >= 0.f && pos.y <= S.y &&
      pos.z >= 0.f && pos.z <= S.z) {
    int cell_index_x = (int)(pos.x / s);
    int cell_index_y = (int)(pos.y / s);
    int cell_index_z = (int)(pos.z / s);

    int start_x = max(0, cell_index_x - 2);
    int end_x   = min(cell_index_x + 2, (int)S.x - 1);
    int start_y = max(0, cell_index_y - 2);
    int end_y   = min(cell_index_y + 2, (int)S.y - 1);
    int start_z = max(0, cell_index_z - 2);
    int end_z   = min(cell_index_z + 2, (int)S.z - 1);

    for (int z = start_z; z <= end_z; z += 1) {
      for (int y = start_y; y <= end_y; y += 1) {
        for (int x = start_x; x <= end_x; x += 1) {
          int  index_of_Block = g[x][y][z] - 1;
          bool empty_block    = index_of_Block == -1;
          if (!empty_block) {
            vec3f dist_ = (p[index_of_Block] - pos);
            if (dist_.x * dist_.x + dist_.y * dist_.y + dist_.z * dist_.z <=
                r) {
              return false;
            }
          }
        }
      }
    }

    return true;
  }
  return false;
}
bool isValid(vec2f candidate, vec2f sampleRegionSize, float cellSize,
    float radius, vector<vec2f> points, vector<vector<int>> grid) {
  if (candidate.x >= 0 && candidate.x < sampleRegionSize.x &&
      candidate.y >= 0 && candidate.y < sampleRegionSize.y) {
    int cellX = (int)(candidate.x / cellSize);
    int cellY = (int)(candidate.y / cellSize);
    int searchStartX = max(0, cellX - 2);
    int searchEndX   = min((int)cellX + 2,(int) grid.at(0).size()-1);
    int searchStartY = max(0, cellY - 2);
    int searchEndY   = min((int)cellY + 2,(int) grid.at(1).size() - 1);

    for (int x = searchStartX; x <= searchEndX; x++) {
      for (int y = searchStartY; y < searchEndY; y++) {
        int pointIndex = grid.at(x).at(y) - 1;  
        if (pointIndex != -1) {
          vec2f a       = candidate - points.at(pointIndex);
          float sqrtDst = (a.x * a.x + a.y * a.y);
          if (sqrtDst < radius * radius) return false;
        }
      }    
    }
    return true;
  }
  return false;
}

vector<vec2f> poissonSampling(float radius, vec2f regionSize, int max_tries) {
  float cellSize = radius / sqrt(2.f);
  vector<vector<int>> grid;
  auto rng = make_rng(110101);


  for (int i = 0; i < ceil((regionSize.x) / cellSize); i++) {
    vector<int>  col;
    grid.push_back(col);
      for (int j = 0; j < ceil((regionSize.y) / cellSize) ;j++) {
      grid.at(i).push_back(-1);

    }
  }
  vector<vec2f> points;
  vector<vec2f> spawnPoints;
  
  spawnPoints.push_back(regionSize / 2);
  while (spawnPoints.size() > 0) {
    int spawnIndex = rand1i(rng , spawnPoints.size());
    vec2f spawnCentre = spawnPoints.at(spawnIndex);
    bool  candidateAccepted = false;
    for (int i = 0; i < max_tries; i++) {
      float angle = rand1f(rng) * 3.1415 * 2;
      auto  dir   = vec2f{sin(angle) , cos(angle)};
      vec2f candidate = spawnCentre + dir * (rand1i(rng, radius)+radius);
      if (isValid(candidate, regionSize, cellSize, radius, points, grid)) {
        points.push_back(candidate);
        spawnPoints.push_back(candidate);
        grid.at((int)candidate.x/cellSize).at((int)candidate.y/cellSize) = points.size();
        candidateAccepted = true;
        break;
      }
      if (!candidateAccepted) spawnPoints.erase(spawnPoints.begin() + spawnIndex);



    }
  }
  return points;

  }

    








vector<vec3f> poissonSampling3d(float radius, vec3f S, int number_of_possibilities) {
  float s = radius / sqrt(3.f);

  int*** Grid = new int**[ceil(S.x / s)];

    

  for (int i = 0; i < ceil(S.x / s); i++) {
    Grid[i] = new int*[ceil(S.y / s)];
    for (int j = 0; j < ceil(S.y / s); j++) {
      Grid[i][j] = new int[S.z / s];
    }
  }
  auto          rng = make_rng(110101);
  vector<vec3f> p;
  vector<vec3f> s_p;

  s_p.push_back(S / 2);

  while (s_p.size() > 0) {
    int   random_index   = rand() % s_p.size();
    vec3f s_c            = s_p.at(random_index);
    bool  point_Accepted = false;

    for (int i = 0; i < number_of_possibilities; i++) {
      vec3f new_p = s_c + rand3f(rng);
      if (Isvalid3d(new_p, S, s, radius, Grid, p)) {
        p.push_back(new_p);
        s_p.push_back(new_p);

        Grid[(int)(new_p.x / s)][(int)(new_p.y / s)][(int)(new_p.z / s)] =
            p.size();
        point_Accepted = true;
      }
    }
    if (!point_Accepted) {
      s_p.erase(s_p.begin() + random_index);
    }
  }

  return p;
}


float cell(const vec3f& p) {
  // distance
  vec3f x      = {floor(p.x), floor(p.y), floor(p.z)};
  vec3f f      = {p.x - floor(p.x), p.y - floor(p.y), p.z - floor(p.z)};
  auto  random = make_rng(172784);

  vec3f mb;
  vec3f mr;

  float res = 8.0;
  for (int k = -1; k<= 1 ; k++){
    for (int j = -1; j <= 1; j++) {
    for (int i = -1; i <= 1; i++) {
        vec3f b = {i, j, k};
        vec3f r = b + rand3f(random) - f;  //
        float d = dot(r, r);
        if (d < res) {
          res = d;
          mr  = r;
          mb  = b;
        }
      }
      }
    }
  res = 8.0;
  for (float k = -2; k <= 2; k++) {
    for (float j = -2; j <= 2; j++) {
      for (float i = -2; i <= 2; i++) {
        vec3f b = mb + vec3f{i, j, k};
        vec3f r = b + rand3f(random) - f;  //
        float d = dot(0.5 * (mr + r), normalize(r - mr));

        res = min(res, d);
      }
    }
  }
  return res;
}
  

 float getBorder(vec3f p) {
  auto d = cell(p);
  return 1.0 - smoothstep(0.f, 0.05f, d); }
    


 float voronoise(vec3f x, float u, float v) {
     auto rng = make_rng(198767);
     vec3f p   = {floor(x.x), floor(x.y), floor(x.z)};
     vec3f f   = {x.x - floor(x.x), x.y - floor(x.y), x.z - floor(x.z)};

     float k  = 1.0 + 63.0 * pow(1.0 - v, 4.0);
     float va = 0.0;
     float wt = 0.0;
     for (int j=-2;j<=2;j++)
       for (int i = -2; i <= 2; i++) {
         for (int kk = -2; kk <= 2; kk++) {
           vec3f g = {float(i), float(j), float(kk)};
           vec3f m = p + g;
           vec3f o = hash3(m) * vec3f{u, u, 1.0};
           vec3f r = g - f + vec3f{o.x, o.y, o.z};
           float d = dot(r, r);
           float w = pow(1.0f - smoothstep(0.0f, 1.414f, sqrt(d)), k);
           va += w * o.z;
           wt += w;
         }
       }
     return va / wt;
     
 }

 float smoothVoronoi(const vec3f& x) {  
   vec3f p = {floor(x.x), floor(x.y), floor(x.z)};
   vec3f f = {x.x - floor(x.x), x.y - floor(x.y), x.z - floor(x.z)};

   float res = 0.0;
   for (float k = -1; k <= 1; k++) {
   for (float j = -1; j <= 1; j++) {    
       for (float i = -1; i <= 1; i++) {
         vec3f b = {i, j, k};
         vec3f r = b - f + hash3(p + b);
         float d = dot(r, r);
         res += 1.0 / pow(d, 8.0f);

        
       }
     }
   }
   return pow(1.0 / res, 1.0 / 16.0);
 }

void add_polyline(shape_data& shape, const vector<vec3f>& positions,
    const vector<vec4f>& colors, float thickness = 0.0001f) {
  auto offset = (int)shape.positions.size();
  shape.positions.insert(
      shape.positions.end(), positions.begin(), positions.end());
  shape.colors.insert(shape.colors.end(), colors.begin(), colors.end());
  shape.radius.insert(shape.radius.end(), positions.size(), thickness);
  for (auto idx = 0; idx < positions.size() - 1; idx++) {
    shape.lines.push_back({offset + idx, offset + idx + 1});
  }
}

void sample_shape(vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, const shape_data& shape, int num) {
  auto triangles  = shape.triangles;
  auto qtriangles = quads_to_triangles(shape.quads);
  triangles.insert(triangles.end(), qtriangles.begin(), qtriangles.end());
  auto cdf = sample_triangles_cdf(triangles, shape.positions);
  auto rng = make_rng(19873991);
  for (auto idx = 0; idx < num; idx++) {
    auto [elem, uv] = sample_triangles(cdf, rand1f(rng), rand2f(rng));
    auto q          = triangles[elem];
    positions.push_back(interpolate_triangle(
        shape.positions[q.x], shape.positions[q.y], shape.positions[q.z], uv));
    normals.push_back(normalize(interpolate_triangle(
        shape.normals[q.x], shape.normals[q.y], shape.normals[q.z], uv)));
    if (!texcoords.empty()) {
      texcoords.push_back(interpolate_triangle(shape.texcoords[q.x],
          shape.texcoords[q.y], shape.texcoords[q.z], uv));
    } else {
      texcoords.push_back(uv);
    }
  }
}

void make_terrain(shape_data& shape, const terrain_params& params) {
  for (int i = 0; i < shape.positions.size(); i++) {
    vec3f position = shape.positions[i];
    position += shape.normals[i] *
                ridge(position * params.scale, params.octaves) * params.height *
                 (1 - length(position - params.center) / params.size);
    float perc = position.y / params.height;
    if (perc <= 0.3)
      shape.colors.push_back(params.bottom);
    else if (perc > 0.6)
      shape.colors.push_back(params.top);
    else
      shape.colors.push_back(params.middle);
    shape.positions[i] = position;
  }
  compute_normals(shape.normals, shape);
}

void make_displacement(shape_data& shape, const displacement_params& params) {
  
  for (int i = 0; i < shape.positions.size(); i++) {
    vec3f position = shape.positions[i];
    vec3f last_position = vec3f(position);
    position += shape.normals[i] *
                turbulence(position * params.scale, params.octaves) *
                params.height;
    shape.colors.push_back(interpolate_line(params.bottom, params.top,
        distance(position, last_position) / params.height));
    shape.positions[i] = position;
  }
  compute_normals(shape.normals, shape);
}

void make_cell(shape_data& shape, const displacement_params& params) {
   for (int i = 0; i < shape.positions.size(); i++) {
    auto oldpos = shape.positions[i];
    float b      = getBorder(oldpos * params.scale);
    shape.positions[i] += shape.normals[i] * (b * params.height);
    shape.colors.push_back(interpolate_line(params.bottom, params.top,
        distance(oldpos, shape.positions[i]) / params.height));
  }
  compute_normals(shape.normals, shape);



}
void make_voronoise(shape_data& shape, const displacement_params& params) {
  for (int i = 0; i < shape.positions.size(); i++) {
    vec3f position      = shape.positions[i];
    vec3f last_position = vec3f(position);
    position += shape.normals[i] *
                voronoise(position * params.scale, 1,1) *
                params.height;
    shape.colors.push_back(interpolate_line(params.bottom, params.top,
        distance(position, last_position) / params.height));
    shape.positions[i] = position;
  }
  compute_normals(shape.normals, shape);
}

void make_smoothvoronoi(shape_data& shape, const displacement_params& params) {
  for (int i = 0; i < shape.positions.size(); i++) {
    vec3f position      = shape.positions[i];
    vec3f last_position = vec3f(position);
    position += shape.normals[i] * smoothVoronoi(position * params.scale) * params.height;

    shape.colors.push_back(interpolate_line(params.bottom, params.top,
        distance(position, last_position) / params.height));

    shape.positions[i] = position;
  }
  compute_normals(shape.normals, shape);
}

void make_multinoises(shape_data& shape, const displacement_params& params) {
  
    displacement_params params2 = {};
  params2.bottom = srgb_to_rgb(vec4f{64, 200, 224, 255} / 255);
  params2.top =  srgb_to_rgb(vec4f{255, 255, 255, 255} / 255);

  displacement_params params3 = {};
  params3.bottom              = srgb_to_rgb(vec4f{255, 0, 0, 255} / 255);
  params3.top                  = srgb_to_rgb(vec4f{0, 0, 255, 255} / 255);



   auto rng = make_rng(110101);
  int  a   = rand1f(rng) * 3;
  for (int i = 0; i < shape.positions.size(); i++) {
    vec3f position = shape.positions[i];
    vec3f last_position = vec3f(position);
    float posy = position.y;
    float posx = position.x;
    int   region = -1;
    if ((posy > 0.10 && posy < 0.15 && posx > -0.075 && posx < -0.025) ||
        (posy > 0.05 && posy < 0.10 && posx > -0.025 && posx < 0.025) ||
        (posy > 0 && posy < 0.05 && posx > 0.025 && posx < 0.075)) {
      position += shape.normals[i] *
                  turbulence(position * params.scale, params.octaves) *
                  params.height;
      region = 1;
    }

    if ((posy > 0.10 && posy < 0.15 && posx > -0.025 && posx < 0.025) ||
        (posy > 0.05 && posy < 0.10 && posx > 0.025 && posx < 0.075) ||
        (posy > 0 && posy < 0.05 && posx > -0.075 && posx < -0.025)) {
      position += shape.normals[i] * smoothVoronoi(position * params.scale) *
                  params.height;
      region = 2;
    }
    if ((posy > 0.10 && posy < 0.15 && posx > 0.025 && posx < 0.075) ||
        (posy > 0.05 && posy < 0.10 && posx > -0.075 && posx < -0.025) ||
        (posy > 0 && posy < 0.05 && posx > -0.025 && posx < 0.025)) {
      position += shape.normals[i] * cell(position * params.scale) *
                  params.height;
      region = 3;
    }
    
  if (region == 1) {
      shape.colors.push_back(interpolate_line(params.bottom, params.top,
          distance(position, last_position) / params.height));
  }
  if (region == 2) {
    shape.colors.push_back(interpolate_line(params2.bottom, params2.top,
        distance(position, last_position) / params.height));
  }
    
   
   if (region == 3) {
      shape.colors.push_back(interpolate_line(params3.bottom, params3.top,
          distance(position, last_position) / params.height));
    }

    shape.positions[i] = position;
  }
 
  compute_normals(shape.normals, shape);
}

void make_hair(
    shape_data& hair, const shape_data& shape, const hair_params& params) {
  auto new_shape = shape;
  auto hr = params.lenght / params.steps;
  auto random    = make_rng(110101);
 
  sample_shape(new_shape.positions , new_shape.normals , new_shape.texcoords , new_shape , params.num);
  for (int i = 0; i < new_shape.positions.size(); i++) {
    auto r = rand1f(random);
    if (params.density > r) {
      vector<vec3f> positions;
      vector<vec4f> colors;
      positions.push_back(new_shape.positions[i]);
      colors.push_back(params.bottom);
      vec3f normal = new_shape.normals[i];
      for (int j = 0; j < params.steps; j++) {
        vec3f next = positions[j] + hr * normal +
                     noise3(positions[j] * params.scale) * params.strength;
        next.y -= params.gravity;
        normal = normalize(next - positions[j]);
        positions.push_back(next);
        colors.push_back(interpolate_line(params.bottom, params.top,
            distance(next, positions[0]) / params.lenght));
      }
      colors[params.steps] = params.top;
      add_polyline(hair, positions, colors);
    }
  }
      vector<vec3f> tangents = lines_tangents(hair.lines, hair.positions);
      for (int i = 0; i < tangents.size(); i++) {
        vec4f tangente = vec4f{tangents[i].x, tangents[i].y, tangents[i].z, 0};
        hair.tangents.push_back(tangente);
      }
  
  
}
void make_grass(scene_data& scene, const instance_data& object,
    const vector<instance_data>& grasses, const grass_params& params) {
  auto random = make_rng(110101);
  instance_data new_object = object;
  auto new_shape = scene.shapes[object.shape];
  auto size = new_shape.positions.size();
  sample_shape(new_shape.positions, new_shape.normals ,new_shape.texcoords, new_shape , params.num);
  for (int i = 0; i < new_shape.positions.size(); i++) {
    auto          r        = rand1f(random);
    if (params.density > r) {
      auto          position = new_shape.positions[i];
      auto          grass    = grasses[rand1i(random, grasses.size())];
      instance_data obj;

      grass.frame.y = new_shape.normals[i];
      grass.frame.x = normalize(
          vec3f{1, 0, 0} - dot(vec3f{1, 0, 0}, grass.frame.y) * grass.frame.y);
      grass.frame.z = cross(grass.frame.x, grass.frame.y);
      grass.frame.o = position;

      float random_scaling = 0.9f + rand1f(random) * 0.1f;
      grass.frame *= scaling_frame(
          vec3f{random_scaling, random_scaling, random_scaling});
      float random_rotation_1 = rand1f(random) * 2 * pif;
      grass.frame *= rotation_frame(grass.frame.y, random_rotation_1);
      float random_rotation_2 = 0.1f + rand1f(random) * 0.1f;
      grass.frame *= rotation_frame(grass.frame.z, random_rotation_2);
      scene.instances.push_back(grass);
    }

  }

}
 
 


 

}  // namespace yocto
