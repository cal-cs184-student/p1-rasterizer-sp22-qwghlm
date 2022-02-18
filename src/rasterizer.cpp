#include "rasterizer.h"
#include <cmath>
using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)

    for (int w = 0; w < sample_rate; w++) {
        sample_buffer[sample_rate * (y * width + x) + w] = c;
    }

  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    // TODO: Task 2: Update to implement super-sampled rasterization
    int max_x = (int) max(max(x0, x1), x2);
    int max_y = (int) max(max(y0, y1), y2);
    int min_x = (int) min(min(x0, x1), x2);
    int min_y = (int) min(min(y0, y1), y2);

    auto samples = (float) (sqrt(sample_rate));
    float x, y;
    for (int i = min_x; i <= max_x; i++) {
        for (int j = min_y; j <= max_y; j++) {
            int count = 0;
            for (int w = 0; w < samples; w++) {
                for (int z = 0; z  < samples; z++) {

                    x = (2*((float) i) + ((float) w)/samples + ((float) w + 1)/samples)/2;
                    y = (2*((float) j) + ((float) z)/samples + ((float) z + 1)/samples)/2;

                    if (insideTri(x0, y0, x1, y1, x2, y2, x, y)) {
                        sample_buffer[sample_rate * (j * width + i) + count] = color;
                    }
                    count++;
                }
            }
        }
    }
  }

  bool RasterizerImp::insideTri(float x0, float y0,
                                float x1, float y1,
                                float x2, float y2,
                                float i, float j) {
    float first_plane = -(i - x0)*(y0-y2) + (j - y0)*(x0-x2);
    float second_plane = -(i - x1)*(y1-y0) + (j - y1)*(x1-x0);
    float third_plane = -(i - x2)*(y2-y1) + (j - y2)*(x2-x1);

    float fourth_plane = -(i - x2)*(y1-y2) + (j - y2)*(x1-x2);
    float fifth_plane = -(i - x1)*(y0-y1) + (j - y1)*(x0-x1);
    float sixth_plane = -(i - x0)*(y2-y0) + (j - y0)*(x2-x0);


    return (first_plane >= 0 && second_plane >= 0 && third_plane >= 0) ||
            (fourth_plane >=0 && fifth_plane >= 0 && sixth_plane >= 0);

  }

  bool RasterizerImp::inside(float tri, float x0, float y0,
              float x1, float y1,
              float x2, float y2,
              float  i, float j) {
    float A1 = area(i, j, x1, y1, x2, y2);
    float A2 = area(x0, y0, i, j, x2, y2);
    float A3 = area(x0, y0, x1, y1, i, j);

    if ((int) tri == (int) (A1 + A2 + A3)) {
        return true;
    }
    return false;
  }

  float RasterizerImp::area(float x1, float y1,
             float x2, float y2,
             float x3, float y3) {

      return (float) std::abs ((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0);
  }

  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
    int max_x = (int) max(max(x0, x1), x2);
    int max_y = (int) max(max(y0, y1), y2);
    int min_x = (int) min(min(x0, x1), x2);
    int min_y = (int) min(min(y0, y1), y2);

    auto samples = (float) (sqrt(sample_rate));
    float alpha, beta, gamma;
    float x, y;

    for (int i = min_x; i <= max_x; i++) {
      for (int j = min_y; j <= max_y; j++) {
          int count = 0;
          for (int w = 0; w < samples; w++) {
              for (int z = 0; z  < samples; z++) {

                  x = (2*((float) i) + ((float) w)/samples + ((float) w + 1)/samples)/2;
                  y = (2*((float) j) + ((float) z)/samples + ((float) z + 1)/samples)/2;
                  float A_b = area(x0, y0, x, y, x2, y2);
                  float A_a = area(x2, y2, x, y, x1, y1);
                  float A_c = area(x0, y0, x, y, x1, y1);
                  float A = A_a + A_b + A_c;
                  alpha = A_a / A;
                  beta = A_b / A;
                  gamma = A_c / A;
                  if (insideTri(x0, y0, x1, y1, x2, y2, x, y)) {
                      sample_buffer[sample_rate * (j * width + i) + count] = alpha * c0 + beta * c1 + gamma * c2;
                  }
                  count++;
              }
          }
      }
    }


  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
    int max_x = (int) max(max(x0, x1), x2);
    int max_y = (int) max(max(y0, y1), y2);
    int min_x = (int) min(min(x0, x1), x2);
    int min_y = (int) min(min(y0, y1), y2);

    auto samples = (float) (sqrt(sample_rate));
    float x, y;

    for (int i = min_x; i <= max_x; i++) {
        for (int j = min_y; j <= max_y; j++) {
            int count = 0;

            for (int w = 0; w < sample_rate; w++) {
                for (int z = 0; z < sample_rate; z++) {
                    x = (2*((float) i) + ((float) w)/samples + ((float) w + 1)/samples)/2;
                    y = (2*((float) j) + ((float) z)/samples + ((float) z + 1)/samples)/2;


                    if (insideTri(x0, y0, x1, y1, x2, y2, x, y)) {
                        Vector2D p_uv = compute_barycentric_coordinates(x0, y0, u0, v0, x1, y1, u1, v1, x2, y2, u2, v2, x, y);
                        Vector2D p_dx_uv = compute_barycentric_coordinates(x0, y0, u0, v0, x1, y1, u1, v1, x2, y2, u2, v2, x+1, y);
                        Vector2D p_dy_uv = compute_barycentric_coordinates(x0, y0, u0, v0, x1, y1, u1, v1, x2, y2, u2, v2, x, y+1);
                        SampleParams sp = SampleParams();
                        sp.lsm = lsm;
                        sp.psm = psm;
                        sp.p_uv = p_uv;
                        sp.p_dy_uv = p_dy_uv;
                        sp.p_dx_uv = p_dx_uv;
                        Color c = tex.sample(sp);

                        sample_buffer[sample_rate * (j * width + i) + count] = c;
                    }
                    count++;
                }
            }
        }
    }



  }

  Vector2D RasterizerImp::compute_barycentric_coordinates(float x0, float y0, float u0, float v0,
                                                          float x1, float y1, float u1, float v1,
                                                          float x2, float y2, float u2, float v2,
                                                          float x, float y) {
      float A_b = area(x0, y0, x, y, x2, y2);
      float A_a = area(x2, y2, x, y, x1, y1);
      float A_c = area(x0, y0, x, y, x1, y1);
      float A = A_a + A_b + A_c;
      float alpha = A_a / A;
      float beta = A_b / A;
      float gamma = A_c / A;

      float u = alpha * u0 + beta * u1 + gamma * u2;
      float v = alpha * v0 + beta * v1 + gamma * v2;

      return Vector2D(u, v);
  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;

    this->sample_buffer.resize(width * height * rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;

    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support
    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = sample_buffer[sample_rate * (y * width) + x];
        float avg_r = 0;
        float avg_g = 0;
        float avg_b = 0;

        for (int w = 0; w < sample_rate; w++) {
            avg_r += sample_buffer[(y * width + x) * sample_rate + w].r * 255;
            avg_g += sample_buffer[(y * width + x) * sample_rate + w].g * 255;
            avg_b += sample_buffer[(y * width + x) * sample_rate + w].b * 255;
        }

        avg_r /= (float) sample_rate;
        avg_g /= (float) sample_rate;
        avg_b /= (float) sample_rate;

        this->rgb_framebuffer_target[3 * (y * width + x)] = (unsigned char) avg_r;
        this->rgb_framebuffer_target[3 * (y * width + x) + 1] = (unsigned char) avg_g;
        this->rgb_framebuffer_target[3 * (y * width + x) + 2] = (unsigned char) avg_b;

      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
