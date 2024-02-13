#include "rasterizer.h"
#include <cmath>
#include "vector2D.h"

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


    sample_buffer[y * width + x] = c;
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

  vector<Vector2D> makeCounterClockwise(float x0, float y0,
    float x1, float y1,
    float x2, float y2) {

      Vector2D v01 = Vector2D(x0 - x1, y0 - y1);
      Vector2D v12 = Vector2D(x1 - x2, y1 - y2);
      Vector2D v20 = Vector2D(x2 - x0, y2 - y0);
      vector<Vector2D> edges;

      if (cross(v01, v20) > cross(v01, v12)) {
        edges = {Vector2D(x0, y0), 
                 Vector2D(x2, y2), 
                 Vector2D(x1, y1)};
      } else {
        edges = {Vector2D(x0, y0), 
                 Vector2D(x1, y1), 
                 Vector2D(x2, y2)};
      }

      return edges;
    }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {

    vector<Vector2D> coords = makeCounterClockwise(x0, y0, x1, y1, x2, y2);
    x0 = coords[0][0]; x1 = coords[1][0]; x2 = coords[2][0];
    y0 = coords[0][1]; y1 = coords[1][1]; y2 = coords[2][1];

    int xmin = floor(min(min(x0, x1), x2));
    int xmax = ceil(max(max(x0, x1), x2));
    int ymin = floor(min(min(y0, y1), y2));
    int ymax = ceil(max(max(y0, y1), y2));

    int dim = ceil(sqrt(sample_rate));
    double step = 1.0 / dim;

    for (double x=xmin; x<xmax; x++) {
      for (double y=ymin; y<ymax; y++) {
        
        double scale = 0;
        double i_min = x + (0.5*step);
        double j_min = y + (0.5*step);
        double iters = 0;

        for (double i=i_min; i<i_min+1-(0.5*step); i=i+step) {
          for (double j=j_min; j<j_min+1-(0.5*step); j=j+step) {
            double l0 = -(i - x0)*(y1-y0) + (j - y0)*(x1-x0);
            double l1 = -(i - x1)*(y2-y1) + (j - y1)*(x2-x1);
            double l2 = -(i - x2)*(y0-y2) + (j - y2)*(x0-x2);
            if (l0 >= 0 && l1 >= 0 && l2 >= 0) {
              scale += 1;
            }
            iters += 1;
          }
        }

        scale = scale / iters;
        if (scale > 0) {
          fill_pixel(x, y, color * scale);
        }
      }
    }
  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {

    vector<Vector2D> coords = makeCounterClockwise(x0, y0, x1, y1, x2, y2);
    x0 = coords[0][0]; x1 = coords[1][0]; x2 = coords[2][0];
    y0 = coords[0][1]; y1 = coords[1][1]; y2 = coords[2][1];

    int xmin = floor(min(min(x0, x1), x2));
    int xmax = ceil(max(max(x0, x1), x2));
    int ymin = floor(min(min(y0, y1), y2));
    int ymax = ceil(max(max(y0, y1), y2));

    double L0 = -(x0-x1)*(y2-y1) + (y0-y1)*(x2-x1);
    double L1 = -(x1-x2)*(y0-y2) + (y1-y2)*(x0-x2);

    for (int x = xmin; x < xmax; x++) {
      for (int y = ymin; y < ymax; y++) {
        
        double i = x + 0.5;
        double j = y + 0.5;
        double l0 = -(i - x0)*(y1-y0) + (j - y0)*(x1-x0);
        double l1 = -(i - x1)*(y2-y1) + (j - y1)*(x2-x1);
        double l2 = -(i - x2)*(y0-y2) + (j - y2)*(x0-x2);
        if (l0 >= 0 && l1 >= 0 && l2 >= 0) {
          double L12 = -(x-x1)*(y2-y1) + (y-y1)*(x2-x1);
          double L20 = -(x-x2)*(y0-y2) + (y-y2)*(x0-x2);
          double a = L12 / L0;
          double b = L20 / L1;
          double c = 1 - a - b;
          Color color = (a*c0) + (b*c1) + (c*c2);
          fill_pixel(x, y, color);
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




  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height, Color::White);
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
        Color col = sample_buffer[y * width + x];

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
